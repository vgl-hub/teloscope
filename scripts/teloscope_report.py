#!/usr/bin/env python3
"""
teloscope_report.py — Publication-ready figures from Teloscope output.

Usage:
    python teloscope_report.py <output_directory> [-o report.pdf]
    python teloscope_report.py <output_directory> --png

Reads Teloscope output files and writes separate terminal and ITS PDFs by default.
  Terminal: assembly overview and per-chromosome terminal zoom figures.
  ITS: observed-row distributions, paginated atlas, candidates, selected loci.
Use --section all for a combined PDF, or terminal/its to select one section.

Requires: Python 3.6+, matplotlib, numpy, pandas
"""

import sys
import os
import glob
import re
import argparse
import textwrap
import inspect
from collections import Counter, defaultdict, OrderedDict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch, ConnectionPatch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection, PatchCollection, PolyCollection
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.patheffects as patheffects
from matplotlib import cbook, transforms
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath

# ---------------------------------------------------------------------------
# Nature-style configuration
# ---------------------------------------------------------------------------

# Nature figure widths: 89 mm (single column), 183 mm (double column)
FIG_WIDTH_SINGLE = 3.50    # inches (89 mm)
FIG_WIDTH_DOUBLE = 7.20    # inches (183 mm)
REPORT_PAGE_HEIGHT = 2.2 + 0.25 * 6  # terminal page with all six tracks
SLIDE_SIZE = (FIG_WIDTH_DOUBLE, REPORT_PAGE_HEIGHT)  # one page size for every report page
ATLAS_PITCH_IN, ATLAS_GROUP_GAP_IN, ATLAS_BAR_IN = 0.10, 0.05, 0.04  # fixed row pitch, extra gap between homolog groups, bar
ATLAS_HEADER_IN, ATLAS_KEY_IN, ATLAS_FOOT_IN, ATLAS_AXIS_IN, ATLAS_TITLE_IN = 0.36, 0.74, 0.06, 0.30, 0.20
ATLAS_COLUMNS, ATLAS_COL_GAP_IN = 2, 0.30
ATLAS_TIER_RATIO = 4.0  # size range per atlas panel, so every scaffold spans at least a quarter of its axis
ATLAS_HEIGHTS = (SLIDE_SIZE[1], 2 * SLIDE_SIZE[1])  # one slide, else a double slide within Nature's 247 mm page
ITS_LENGTH_COLOR = "#91C9C1"  # teal, outside the double-key hues
ATLAS_BACKBONE = "#E8E8E8"  # light scaffold backbone behind atlas ITS cells

COLORS = {
    # Classification palette (colorblind-friendly, quality-graduated)
    "T2T":                 "#08519C",
    "Gapped T2T":          "#9ECAE1",
    "Incomplete":          "#66A61E",
    "Gapped Incomplete":   "#B3D38F",
    "Fragmented":          "#E6AB02",
    "Gapped Fragmented":   "#F2D580",
    "Discordant":          "#D7191C",
    "Gapped Discordant":   "#EB8C8D",
    "No telomeres":        "#762A83",
    "Gapped No telomeres": "#BB95C1",
    # Arm / track colors
    "p":            "#0072B2",
    "q":            "#D55E00",
    "b":            "#FFC754",
    "density":      "#2A9D59",
    "canonical":    "#1B9ECA",
    "strand_bias":  "#8D80C6",
    "gc":           "#E69F00",
    "entropy":      "#CC79A7",
    "terminal":     "#4A4A4A",    # dark matte grey
    "its":          "#8C8C8C",    # medium grey
    "gap":          "#D6D6D6",    # light grey
    "no_data":      "#e0e0e0",
}

# ITS junction classes, fixed order.
CLASS_ORDER = ["fusion", "tail_to_tail", "fragmentation", "single"]

# Text glyphs keep the directional symbols light while remaining editable in PDF export.
BLOCK_GLYPHS = {
    "p": "<",
    "q": ">",
    "b": "<>",
}

# Orientation palette shared by terminal telomeres, ITS report pages and plot_its.py (moved here so plot_its can import it without a circular import).
ITS_ORIENT_COLORS = {k: COLORS[k] for k in ("p", "q", "b")}
ITS_ORIENT_COLORS["unknown"] = COLORS["its"]
ITS_ORIENT_LABELS = (("p", "p (forward)"), ("q", "q (reverse)"),
                     ("b", "balanced"), ("unknown", "unknown"))
ZOOM_COLOR = COLORS["Discordant"]  # red box/funnel marking a zoom region
# Double key: hue = strand (rev, both, fwd), opacity = canonicity (nonCan, mixed, can).
DOUBLE_KEY_BASE = ("#D55E00", "#88307A", "#0072B2")
DOUBLE_KEY_ALPHA = (0.40, 0.70, 1.0)
DOUBLE_KEY_RGBA = tuple(tuple(mcolors.to_rgba(b, a) for b in DOUBLE_KEY_BASE) for a in DOUBLE_KEY_ALPHA)
DOUBLE_KEY_COLORS = tuple(tuple(mcolors.to_hex([a * v + 1 - a for v in mcolors.to_rgb(b)]) for b in DOUBLE_KEY_BASE)
                          for a in DOUBLE_KEY_ALPHA)  # the same cells flattened on white, for opaque marks
DOUBLE_KEY_STRAND = ("rev", "both", "fwd")
DOUBLE_KEY_CANON = ("nonCan", "mixed", "can")
LABEL_THRESHOLD = 0.667  # engine default for symmetric thirds
ITS_CLUSTER_MERGE_GAP = 50_000
MAX_COORD = 2**53 - 1  # exact-integer limit for float64 BED coordinates
ITS_CLUSTER_MIN_ROWS = 3

OVERVIEW_DASH_STYLE = (0, (2.2, 2.2))
OVERVIEW_DASH_WIDTH = 0.35
FLAGGED_SCAFFOLD_CATEGORIES = (
    "Fragmented",
    "Gapped Fragmented",
    "Discordant",
    "Gapped Discordant",
)

FIGURE_SUMMARY_SIZE = 6.0
PANEL_LABEL_SIZE = 8.0
PANEL_TITLE_SIZE = 6.9
AXIS_LABEL_SIZE = 6.2
AXIS_TICK_SIZE = 5.5
LEGEND_TEXT_SIZE = 5.4
MIN_TEXT_SIZE = 5.0
ANNOTATION_TEXT_SIZE = MIN_TEXT_SIZE
PLACEHOLDER_TEXT_SIZE = 6.2
BLOCK_GLYPH_SIZE = MIN_TEXT_SIZE
BLOCK_GLYPH_MIN_FRACTION = {
    "p": 0.018,
    "q": 0.018,
    "b": 0.030,
}
TRACK_LABEL_X = -0.10
PANEL_LETTER_X_IN = 0.24  # leftmost panel letter, as on the terminal page
PANEL_LETTER_DX_IN = 0.36  # other letters sit this far left of their axes

def _apply_nature_style():
    """Apply Nature journal rcParams globally."""
    plt.rcParams.update({
        # Fonts — Nature requires sans-serif (Arial / Helvetica)
        "font.family":        "sans-serif",
        "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":          7,
        # Axes
        "axes.titlesize":     PANEL_TITLE_SIZE,
        "axes.labelsize":     AXIS_LABEL_SIZE,
        "axes.linewidth":     0.5,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        # Ticks
        "xtick.labelsize":    AXIS_TICK_SIZE,
        "ytick.labelsize":    AXIS_TICK_SIZE,
        "xtick.major.width":  0.5,
        "ytick.major.width":  0.5,
        "xtick.major.size":   3,
        "ytick.major.size":   3,
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        # Lines
        "lines.linewidth":    0.8,
        "lines.markersize":   3,
        # Legend
        "legend.fontsize":    LEGEND_TEXT_SIZE,
        "legend.frameon":     False,
        # Figure
        "figure.dpi":         150,
        "savefig.dpi":        450,
        "savefig.bbox":       None,
        "savefig.pad_inches": 0.02,
        "savefig.transparent": False,
        # PDF — TrueType embedding (required by Nature)
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
    })

_apply_nature_style()


def _warn(message):
    """Emit a warning message to stderr."""
    print(f"Warning: {message}", file=sys.stderr)


def _bed_header(line):
    """Recognize directives without dropping scaffold IDs such as track_001."""
    fields = line.split()
    return (not fields or line.startswith("#") or
            fields[0] in ("track", "browser") and
            (len(fields) < 2 or not fields[1].lstrip("-").isdigit()))


def _bad_interval(start, end):
    """True when a BED start/end pair is negative, empty, or exceeds exact float64 precision."""
    return start < 0 or end <= start or end > MAX_COORD

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def find_files(directory):
    """Auto-detect Teloscope output files in *directory*."""
    files = {}
    patterns = {
        "terminal":        "*_terminal_telomeres.bed",
        "interstitial":    "*_interstitial_telomeres.bed",
        "gaps":            "*_gaps.bed",
        "density":         "*_window_repeat_density.bedgraph",
        "canonical_ratio": "*_window_canonical_ratio.bedgraph",
        "strand_ratio":    "*_window_strand_ratio.bedgraph",
        "gc":              "*_window_gc.bedgraph",
        "entropy":         "*_window_entropy.bedgraph",
        "report":          "*_report.tsv",
    }
    for key, pat in patterns.items():
        hits = sorted(glob.glob(os.path.join(directory, pat)))
        if hits:
            if len(hits) > 1:
                _warn(f"Multiple matches for '{pat}' in '{directory}'; using '{hits[0]}'.")
            files[key] = hits[0]
    return files


def parse_terminal_bed(path):
    """
    Parse a terminal or interstitial telomere BED -> dict[chrom -> list of block dicts].
    12-column schema: chr start end teloLen teloLabel closestEnd fwdCan revCan
    fwdNonCan revNonCan chrSize teloType.
    """
    blocks = defaultdict(list)
    malformed = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if _bed_header(line):
                continue
            parts = line.split("\t")
            if len(parts) != 12:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: expected 12 BED columns, found {len(parts)}; skipping.")
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                telo_len = int(parts[3])
                label = parts[4]
                closest_end = parts[5]
                fwd_can = int(parts[6])
                rev_can = int(parts[7])
                fwd_noncan = int(parts[8])
                rev_noncan = int(parts[9])
                path_size = int(parts[10]) if parts[10] else 0
                terminality = parts[11]
            except ValueError as exc:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid numeric field ({exc}); skipping.")
                continue

            numeric = (start, end, telo_len, fwd_can, rev_can, fwd_noncan, rev_noncan, path_size)
            if (_bad_interval(start, end) or telo_len <= 0 or
                    min(fwd_can, rev_can, fwd_noncan, rev_noncan, path_size) < 0 or
                    max(numeric) > MAX_COORD):
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid interval start={start} end={end} teloLen={telo_len}; skipping.")
                continue

            chrom = parts[0]
            blocks[chrom].append({
                "start":     start,
                "end":       end,
                "length":    telo_len,
                "label":     label,
                "closestEnd": closest_end,
                "fwd":       fwd_can + fwd_noncan,
                "rev":       rev_can + rev_noncan,
                "can":       fwd_can + rev_can,
                "noncan":    fwd_noncan + rev_noncan,
                "fwdCan":    fwd_can,
                "revCan":    rev_can,
                "fwdNonCan": fwd_noncan,
                "revNonCan": rev_noncan,
                "pathSize":  path_size,
                "term":      terminality,
            })
    if malformed:
        suffix = " (first 3 shown above)" if malformed > 3 else ""
        _warn(f"Skipped {malformed} malformed BED line(s) from '{path}'{suffix}.")
    return dict(blocks)


def parse_interval_bed(path):
    """Parse a simple BED file -> dict[chrom -> list of interval dicts]."""
    intervals = defaultdict(list)
    malformed = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if _bed_header(line):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: expected at least 3 BED columns, found {len(parts)}; skipping.")
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as exc:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid BED coordinate ({exc}); skipping.")
                continue
            if _bad_interval(start, end):
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid interval start={start} end={end}; skipping.")
                continue
            intervals[parts[0]].append({
                "start": start,
                "end": end,
                "label": parts[3] if len(parts) > 3 else "",
            })
    if malformed:
        suffix = " (first 3 shown above)" if malformed > 3 else ""
        _warn(f"Skipped {malformed} malformed BED interval line(s) from '{path}'{suffix}.")
    return dict(intervals)


def _parse_bedgraph_fallback(path):
    """Robust line-by-line BEDgraph parser used if pandas parsing fails."""
    per_chrom = OrderedDict()
    malformed = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if _bed_header(line):
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: expected 4 BEDgraph columns, found {len(parts)}; skipping.")
                continue

            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                value = float(parts[3])
            except ValueError as exc:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid BEDgraph value ({exc}); skipping.")
                continue

            if _bad_interval(start, end) or not np.isfinite(value):
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid interval/value start={start} end={end} value={value}; skipping.")
                continue

            if chrom not in per_chrom:
                per_chrom[chrom] = [[], [], []]
            per_chrom[chrom][0].append(start)
            per_chrom[chrom][1].append(end)
            per_chrom[chrom][2].append(value)

    if malformed:
        suffix = " (first 3 shown above)" if malformed > 3 else ""
        _warn(f"Skipped {malformed} malformed BEDgraph line(s) from '{path}'{suffix}.")

    data = OrderedDict()
    for chrom, (starts, ends, values) in per_chrom.items():
        order = np.argsort(starts, kind="stable")
        data[chrom] = tuple(a[order] for a in (
            np.asarray(starts, dtype=np.int64),
            np.asarray(ends, dtype=np.int64),
            np.asarray(values, dtype=np.float64),
        ))
    return data


def parse_bedgraph(path):
    """
    Parse a BEDgraph file -> OrderedDict[chrom -> (starts, ends, values)].

    Returns numpy arrays per chromosome for fast downstream processing.
    """
    # Count header lines to skip (track/comment lines)
    skip = 0
    with open(path) as fh:
        for line in fh:
            if _bed_header(line):
                skip += 1
            else:
                break

    try:
        df = pd.read_csv(path, sep="\t", header=None, skiprows=skip,
                         names=["chrom", "start", "end", "value"],
                         dtype={"chrom": "category", "start": np.int64,
                                "end": np.int64, "value": np.float64},
                         engine="c", on_bad_lines="skip")
    except pd.errors.EmptyDataError:
        _warn(f"BEDgraph file '{path}' is empty; skipping.")
        return OrderedDict()
    except Exception as exc:
        _warn(f"Fast BEDgraph parse failed for '{path}' ({exc}); retrying line-by-line.")
        return _parse_bedgraph_fallback(path)

    valid = ((df["start"] >= 0) & (df["end"] > df["start"]) &
             (df["end"] <= MAX_COORD) & np.isfinite(df["value"])).to_numpy()
    if not valid.all():
        _warn(f"Skipped {int((~valid).sum())} invalid BEDgraph row(s) from '{path}'.")
    # Plain arrays, filtered only when needed; each scaffold becomes a slice (a view) of them.
    codes, names = df["chrom"].cat.codes.to_numpy(), df["chrom"].cat.categories
    starts, ends, values = (df[c].to_numpy() for c in ("start", "end", "value"))
    if not valid.all():
        codes, starts, ends, values = codes[valid], starts[valid], ends[valid], values[valid]
    del df
    cuts = np.flatnonzero(codes[1:] != codes[:-1]) + 1
    if len(cuts) + 1 > len(np.unique(codes)):  # a scaffold's windows are split up: group them, keeping file order
        order = np.argsort(codes, kind="stable")
        codes, starts, ends, values = codes[order], starts[order], ends[order], values[order]
        cuts = np.flatnonzero(codes[1:] != codes[:-1]) + 1
    data = {}
    for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, len(codes)]):
        s, e, v = starts[lo:hi], ends[lo:hi], values[lo:hi]
        if np.any(s[1:] < s[:-1]):  # engine output is already ordered; sort only when it is not
            order = np.argsort(s, kind="stable")
            s, e, v = s[order], e[order], v[order]
        data[str(names[codes[lo]])] = (s, e, v)
    return OrderedDict(sorted(data.items()))


_ITS_COLUMNS = ["chr", "start", "end", "teloLen", "teloLabel", "closestEnd",
                "fwdCan", "revCan", "fwdNonCan", "revNonCan", "chrSize", "teloType"]
_ITS_DTYPES = {"chr": str, "start": np.int64, "end": np.int64, "teloLen": np.int64,
               "teloLabel": str, "closestEnd": str, "fwdCan": np.int64, "revCan": np.int64,
               "fwdNonCan": np.int64, "revNonCan": np.int64, "chrSize": np.int64, "teloType": str}


def load_its_frame(path, motif_len=6, blocks=None):
    """Read the 12-column interstitial BED into a vectorised DataFrame with derived columns."""
    if not isinstance(motif_len, int) or motif_len <= 0:
        raise ValueError("motif_len must be a positive integer")
    parsed = parse_terminal_bed(path) if blocks is None else blocks
    rows = [(chrom, b["start"], b["end"], b["length"], b["label"], b["closestEnd"],
             b["fwdCan"], b["revCan"], b["fwdNonCan"], b["revNonCan"], b["pathSize"], b["term"])
            for chrom, blist in parsed.items() for b in blist]
    df = pd.DataFrame(rows, columns=_ITS_COLUMNS).astype(_ITS_DTYPES)

    df["canonical_bp"] = (df["fwdCan"] + df["revCan"]) * motif_len
    df["can_prop"] = np.where(df["teloLen"] > 0, df["canonical_bp"] / df["teloLen"], np.nan)
    # Missing/inconsistent scaffold sizes are not evidence of a q end; leave relative coordinates undefined instead.
    sizes = df.groupby("chr")["chrSize"].transform("max")
    ends = df.groupby("chr")["end"].transform("max")
    distinct = df["chrSize"].where(df["chrSize"] > 0).groupby(df["chr"]).transform("nunique")
    known_size = (sizes >= ends) & (sizes > 0) & (distinct == 1)
    uncertain = int(df.loc[~known_size, "chr"].nunique())
    if uncertain:
        _warn(f"{uncertain} ITS scaffold(s) have missing or inconsistent sizes; relative positions are undefined.")
    df["chrSize"] = sizes.where(known_size, 0).astype(np.int64)
    df["pos_frac"] = ((df["start"] + df["end"]) / 2.0 / df["chrSize"].replace(0, np.nan))
    df["end_dist"] = np.minimum(df["start"], df["chrSize"] - df["end"]).where(known_size)
    return df


def load_gaps_frame(path):
    """Read a BED3 gaps file into a plain chr/start/end DataFrame."""
    cols = ["chr", "start", "end"]
    dtypes = {"chr": str, "start": np.int64, "end": np.int64}
    intervals = parse_interval_bed(path)
    return pd.DataFrame([(c, b["start"], b["end"]) for c, rows in intervals.items()
                         for b in rows], columns=cols).astype(dtypes)


def read_params(report_tsv):
    """Parse the report's #params line for max_block_dist, terminal_limit, ultra_fast,
    the canonical motif length and min_canonical_count (the ITS composition floor).

    Each key is converted on its own, so one malformed value cannot abort the whole
    line and silently leave the later keys at their defaults.
    """
    result = {"max_block_dist": 1000, "terminal_limit": None, "ultra_fast": None,
              "motif_len": 6, "min_canonical_count": 4, "label_threshold": LABEL_THRESHOLD,
              "known_params": set()}
    if not report_tsv:
        return result
    try:
        with open(report_tsv) as fh:
            for line in fh:
                if not line.startswith("#params"):
                    continue
                tokens = dict(tok.split("=", 1) for tok in line.split() if "=" in tok)
                if "max_block_dist" in tokens:
                    try:
                        value = int(tokens["max_block_dist"])
                        if value >= 0:
                            result["max_block_dist"] = value
                            result["known_params"].add("max_block_dist")
                    except ValueError:
                        pass
                if "terminal_limit" in tokens:
                    try:
                        value = int(tokens["terminal_limit"])
                        if value > 0:
                            result["terminal_limit"] = value
                            result["known_params"].add("terminal_limit")
                    except ValueError:
                        pass
                if "ultra_fast" in tokens:
                    result["ultra_fast"] = {"true": True, "false": False}.get(tokens["ultra_fast"].lower())
                if "canonical" in tokens:
                    fwd_motif = tokens["canonical"].split("/", 1)[0]
                    if fwd_motif:
                        result["motif_len"] = len(fwd_motif)
                        result["known_params"].add("motif_len")
                if "label_threshold" in tokens:
                    try:
                        value = float(tokens["label_threshold"])
                        if 0.5 <= value <= 1.0:
                            result["label_threshold"] = value
                            result["known_params"].add("label_threshold")
                    except ValueError:
                        pass
                if "min_canonical_count" in tokens:
                    try:
                        value = int(tokens["min_canonical_count"])
                        if value >= 0:
                            result["min_canonical_count"] = value
                            result["known_params"].add("min_canonical_count")
                    except ValueError:
                        pass
                break
    except OSError:
        pass
    return result


_ANOMALY_OF = {
    "discordant_p": "Discordant",
    "discordant_q": "Discordant",
    "fragmented_p": "Fragmented",
    "fragmented_q": "Fragmented",
}

_GAPPED_OF = {
    "T2T": "Gapped T2T",
    "Incomplete": "Gapped Incomplete",
    "Fragmented": "Gapped Fragmented",
    "Discordant": "Gapped Discordant",
    "No telomeres": "Gapped No telomeres",
}

_TYPE_MAP = OrderedDict([
    ("t2t",                "T2T"),
    ("gapped_t2t",         "Gapped T2T"),
    ("incomplete",         "Incomplete"),
    ("gapped_incomplete",  "Gapped Incomplete"),
    ("discordant",         "Discordant"),
    ("gapped_discordant",  "Gapped Discordant"),
    ("none",               "No telomeres"),
    ("gapped_none",        "Gapped No telomeres"),
])


def parse_report(path):
    """Parse *_report.tsv -> OrderedDict[category -> list of chrom names].

    Reads the type column from the Path Summary table written by the C++ tool.
    Skips the Assembly Summary sections (lines starting with +++).
    """
    cats = OrderedDict([
        ("T2T",                 []),
        ("Gapped T2T",          []),
        ("Incomplete",          []),
        ("Gapped Incomplete",   []),
        ("Fragmented",          []),
        ("Gapped Fragmented",   []),
        ("Discordant",          []),
        ("Gapped Discordant",   []),
        ("No telomeres",        []),
        ("Gapped No telomeres", []),
    ])
    header_idx = None
    type_col = None
    header_col = None
    gaps_col = None
    anomaly_col = None
    parsed_rows = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("+++"):
                header_idx = None  # reset on section break
                continue
            parts = line.split("\t")
            if header_idx is None:
                # First non-empty, non-+++ line is the TSV header
                if parts[0] == "pos":
                    try:
                        header_col = parts.index("header")
                        type_col = parts.index("type")
                        gaps_col = parts.index("gaps") if "gaps" in parts else None
                        anomaly_col = parts.index("anomaly") if "anomaly" in parts else None
                    except ValueError:
                        _warn(f"{path}:{lineno}: report header is missing required 'header'/'type' columns; skipping section.")
                        header_col = None
                        type_col = None
                        gaps_col = None
                        anomaly_col = None
                        continue
                    header_idx = 0
                continue
            if type_col is None or header_col is None or len(parts) <= max(header_col, type_col):
                continue
            chrom = parts[header_col]
            raw_type = parts[type_col].lower()
            cat = _TYPE_MAP.get(raw_type)
            # v0.1.6 splits gappedness out of the type, so re-attach it from the gaps column; older reports already carry it in the type string.
            if cat and gaps_col is not None and not raw_type.startswith("gapped_"):
                if len(parts) > gaps_col and parts[gaps_col].isdigit() and int(parts[gaps_col]) > 0:
                    cat = _GAPPED_OF.get(cat, cat)
            gapped = (gaps_col is not None and len(parts) > gaps_col
                      and parts[gaps_col].isdigit() and int(parts[gaps_col]) > 0)
            if anomaly_col is not None and len(parts) > anomaly_col:
                for token in parts[anomaly_col].split(","):
                    flag = _ANOMALY_OF.get(token.strip())
                    if not flag:
                        continue
                    if gapped:
                        flag = _GAPPED_OF.get(flag, flag)
                    if flag in cats and chrom not in cats[flag]:
                        cats[flag].append(chrom)
            if cat and cat in cats:
                cats[cat].append(chrom)
                parsed_rows += 1

    if parsed_rows == 0:
        _warn(f"No usable classification rows were parsed from '{path}'; the overview classification panel will show no data.")

    return OrderedDict((k, v) for k, v in cats.items() if v)


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------


def get_chrom_sizes(blocks, *bedgraph_datasets):
    """Get chromosome sizes from BED pathSize, block ends, or BEDgraph extents."""
    sizes = {}
    for chrom, blist in blocks.items():
        if blist:
            path_sizes = [b["pathSize"] for b in blist if b.get("pathSize", 0) > 0]
            max_end = max(b["end"] for b in blist)
            sizes[chrom] = max(path_sizes) if path_sizes else max_end
            sizes[chrom] = max(sizes[chrom], max_end)
    for bedgraph_data in bedgraph_datasets:
        if bedgraph_data:
            for chrom, (_, ends, _) in bedgraph_data.items():
                if len(ends) > 0:
                    sizes[chrom] = max(sizes.get(chrom, 0), int(ends.max()))
    return sizes


# ---------------------------------------------------------------------------
# ITS analysis: candidate fusion pairing and top hits
# ---------------------------------------------------------------------------

_PAIR_COLUMNS = ["chr", "start", "end", "q_bp", "p_bp", "min_arm", "combined_bp",
                 "spacer_bp", "q_can_prop", "p_can_prop", "pos_frac"]


def pair_fusions(df, gaps, d):
    """Rebuild q->p pairs within d bp (no N-gap between), ranked by the shorter arm.

    A candidate is kept only when at least one of the two rows already carries
    teloType=="fusion" (src/teloscope.cpp:300 assignJunctions), so a pair never
    contradicts the class the engine itself assigned; its partner can still be
    labelled fragmentation/single if that row's *own* nearer neighbour differs.
    """
    if df.empty:
        return pd.DataFrame(columns=_PAIR_COLUMNS)

    ordered = df.sort_values(["chr", "start", "end", "teloLabel", "teloType", "teloLen", "canonical_bp"],
                             kind="mergesort").reset_index(drop=True)
    nxt = ordered.groupby("chr", sort=False).shift(-1)

    spacer = nxt["start"] - ordered["end"]
    mask = (
        (ordered["teloLabel"] == "q")
        & (nxt["teloLabel"] == "p")
        & nxt["start"].notna()
        & (spacer <= d)
        & ((ordered["teloType"] == "fusion") | (nxt["teloType"] == "fusion"))
    )
    if not mask.any():
        return pd.DataFrame(columns=_PAIR_COLUMNS)

    cand = pd.DataFrame({
        "chr":          ordered.loc[mask, "chr"].to_numpy(),
        "start":        ordered.loc[mask, "start"].to_numpy(),
        "q_end":        ordered.loc[mask, "end"].to_numpy(),
        "p_start":      nxt.loc[mask, "start"].to_numpy().astype(np.int64),
        "end":          nxt.loc[mask, "end"].to_numpy().astype(np.int64),
        "q_bp":         ordered.loc[mask, "teloLen"].to_numpy(),
        "p_bp":         nxt.loc[mask, "teloLen"].to_numpy().astype(np.int64),
        "q_can_prop":   ordered.loc[mask, "can_prop"].to_numpy(),
        "p_can_prop":   nxt.loc[mask, "can_prop"].to_numpy(),
        "chrSize":      ordered.loc[mask, "chrSize"].to_numpy(),
    })

    # drop pairs whose gap spans a real N-gap; searchsorted per chromosome, not per pair
    drop = np.zeros(len(cand), dtype=bool)
    if gaps is not None and len(gaps):
        gap_intervals = {}
        for chrom, grp in gaps.groupby("chr", sort=False):
            grp = grp.sort_values("start")
            gap_intervals[chrom] = (grp["start"].to_numpy(), np.maximum.accumulate(grp["end"].to_numpy()))
        for chrom, grp in cand.groupby("chr", sort=False):
            interval = gap_intervals.get(chrom)
            if interval is None:
                continue
            starts, ends = interval
            idx = np.searchsorted(starts, grp["p_start"].to_numpy(), side="left") - 1
            overlap = ((idx >= 0) & (ends[np.maximum(idx, 0)] > grp["q_end"].to_numpy())
                       & (grp["p_start"].to_numpy() > grp["q_end"].to_numpy()))
            drop[grp.index[overlap]] = True
    cand = cand[~drop]

    cand["spacer_bp"] = np.clip(cand["p_start"] - cand["q_end"], 0, None).astype(np.int64)
    cand["min_arm"] = np.minimum(cand["q_bp"], cand["p_bp"])
    cand["combined_bp"] = cand["q_bp"] + cand["p_bp"]
    cand["pos_frac"] = np.where(cand["chrSize"] > 0,
                                (cand["start"] + cand["end"]) / 2.0 / cand["chrSize"], np.nan)

    cand = cand.sort_values(["min_arm", "combined_bp", "chr", "start", "end"],
                            ascending=[False, False, True, True, True]).reset_index(drop=True)
    return cand[_PAIR_COLUMNS]


_LONG_ITS_COLUMNS = ["chr", "start", "end", "teloLen", "canonical_bp", "can_prop",
                     "teloLabel", "teloType", "pos_frac"]
_CLUSTER_COLUMNS = ["chr", "start", "end", "rows", "span", "its_bp", "canonical_bp"]


def rank_long_its(df, top_n=25):
    """Long ITS rows ranked by canonical bp descending, ties broken by length descending."""
    if df.empty:
        return df[_LONG_ITS_COLUMNS] if set(_LONG_ITS_COLUMNS) <= set(df.columns) else df
    ranked = df.sort_values(["canonical_bp", "teloLen", "chr", "start", "end", "teloLabel", "teloType"],
                           ascending=[False, False, True, True, True, True, True]).reset_index(drop=True)
    return ranked.head(top_n)[_LONG_ITS_COLUMNS]


def its_top_hits(out_path, df, pairs, clusters, top_longest=25):
    """Write <prefix>_its_top_hits.tsv: candidate fusions, long ITS by canonical bp, then clusters."""
    with open(out_path, "w") as fh:
        fh.write("# teloscope ITS top hits\n")
        fh.write("# section 1: candidate fusion pairs (q->p), ranked by min(q_bp, p_bp) "
                "descending, ties broken by combined_bp descending\n")
        fh.write("#" + "\t".join(_PAIR_COLUMNS) + "\n")
        pairs[_PAIR_COLUMNS].to_csv(fh, sep="\t", header=False, index=False, float_format="%.4f")

        fh.write(f"# section 2: {top_longest} longest ITS by canonical bp "
                "descending, ties broken by teloLen descending\n")
        fh.write("#chr\tstart\tend\tteloLen\tcanonical_bp\tcan_prop\tlabel\tclass\tpos_frac\n")
        rank_long_its(df, top_longest).to_csv(fh, sep="\t", header=False, index=False, float_format="%.4f")

        fh.write("# section 3: ITS clusters (scaffold rows within 50 kb of each other) with "
                ">= 3 rows, ranked by summed ITS bp descending\n")
        fh.write("#chr\tstart\tend\trows\tspan\tits_bp\tcanonical_bp\n")
        clusters[_CLUSTER_COLUMNS].to_csv(fh, sep="\t", header=False, index=False, float_format="%.4f")


def compute_its_clusters(df, merge_gap=ITS_CLUSTER_MERGE_GAP):
    """Assign a cluster id to each ITS row: same scaffold, merged while the gap to the
    running-max end of earlier rows on that scaffold is <= merge_gap (vectorised; shared
    by the ITS-1 ideogram, the ITS-2/TSV cluster tables, and plot_its.py's auto-window)."""
    if df.empty:
        out = df.copy()
        out["cluster_id"] = pd.Series(dtype=np.int64)
        return out
    ordered = df.sort_values(["chr", "start"], kind="mergesort").reset_index(drop=True)
    running_end = ordered.groupby("chr")["end"].cummax().shift(1)
    new_chr = ordered["chr"] != ordered["chr"].shift(1)
    gap = ordered["start"] - running_end
    new_cluster = new_chr | gap.isna() | (gap > merge_gap)
    ordered["cluster_id"] = new_cluster.cumsum()
    return ordered


def summarize_its_clusters(df, merge_gap=ITS_CLUSTER_MERGE_GAP, min_rows=ITS_CLUSTER_MIN_ROWS):
    """Per-cluster chr/start/end/rows/span/its_bp/canonical_bp, ranked by its_bp descending."""
    if df.empty:
        return pd.DataFrame(columns=_CLUSTER_COLUMNS)
    clustered = compute_its_clusters(df, merge_gap)
    if "canonical_bp" not in clustered.columns:
        clustered = clustered.assign(canonical_bp=0)
    agg = clustered.groupby("cluster_id").agg(
        chr=("chr", "first"), start=("start", "min"), end=("end", "max"),
        rows=("chr", "size"), its_bp=("teloLen", "sum"), canonical_bp=("canonical_bp", "sum"))
    agg = agg[agg["rows"] >= min_rows].copy()
    agg["span"] = agg["end"] - agg["start"]
    agg = agg.sort_values(["its_bp", "chr", "start", "end"],
                          ascending=[False, True, True, True]).reset_index(drop=True)
    return agg[_CLUSTER_COLUMNS]


def pad_window(start, end, chrom_size, pad=None):
    """Symmetric bp padding around [start, end), clamped to the scaffold.

    Default padding is 2x the span (floored at 10 kb), matching the auto-centred single-
    locus zoom windows in plot_its.py.
    """
    span = max(end - start, 1)
    use_pad = pad if pad is not None else max(2 * span, 10_000)
    return max(0, int(start - use_pad)), min(int(chrom_size), int(end + use_pad))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _format_bp_axis(ax, max_bp):
    """Set x-axis labels in bp / kb / Mb depending on scale."""
    if max_bp >= 5_000_000:
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f}"))
        ax.set_xlabel("Position (Mb)")
    elif max_bp >= 5_000:
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f}"))
        ax.set_xlabel("Position (kb)")
    else:
        ax.set_xlabel("Position (bp)")


def _pick_bp_unit(max_bp):
    """Pick one bp/kb/Mb/Gb unit for a value or panel from its largest value."""
    if max_bp >= 1_000_000_000:
        return 1e9, "Gb"
    if max_bp >= 1_000_000:
        return 1e6, "Mb"
    if max_bp >= 1_000:
        return 1e3, "kb"
    return 1.0, "bp"


def _mb_tick_decimals(span_mb):
    """Decimal places for an Mb-axis tick label, so nearby ticks never round to duplicates."""
    return max(0, int(np.ceil(np.log10(6.0 / max(span_mb, 1e-12)))))


def _fmt_bp(bp, _pos=None):
    """Format a bp value in its own best-fitting unit, trimming trailing zeros."""
    if bp <= 0:
        return "0"
    divisor, unit = _pick_bp_unit(bp)
    if unit == "bp":
        return f"{int(round(bp)):,} bp"
    decimals = 2 if unit in ("Mb", "Gb") else 1
    text = f"{bp / divisor:.{decimals}f}".rstrip("0").rstrip(".")
    return f"{text} {unit}"


def _fmt_bp_fixed(bp):
    """Format a bp value in its own unit with fixed decimals per unit, so labels line up."""
    if bp <= 0:
        return "0"
    divisor, unit = _pick_bp_unit(bp)
    decimals = {"Gb": 2, "Mb": 2, "kb": 1}.get(unit, 0)
    return f"{bp / divisor:,.{decimals}f} {unit}"


def _plain_number(value):
    """Tick text without float noise or trailing zeros: 120, 2.5, 0."""
    return np.format_float_positional(round(float(value), 6), trim="-")


def _fmt_kbp(bp):
    """Format a base-pair distance as compact kbp text."""
    kbp = bp / 1e3
    if kbp >= 100:
        return f"{kbp:.0f}"
    if kbp >= 10:
        return f"{kbp:.1f}".rstrip("0").rstrip(".")
    return f"{kbp:.2f}".rstrip("0").rstrip(".")


def _median(values):
    """Median via numpy."""
    return float(np.median(values))


def _bedgraph_to_step(starts, ends, values):
    """Convert BEDgraph arrays to step-plot coordinates."""
    n = len(starts)
    if n == 0:
        return np.array([]), np.array([])
    xs = np.empty(n + 1)
    ys = np.empty(n + 1)
    xs[:n] = starts
    ys[:n] = np.clip(values, 0, None)
    xs[n] = ends[-1]
    ys[n] = values[-1]
    return xs, ys


def _panel_label(ax, label, x_in=None, y=1):
    """Bold panel letter at a page x in inches (default: a fixed gap left of the axes), on the title baseline."""
    fig = ax.figure
    if x_in is None:
        x_in = ax.get_position().x0 * fig.get_size_inches()[0] - PANEL_LETTER_DX_IN
    trans = (transforms.blended_transform_factory(fig.dpi_scale_trans, ax.transAxes)
             + transforms.ScaledTranslation(0, 3 / 72, fig.dpi_scale_trans))
    ax.text(x_in, y, label, transform=trans, fontsize=PANEL_LABEL_SIZE, fontweight="bold",
            va="baseline", ha="left")


def _panel_title(ax, title, label=None, x_in=None, y=None):
    """Axis-attached panel title with an optional panel letter, as on the terminal page; y lifts both in axes units."""
    ax.set_title(title, fontsize=PANEL_TITLE_SIZE, pad=3, y=y)
    if label:
        _panel_label(ax, label, x_in=x_in, y=1 if y is None else y)


def _page_title(fig, title):
    """One bold page title at a fixed physical distance from the top edge; no subtitle."""
    fig.suptitle(title, fontsize=PANEL_LABEL_SIZE, fontweight="bold",
                 y=1 - 0.155 / fig.get_size_inches()[1])


def _thirds(num, den, threshold=LABEL_THRESHOLD):
    """Engine symmetric-thirds rule (computeStrandLabel): 2 above t, 0 below 1-t, else 1."""
    num = np.asarray(num, dtype=np.int64)
    den = np.asarray(den, dtype=np.int64)
    scale = 1_000_000
    t = int(round(threshold * scale))
    out = np.ones(np.broadcast(num, den).shape, dtype=np.int64)
    out = np.where(num * scale > den * t, 2, out)
    out = np.where(num * scale < den * (scale - t), 0, out)
    return np.where(den == 0, 1, out)


def composition_class(fwd_can, rev_can, fwd_noncan, rev_noncan, threshold=LABEL_THRESHOLD):
    """Double-key class per ITS: (strand 0 rev/1 both/2 fwd, canonicity 0 nonCan/1 mixed/2 can)."""
    fwd_can, rev_can, fwd_noncan, rev_noncan = (np.asarray(v, dtype=np.int64)
                                                for v in (fwd_can, rev_can, fwd_noncan, rev_noncan))
    total = fwd_can + rev_can + fwd_noncan + rev_noncan
    return (_thirds(fwd_can + fwd_noncan, total, threshold),
            _thirds(fwd_can + rev_can, total, threshold))


def double_key_colors(df, threshold=LABEL_THRESHOLD):
    """Hex colour per ITS row from its four match counts."""
    strand, canon = composition_class(df["fwdCan"], df["revCan"], df["fwdNonCan"], df["revNonCan"], threshold)
    return [DOUBLE_KEY_COLORS[c][s] for s, c in zip(strand, canon)]


def _draw_double_key(ax, counts=None, bp=None, fontsize=MIN_TEXT_SIZE):
    """3x3 composition key, fwd drawn left (5' end); with counts[c][s] (and bp) each cell carries n and % of ITS bp."""
    total_bp = float(np.sum(bp)) if bp is not None else 0.0
    for c in range(3):
        for s in range(3):
            ax.add_patch(Rectangle((2 - s, c), 1, 1, facecolor=DOUBLE_KEY_RGBA[c][s], edgecolor="white", linewidth=0.8))
            if counts is None:
                continue
            n = int(counts[c][s])
            share = 100 * bp[c][s] / total_bp if total_bp > 0 else 0
            text = f"{n:,}" if total_bp <= 0 or not n else f"{n:,}\n{share:.0f}%" if share >= 1 else f"{n:,}\n<1%"
            ax.text(2.5 - s, c + 0.5, text, ha="center", va="center", fontsize=fontsize,
                    color="white" if c == 2 else "#222222", linespacing=0.95)
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 3)
    ax.set_aspect("equal")
    ax.set_xticks([0.5, 1.5, 2.5], DOUBLE_KEY_STRAND[::-1])
    ax.set_yticks([0.5, 1.5, 2.5], DOUBLE_KEY_CANON)
    ax.tick_params(length=0, labelsize=fontsize, pad=1.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel("← Forward share", fontsize=fontsize, labelpad=1.5)
    ax.set_ylabel("Canonical share →", fontsize=fontsize, labelpad=1.5)


def _sanitize_filename(name):
    """Convert a chromosome name into a filesystem-safe stem."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return safe.strip("._") or "unnamed"


def _hide_x_axis(ax):
    """Hide x-axis ticks and labels for non-bottom tracks."""
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines["bottom"].set_visible(False)


def _set_track_label(ax, label, x=TRACK_LABEL_X):
    """Keep labels centered to their track while right-aligning multiline text."""
    if label:
        ax.set_ylabel(label, fontsize=AXIS_LABEL_SIZE, rotation=0,
                      ha="right", va="center")
        ax.yaxis.label.set_multialignment("right")
        ax.yaxis.label.set_linespacing(0.95)
        ax.yaxis.set_label_coords(x, 0.5)
    else:
        ax.set_ylabel("")


def _block_symbol_fits(block_width, view_span, label):
    """Return whether a block is wide enough for its direction glyph."""
    min_fraction = BLOCK_GLYPH_MIN_FRACTION.get(label)
    if min_fraction is None or view_span <= 0:
        return False
    return (block_width / float(view_span)) >= min_fraction


def _draw_block_symbol(ax, x_pos, y_pos, label, zorder):
    """Draw a thin, editable text glyph for the block direction."""
    glyph = BLOCK_GLYPHS.get(label)
    if glyph is None:
        return
    ax.text(
        x_pos,
        y_pos,
        glyph,
        ha="center",
        va="center",
        fontsize=BLOCK_GLYPH_SIZE,
        color="white",
        zorder=zorder,
        clip_on=True,
    )


def _style_fraction_axis(ax, label=None, y_max=1.0, y_min=0.0):
    """Style quantitative tracks with minimal ticks."""
    margin = (y_max - y_min) * 0.02
    ax.set_ylim(y_min - margin, y_max + margin)
    mid = (y_min + y_max) / 2
    ax.set_yticks([y_min, mid, y_max])
    min_str = f"{y_min:g}"
    mid_str = f"{mid:g}"
    max_str = f"{y_max:g}"
    ax.set_yticklabels([min_str, mid_str, max_str])
    ax.tick_params(axis="y", length=2.0, width=0.45, pad=2.2,
                   labelsize=AXIS_TICK_SIZE, colors="#555555")
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.45)
    ax.spines["left"].set_color("#bcbcbc")
    _set_track_label(ax, label)
    _hide_x_axis(ax)


def _terminal_distance_interval(start, end, chrom_size, arm):
    """Return interval coordinates expressed as distance from the relevant end."""
    if arm == "p":
        return max(start, 0), max(end, 0)
    if arm == "q":
        return max(chrom_size - end, 0), max(chrom_size - start, 0)
    return start, end


def _project_terminal_interval(start, end, chrom_size, arm):
    """Project genomic coordinates into panel coordinates for terminal plots."""
    if arm in {"p", "q"}:
        return _terminal_distance_interval(start, end, chrom_size, arm)
    return start, end


def _project_terminal_series(starts, ends, values, chrom_size, arm):
    """Project BEDgraph intervals into panel coordinates and sort left-to-right."""
    starts = np.asarray(starts, dtype=np.float64)
    ends = np.asarray(ends, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if arm != "q":
        return starts, ends, values

    proj_starts = chrom_size - ends
    proj_ends = chrom_size - starts
    order = np.argsort(proj_starts)
    return proj_starts[order], proj_ends[order], values[order]


def _format_terminal_tick_value(value_kbp):
    """Format terminal-axis tick labels in kbp without spurious rounding."""
    return f"{value_kbp:.2f}".rstrip("0").rstrip(".")


def _get_terminal_axis_spec(view_start, view_end, arm):
    """Return x-axis limits, tick positions, labels, and xlabel for terminal plots."""
    if view_end <= view_start:
        return None

    n_ticks = 5
    if arm in ("p", "q"):
        axis_span_bp = max(float(view_end - view_start), 1.0)
        axis_span_kbp = axis_span_bp / 1e3
        axis_start = view_start
        axis_end = view_end
        tick_pos = np.linspace(axis_start, axis_end, n_ticks)
        tick_values = np.linspace(0.0, float(axis_span_kbp), n_ticks)
        tick_labels = [_format_terminal_tick_value(value) for value in tick_values]
        xlabel = "Distance to end (kbp)"
        return {
            "axis_start": axis_start,
            "axis_end": axis_end,
            "tick_pos": tick_pos,
            "tick_labels": tick_labels,
            "xlabel": xlabel,
        }

    span = view_end - view_start
    if span >= 2_000_000:
        fmt = lambda x: f"{x / 1e6:.1f}"
        unit = "Mbp"
    elif span >= 2_000:
        fmt = lambda x: f"{int(round(max(x, 0) / 1e3))}"
        unit = "kbp"
    else:
        fmt = lambda x: f"{max(x, 0):.0f}"
        unit = "bp"

    tick_pos = np.linspace(view_start, view_end, n_ticks)
    tick_labels = [fmt(pos) for pos in tick_pos]

    xlabel = f"Position along scaffold ({unit})"

    return {
        "axis_start": view_start,
        "axis_end": view_end,
        "tick_pos": tick_pos,
        "tick_labels": tick_labels,
        "xlabel": xlabel,
    }


def _apply_terminal_x_axis(ax, axis_spec, arm):
    """Apply end-relative x-axis with adaptive units (bp / kbp / Mbp)."""
    if axis_spec is None:
        return

    if arm == "q":
        ax.set_xlim(axis_spec["axis_end"], axis_spec["axis_start"])
    else:
        ax.set_xlim(axis_spec["axis_start"], axis_spec["axis_end"])

    ax.set_xticks(axis_spec["tick_pos"])
    ax.set_xticklabels(axis_spec["tick_labels"])
    ax.tick_params(axis="x", bottom=True, labelbottom=True,
                   length=2.3, width=0.45, pad=1.5)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(0.45)
    ax.spines["bottom"].set_color("#bcbcbc")
    ax.set_xlabel(axis_spec["xlabel"], fontsize=AXIS_LABEL_SIZE, labelpad=2.5)
    ax.minorticks_off()


def _iter_true_runs(mask):
    """Yield contiguous [start, stop) runs where mask is true."""
    run_start = None
    for idx, flag in enumerate(mask):
        if flag and run_start is None:
            run_start = idx
        elif not flag and run_start is not None:
            yield run_start, idx
            run_start = None
    if run_start is not None:
        yield run_start, len(mask)


def _block_end_distance(block, chrom_size):
    """Return the relevant distance from a telomere block to the scaffold end in bp."""
    closest_end = block.get("closestEnd")
    if closest_end in {"p", "q"}:
        start_dist, end_dist = _terminal_distance_interval(
            int(block["start"]), int(block["end"]), int(chrom_size), closest_end)
        return min(start_dist, end_dist)
    left_gap = max(int(block["start"]), 0)
    right_gap = max(int(chrom_size) - int(block["end"]), 0)
    return min(left_gap, right_gap)


def _classify_outliers(values):
    """Return a boolean mask of Tukey outliers; conservative for small samples."""
    values = np.asarray(values, dtype=np.float64)
    if len(values) < 4 or np.ptp(values) == 0:
        return np.zeros(len(values), dtype=bool)

    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    if iqr <= 0:
        return np.zeros(len(values), dtype=bool)

    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    return (values < lo) | (values > hi)


_ORIENT_KW_CACHE = {}

def _orient_kw(fn, vert):
    """Orientation kwargs: 'orientation' on matplotlib >= 3.10, 'vert' before."""
    fn_name = fn.__name__
    if fn_name not in _ORIENT_KW_CACHE:
        has_orientation = "orientation" in inspect.signature(fn).parameters
        _ORIENT_KW_CACHE[fn_name] = has_orientation

    if _ORIENT_KW_CACHE[fn_name]:
        return {"orientation": "vertical" if vert else "horizontal"}
    else:
        return {"vert": vert}


def _draw_raincloud_group(ax, values, position, color, rng, vert=False):
    """Draw a half-violin + boxplot + jittered points group.

    vert=False (default): horizontal — values on x, group position on y.
    vert=True:            vertical   — values on y, group position on x.
    """
    values = np.asarray(values, dtype=np.float64)
    if len(values) == 0:
        return

    outliers = _classify_outliers(values)
    core_values = values[~outliers]
    if len(core_values) < 2:
        core_values = values

    if len(core_values) >= 2 and np.ptp(core_values) > 0:
        violin_kw = _orient_kw(ax.violinplot, vert)
        violin = ax.violinplot(
            [core_values], positions=[position],
            widths=0.28, showmeans=False, showmedians=False,
            showextrema=False,
            **violin_kw
        )
        body = violin["bodies"][0]
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_linewidth(0.5)
        body.set_alpha(0.15)

        verts = body.get_paths()[0].vertices
        if vert:
            # Clip to right half so violin sits to the right of the scatter.
            verts[:, 0] = np.maximum(verts[:, 0], position)
        else:
            # Clip to upper half so violin sits above the scatter.
            verts[:, 1] = np.maximum(verts[:, 1], position)

    jitter = rng.uniform(-0.13, -0.07, size=len(values))
    filled_mask = ~outliers
    if np.any(filled_mask):
        xs = (np.full(np.sum(filled_mask), position) + jitter[filled_mask]
              if vert else values[filled_mask])
        ys = (values[filled_mask]
              if vert else np.full(np.sum(filled_mask), position) + jitter[filled_mask])
        ax.scatter(xs, ys, s=8.5, color=color, alpha=0.48,
                   edgecolors="white", linewidths=0.22, rasterized=True, zorder=3)
    if np.any(outliers):
        xs = (np.full(np.sum(outliers), position) + jitter[outliers]
              if vert else values[outliers])
        ys = (values[outliers]
              if vert else np.full(np.sum(outliers), position) + jitter[outliers])
        ax.scatter(xs, ys, s=8.5, color="#e6e6e6",
                   edgecolors="none", linewidths=0.0,
                   rasterized=True, zorder=4)

    boxplot_kw = _orient_kw(ax.boxplot, vert)
    box = ax.boxplot(
        [core_values],
        positions=[position],
        widths=0.042,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        manage_ticks=False,
        **boxplot_kw
    )
    for patch in box["boxes"]:
        patch.set_facecolor("white")
        patch.set_edgecolor(color)
        patch.set_linewidth(0.55)
    for key in ("whiskers", "caps", "medians"):
        for artist in box[key]:
            artist.set_color(color)
            artist.set_linewidth(0.55)


def _draw_balanced_legend(ax, handles, y_anchor=0.50, ncol=None, fontsize=LEGEND_TEXT_SIZE,
                          loc="center", bbox_to_anchor=None, alignment=None):
    """Render a compact two-row legend within a dedicated legend axis."""
    ax.axis("off")
    if not handles:
        return None

    ncol = ncol or max(1, int(np.ceil(len(handles) / 2.0)))
    if bbox_to_anchor is None:
        bbox_to_anchor = (0.5, y_anchor)
    leg = ax.legend(
        handles=handles,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        ncol=ncol,
        fontsize=fontsize,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        framealpha=1.0,
        fancybox=False,
        handlelength=0.88,
        handleheight=0.7,
        columnspacing=0.90,
        handletextpad=0.28,
        borderaxespad=0.0,
        labelspacing=0.42,
    )
    leg.get_frame().set_linewidth(0.4)
    if alignment is not None and hasattr(leg, "_legend_box"):
        leg._legend_box.align = alignment
    return leg


def _draw_centered_legend_pair(ax, left_handles, right_handles, gap=0.016,
                               left_legend_kwargs=None, right_legend_kwargs=None):
    """Draw two boxed legends as a compact centered pair within one shared strip."""
    ax.axis("off")
    left_legend_kwargs = dict(left_legend_kwargs or {})
    right_legend_kwargs = dict(right_legend_kwargs or {})

    left_legend_kwargs.setdefault("loc", "center")
    left_legend_kwargs.setdefault("bbox_to_anchor", (0.5, 0.50))
    left_legend_kwargs.setdefault("alignment", "center")
    left_legend = _draw_balanced_legend(ax, left_handles, **left_legend_kwargs)

    right_legend_kwargs.setdefault("loc", "center")
    right_legend_kwargs.setdefault("bbox_to_anchor", (0.5, 0.50))
    right_legend_kwargs.setdefault("alignment", "center")
    right_legend = _draw_balanced_legend(ax, right_handles, **right_legend_kwargs)

    if left_legend is None or right_legend is None:
        if left_legend is not None:
            ax.add_artist(left_legend)
        return left_legend, right_legend

    ax.add_artist(left_legend)

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    left_box = left_legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    right_box = right_legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    axis_box = ax.get_position()

    group_width = left_box.width + gap + right_box.width
    group_left = axis_box.x0 + ((axis_box.width - group_width) / 2.0)
    left_center = group_left + (left_box.width / 2.0)
    right_center = group_left + left_box.width + gap + (right_box.width / 2.0)

    left_anchor = (left_center - axis_box.x0) / axis_box.width
    right_anchor = (right_center - axis_box.x0) / axis_box.width
    left_legend.set_bbox_to_anchor((left_anchor, left_legend_kwargs["bbox_to_anchor"][1]), transform=ax.transAxes)
    right_legend.set_bbox_to_anchor((right_anchor, right_legend_kwargs["bbox_to_anchor"][1]), transform=ax.transAxes)
    return left_legend, right_legend


def _draw_ranked_bar_panel(ax, rows, x_label, xticks, title=None, xticklabels=None):
    """Ranked horizontal bars in ten fixed slots, so thickness and pitch never depend on n."""
    label_in = 0.80
    n_slots = max(len(rows), 10)
    y_pos = np.arange(len(rows), dtype=np.float64)
    values = [row["value"] for row in rows]
    colors = [row["color"] for row in rows]

    ax.barh(y_pos, values, height=0.6, color=colors, edgecolor="none", zorder=2)
    for y_pos_i, row in zip(y_pos, rows):
        lines = row["label"].split("\n")
        if len(lines) > 2:
            lines = [lines[0], lines[1].rstrip("._-") + "…"]
        ax.annotate("\n".join(lines), (0.0, y_pos_i), xytext=(-3, 0), textcoords="offset points",
                    ha="right", multialignment="right", va="center",
                    fontsize=MIN_TEXT_SIZE, color="#222222", linespacing=0.95)

    x_max = max(max(xticks), max(values) + 0.12) if values else max(xticks)
    # Fixed physical label gutter inside the axes, whatever the panel width.
    width_in = ax.get_position().width * ax.figure.get_size_inches()[0]
    ax.set_xlim(-x_max * label_in / max(width_in - label_in, 0.1), x_max)
    ax.set_ylim(n_slots - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xlabel(x_label, fontsize=AXIS_LABEL_SIZE, labelpad=2.5)
    ax.set_xticks([tick for tick in xticks if tick <= x_max + 0.1])
    if xticklabels:
        ax.set_xticklabels(xticklabels[:len(ax.get_xticks())])
    ax.tick_params(axis="x", length=2.3, width=0.45, pad=1.5, labelsize=AXIS_TICK_SIZE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(0.45)
    ax.spines["bottom"].set_color("#bcbcbc")
    ax.spines["bottom"].set_bounds(0.0, x_max)
    ax.axvline(0.0, color="#bcbcbc", linewidth=0.45, zorder=1)
    ax.minorticks_off()
    ax.xaxis.label.set_x(1 - x_max / 2 / (x_max - ax.get_xlim()[0]))  # centred on the bars, not the label gutter
    if title:
        _panel_title(ax, title)
        ax.title.set_x(ax.xaxis.label.get_position()[0])


def _draw_dual_summary_panel(ax, all_segments, telo_segments, title=None):
    """Two adjacent stacked percentage bars sharing one axis, pitched like the ranked bars."""
    plot_rows = [(0, "All", all_segments), (1, "With telomeres", telo_segments)]
    for y_pos, _, segments in plot_rows:
        left = 0.0
        for segment in segments:
            _, value, color = segment[:3]
            cnt = segment[3] if len(segment) > 3 else None
            if value <= 0:
                continue
            ax.barh(y_pos, value, left=left, height=0.5, color=color,
                    edgecolor="white", linewidth=0.6, zorder=2)
            if cnt is not None and value >= 5.0:
                r, g, b = matplotlib.colors.to_rgb(color)
                ink = "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.45 else "#222222"
                ax.text(left + value / 2.0, y_pos, str(cnt), ha="center", va="center",
                        fontsize=ANNOTATION_TEXT_SIZE, color=ink, zorder=3)
            left += value

    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(2.0, -1.0)
    ax.set_yticks([row[0] for row in plot_rows])
    ax.set_yticklabels([row[1] for row in plot_rows], fontsize=AXIS_TICK_SIZE)
    ax.tick_params(axis="y", length=0, pad=3.0)
    ax.tick_params(axis="x", length=2.3, width=0.45, pad=1.5, labelsize=AXIS_TICK_SIZE)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}%"))
    ax.set_xlabel("Sequences (%)", fontsize=AXIS_LABEL_SIZE, labelpad=2.5)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#bcbcbc")
    ax.spines["bottom"].set_linewidth(0.45)
    ax.minorticks_off()
    if title:
        _panel_title(ax, title)


def _wrap_scaffold_name(name, width=18):
    """Wrap long scaffold names across multiple lines while preserving delimiters."""
    parts = re.findall(r"[^._-]+[._-]*", name) or [name]
    lines = []
    current = ""

    for part in parts:
        while len(part) > width:
            if current:
                lines.append(current)
                current = ""
            lines.append(part[:width])
            part = part[width:]
        if not current:
            current = part
        elif len(current) + len(part) <= width:
            current += part
        else:
            lines.append(current)
            current = part

    if current:
        lines.append(current)
    return "\n".join(lines)


def _short_exception(exc):
    """Compact exception string for warnings and placeholder figures."""
    message = str(exc).strip()
    return f"{type(exc).__name__}: {message}" if message else type(exc).__name__


def _placeholder_figure(title, message):
    """Simple fallback figure used when a panel cannot be rendered."""
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, REPORT_PAGE_HEIGHT))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.01, 0.92, title, transform=ax.transAxes,
            fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="top", ha="left")
    ax.text(0.01, 0.72, textwrap.fill(message, width=95), transform=ax.transAxes,
            fontsize=FIGURE_SUMMARY_SIZE, va="top", ha="left", color="#444444")
    fig.subplots_adjust(left=0.03, right=0.98, top=0.95, bottom=0.08)
    return fig


def _save_figure_with_fallback(save_figure, build_figure, title, message_prefix):
    """Build and save a figure, falling back to a placeholder page on render errors."""
    fig = None
    existing_figures = set(plt.get_fignums())
    try:
        fig = build_figure()
        save_figure(fig)
        return True, None
    except Exception as exc:
        error_text = _short_exception(exc)
        _warn(f"{title}: {error_text}")
        if fig is not None:
            plt.close(fig)
            fig = None
        fallback = _placeholder_figure(title, f"{message_prefix}\n\n{error_text}")
        try:
            save_figure(fallback)
        finally:
            plt.close(fallback)
        return False, error_text
    finally:
        if fig is not None:
            plt.close(fig)
        # A builder can fail after allocating a figure but before returning it.
        for number in set(plt.get_fignums()) - existing_figures:
            plt.close(number)


# ---------------------------------------------------------------------------
# Tier 1: Assembly overview (two pages)
# ---------------------------------------------------------------------------

def _draw_flagged_scaffolds_panel(ax, classifications, chrom_sizes, top_n=10, title=None):
    """Horizontal bars of flagged scaffolds ranked by scaffold size."""
    flagged = []
    for cat in FLAGGED_SCAFFOLD_CATEGORIES:
        for chrom in classifications.get(cat, []):
            size = chrom_sizes.get(chrom, 0)
            if size > 0:
                flagged.append({"chrom": chrom, "category": cat, "size": size})
    flagged.sort(key=lambda x: x["size"], reverse=True)
    flagged_top = flagged[:top_n]

    if not flagged_top:
        _draw_empty_state(ax, "No flagged scaffolds\n(fragmented / discordant)")
        if title:
            _panel_title(ax, title)
        return

    # Log bars from a 1 kbp baseline, ticked at plain Mbp decades.
    rows = [{
        "label": _wrap_scaffold_name(entry["chrom"], width=14),
        "value": max(np.log10(entry["size"] / 1e3), 0.05),  # sub-kbp scaffolds keep a sliver
        "color": COLORS.get(entry["category"], "#aaaaaa"),
    } for entry in flagged_top]
    xticks = list(range(max(int(np.ceil(rows[0]["value"] + 0.12)), 3) + 1))
    _draw_ranked_bar_panel(ax, rows, x_label="Scaffold size (Mbp)", xticks=xticks, title=title,
                           xticklabels=[f"{10.0 ** (tick - 3):g}" for tick in xticks])


def _compute_block_rows(blocks, chrom_sizes):
    """Build the flat block_rows list shared by both overview pages."""
    rows = []
    for chrom, blist in blocks.items():
        chrom_size = chrom_sizes.get(chrom, 0)
        if chrom_size <= 0 and blist:
            chrom_size = max(
                max((b.get("pathSize", 0) for b in blist), default=0),
                max(b["end"] for b in blist),
            )
        for block in blist:
            rows.append({
                "chrom": chrom,
                "label": block["label"],
                "closestEnd": block.get("closestEnd"),
                "length": block["length"],
                "distance": _block_end_distance(block, chrom_size),
            })
    return rows


def _add_axes_in(fig, left, top, width, height):
    """Add axes placed in inches from the figure's top-left corner."""
    fig_w, fig_h = fig.get_size_inches()
    return fig.add_axes([left / fig_w, 1 - (top + height) / fig_h, width / fig_w, height / fig_h])


def _draw_empty_state(ax, message):
    """Blank the axes and centre a muted message in it."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center",
            fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999", linespacing=1.4)


def _draw_stat_tiles(ax, tiles):
    """Row of headline numbers over muted labels, split by hairlines."""
    ax.axis("off")
    ax.set_xlim(0, len(tiles))
    ax.set_ylim(0, 1)
    for i, (value, label) in enumerate(tiles):
        ax.text(i + 0.5, 0.50, value, ha="center", va="bottom",
                fontsize=PANEL_LABEL_SIZE, color="#222222")
        ax.text(i + 0.5, 0.38, label, ha="center", va="top",
                fontsize=MIN_TEXT_SIZE, color="#666666")
        if i:
            ax.plot([i, i], [0.08, 0.92], color="#d0d0d0", linewidth=0.4,
                    solid_capstyle="butt")


def plot_overview_page1(classifications, blocks, chrom_sizes):
    """Page 1: headline tiles, scaffold classification (a) and flagged scaffolds (b)."""
    FLAGGED_DISTANCE_BP = 1000

    cat_labels = list(classifications.keys())
    cat_counts = [len(v) for v in classifications.values()]
    cat_colors = [COLORS.get(c, "#aaaaaa") for c in cat_labels]
    total = sum(cat_counts)

    block_rows = _compute_block_rows(blocks, chrom_sizes)
    flagged_scaffolds = {
        chrom
        for category in FLAGGED_SCAFFOLD_CATEGORIES
        for chrom in classifications.get(category, [])
    }
    all_len = [r["length"] for r in block_rows if r["label"] in ("p", "q", "b")]
    total_paths = total if total > 0 else max(len(chrom_sizes), len(blocks))

    # Inch layout: a centred tile strip, then panel a with the legend strip docked under its
    # x-label; b holds ten bar slots and ends level with the strip, so no band is left empty.
    fig = plt.figure(figsize=SLIDE_SIZE)
    left, right = 0.145 * FIG_WIDTH_DOUBLE, 0.97 * FIG_WIDTH_DOUBLE
    top, tile_w = 1.05, 1.12
    b_left = 3.84 + PANEL_LETTER_DX_IN  # letter b where the terminal page puts it
    ax_tiles = _add_axes_in(fig, (FIG_WIDTH_DOUBLE - 5 * tile_w) / 2, 0.36, 5 * tile_w, 0.34)
    ax_summary = _add_axes_in(fig, left, top, 3.55 - left, 1.35)
    ax_flagged_scaffolds = _add_axes_in(fig, b_left, top, right - b_left, 2.10)
    ax_legends = _add_axes_in(fig, left, 2.72, 3.55 - left, 0.43)  # legends hang from its top

    _draw_stat_tiles(ax_tiles, [
        (f"{_median(all_len):,.0f}" if all_len else "NA", "Median assembled array (bp)"),
        (f"{total_paths:,}", "Scaffolds"),
        (f"{len(flagged_scaffolds):,}", "Flagged scaffolds"),
        (f"{len(block_rows):,}", "Terminal telomeres"),
        (f"{sum(r['distance'] > FLAGGED_DISTANCE_BP for r in block_rows):,}",
         "Distance-flagged telomeres"),
    ])

    # ---- Panel a: Classification bars ----
    absolute_segments = [
        (lab, cnt, color)
        for lab, cnt, color in zip(cat_labels, cat_counts, cat_colors)
        if cnt > 0
    ]
    excluded = {"No telomeres", "Gapped No telomeres"}
    telomeric_segments = [
        (lab, cnt, color)
        for lab, cnt, color in absolute_segments
        if lab not in excluded
    ]
    telomeric_total = sum(cnt for _, cnt, _ in telomeric_segments)
    relative_segments = [
        (lab, (100.0 * cnt / telomeric_total) if telomeric_total else 0.0, color, cnt)
        for lab, cnt, color in telomeric_segments
    ]
    absolute_pct_segments = [
        (lab, (100.0 * cnt / total) if total else 0.0, color, cnt)
        for lab, cnt, color in absolute_segments
    ]
    _draw_dual_summary_panel(ax_summary, absolute_pct_segments, relative_segments)
    _panel_title(ax_summary, "Scaffold classification", "a", x_in=PANEL_LETTER_X_IN)

    # ---- Panel b: Flagged scaffolds (discordant / fragmented by size) ----
    _draw_flagged_scaffolds_panel(ax_flagged_scaffolds, classifications, chrom_sizes)
    _panel_title(ax_flagged_scaffolds, "Flagged scaffolds", "b")

    # ---- Shared legend strip: gap key + quality legend ----
    quality_handles = [
        Patch(facecolor=COLORS["T2T"],         label="T2T"),
        Patch(facecolor=COLORS["Incomplete"],   label="Incomplete"),
        Patch(facecolor=COLORS["Fragmented"],   label="Fragmented"),
        Patch(facecolor=COLORS["Discordant"],   label="Discordant"),
        Patch(facecolor=COLORS["No telomeres"], label="No telomeres"),
    ]
    gap_handles = [
        Patch(facecolor="#555555", label="Gapless"),
        Patch(facecolor="#cccccc", label="Gapped"),
    ]
    _draw_centered_legend_pair(
        ax_legends, gap_handles, quality_handles, gap=0.016,
        left_legend_kwargs={"ncol": 1, "fontsize": LEGEND_TEXT_SIZE, "loc": "upper center",
                            "bbox_to_anchor": (0.5, 1.0)},
        right_legend_kwargs={"ncol": 3, "fontsize": LEGEND_TEXT_SIZE, "loc": "upper center",
                             "bbox_to_anchor": (0.5, 1.0)},
    )

    _page_title(fig, "Assembly telomere summary")
    return fig


def plot_overview_page2(blocks, chrom_sizes):
    """Page 2: Length by arm (c), telomere positioning (d), flagged telomeres (e)."""
    FLAGGED_DISTANCE_BP = 1000
    FLAGGED_TOP_N = 10

    block_rows = _compute_block_rows(blocks, chrom_sizes)
    flagged_rows = [r for r in block_rows if r["distance"] > FLAGGED_DISTANCE_BP]
    flagged_by_scaffold = defaultdict(lambda: {"length": 0, "count": 0, "arms": Counter()})
    for fr in flagged_rows:
        entry = flagged_by_scaffold[fr["chrom"]]
        entry["length"] += fr["length"]
        entry["count"] += 1
        entry["arms"][fr["closestEnd"]] += 1
    flagged_ranked = sorted(flagged_by_scaffold.items(),
                            key=lambda kv: kv[1]["length"], reverse=True)
    flagged_top = flagged_ranked[:FLAGGED_TOP_N]

    p_len = [row["length"] for row in block_rows if row["closestEnd"] == "p"]
    q_len = [row["length"] for row in block_rows if row["closestEnd"] == "q"]
    b_len = [row["length"] for row in block_rows if row["label"] == "b"]
    all_len = p_len + q_len + b_len

    # Inch layout: three square axes sharing one top edge, so titles share a baseline.
    fig = plt.figure(figsize=SLIDE_SIZE)
    top, side = 1.05, 1.66  # plots sit below the independent panel headings
    right = 0.97 * FIG_WIDTH_DOUBLE
    ax_rain = _add_axes_in(fig, 0.66, top, side, side)
    ax_scatter = _add_axes_in(fig, 2.98, top, side, side)
    ax_flagged = _add_axes_in(fig, 5.21, top, right - 5.21, side)
    # c and d share one linear kbp length axis.
    y_top = max(max(all_len, default=0) / 1e3 * 1.06, 1.2)
    for ax in (ax_rain, ax_scatter):
        ax.set_ylim(-0.03 * y_top, y_top)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    tick_style = dict(length=2.3, width=0.45, pad=1.5, labelsize=AXIS_TICK_SIZE)
    for ax in (ax_rain, ax_scatter):
        for spine in ("left", "bottom"):
            ax.spines[spine].set_linewidth(0.45)
            ax.spines[spine].set_color("#bcbcbc")
        ax.tick_params(axis="both", **tick_style)

    # ---- Panel c: Vertical raincloud ----
    groups = [
        ("p-arm", p_len, COLORS["p"]),
        ("q-arm", q_len, COLORS["q"]),
        ("balanced", b_len, COLORS["b"]),
    ]
    groups = [(lbl, vals, col) for lbl, vals, col in groups if vals]
    ax_rain.set_ylabel("Telomere length (kbp)", fontsize=AXIS_LABEL_SIZE)

    if groups:
        rng = np.random.default_rng(0)
        spacing = 0.36
        offsets = np.arange(len(groups), dtype=np.float64)
        offsets -= (len(groups) - 1) / 2.0
        positions = 0.75 + (spacing * offsets)
        for position, (lbl, vals, col) in zip(positions, groups):
            _draw_raincloud_group(ax_rain, np.asarray(vals, dtype=np.float64) / 1e3, position, col, rng, vert=True)

        x_labels = [f"{lbl}\n(n={len(vals)})" for lbl, vals, _ in groups]
        ax_rain.set_xticks(positions)
        ax_rain.set_xticklabels(x_labels, fontsize=AXIS_TICK_SIZE)
        ax_rain.set_xlim(positions.min() - 0.34, positions.max() + 0.30)

        median_kbp = _median(all_len) / 1e3
        ax_rain.axhline(median_kbp, color="#555555",
                        linestyle=OVERVIEW_DASH_STYLE, linewidth=OVERVIEW_DASH_WIDTH)
        median_trans = transforms.blended_transform_factory(ax_rain.transAxes, ax_rain.transData)
        ax_rain.text(1.01, median_kbp, "median",
                     transform=median_trans, ha="left", va="center",
                     rotation=90, fontsize=LEGEND_TEXT_SIZE, color="#555555")
    else:
        _draw_empty_state(ax_rain, "No labeled terminal telomeres" if block_rows else "No terminal telomeres")

    # ---- Panel d: Scatter (terminal offset vs length) ----
    scatter_groups = [
        ("p-arm", "p", COLORS["p"]),
        ("q-arm", "q", COLORS["q"]),
        ("balanced", "b", COLORS["b"]),
    ]
    scatter_rows = [row for row in block_rows if row["label"] in {"p", "q", "b"}]
    ax_scatter.set_xlabel("Distance to end (bp)", fontsize=AXIS_LABEL_SIZE, labelpad=2.5)
    ax_scatter.set_ylabel("Telomere length (kbp)", fontsize=AXIS_LABEL_SIZE)
    if scatter_rows:
        max_x = 1.0
        for leg_lbl, arm_key, col in scatter_groups:
            field = "label" if arm_key == "b" else "closestEnd"
            arm_rows = [row for row in scatter_rows if row[field] == arm_key]
            if not arm_rows:
                continue
            # Blocks touching the end sit on the 1 bp tick of the log axis.
            x = np.maximum(np.asarray([row["distance"] for row in arm_rows], dtype=np.float64), 1.0)
            y = np.asarray([row["length"] / 1e3 for row in arm_rows], dtype=np.float64)
            max_x = max(max_x, float(np.max(x)))
            ax_scatter.scatter(x, y, s=11, color=col, alpha=0.55,
                               edgecolors="white", linewidths=0.28,
                               rasterized=True, label=leg_lbl)

        x_right = max(10 ** 4.12, max_x * 10 ** 0.18)
        ax_scatter.set_xscale("log")
        ax_scatter.axvspan(FLAGGED_DISTANCE_BP, x_right, facecolor="#ededed", linewidth=0, zorder=0)
        ax_scatter.set_xlim(10 ** -0.10, x_right)
        decades = [10 ** k for k in range(int(np.log10(x_right)) + 1)]
        ax_scatter.xaxis.set_major_locator(ticker.FixedLocator(decades))
        ax_scatter.xaxis.set_major_formatter(ticker.FixedFormatter([f"{d:d}" if d <= 10_000 else "" for d in decades]))
        ax_scatter.xaxis.set_minor_locator(ticker.NullLocator())
        leg_d = ax_scatter.legend(loc="upper right", fontsize=LEGEND_TEXT_SIZE,
                                  frameon=True, facecolor="white", edgecolor="black",
                                  framealpha=1.0, fancybox=False, handletextpad=0.35,
                                  borderaxespad=0.3, markerscale=0.8)
        leg_d.get_frame().set_linewidth(0.4)
        ax_scatter.axvline(FLAGGED_DISTANCE_BP, color="#555555", linewidth=OVERVIEW_DASH_WIDTH,
                           linestyle=OVERVIEW_DASH_STYLE, zorder=1)
    else:
        _draw_empty_state(ax_scatter, "No terminal telomeres")

    # ---- Panel e: Flagged telomeres ----
    if not block_rows:
        _draw_empty_state(ax_flagged, "No terminal telomeres")
    elif not flagged_top:
        _draw_empty_state(ax_flagged, f"No flagged telomeres\n(distance > {FLAGGED_DISTANCE_BP // 1000} kbp)")
    else:
        rows = []
        for scaffold, info in flagged_top:
            dominant_arm = info["arms"].most_common(1)[0][0]
            rows.append({
                "label": f"{_wrap_scaffold_name(scaffold, width=14)} (n={info['count']})",
                "value": info["length"] / 1e3,
                "color": COLORS.get(dominant_arm, "#aaaaaa"),
            })
        _draw_ranked_bar_panel(ax_flagged, rows, x_label="Flagged length (kbp)",
                               xticks=list(ticker.MaxNLocator(nbins=4).tick_values(0, rows[0]["value"])))

    _panel_title(ax_rain, "Length by arm", "c", x_in=PANEL_LETTER_X_IN)
    _panel_title(ax_scatter, "Telomere positioning", "d", x_in=2.48)
    _panel_title(ax_flagged, "Flagged telomeres", "e", x_in=4.71)
    for ax in (ax_rain, ax_scatter, ax_flagged):
        ax.set_title(ax.get_title(), pad=24.6)  # keep headings 0.3 in above the plots
        letter = ax.texts[-1]
        letter.set_transform(letter.get_transform() + transforms.ScaledTranslation(0, 0.3, fig.dpi_scale_trans))
    _page_title(fig, "Terminal telomeres")
    return fig


# ---------------------------------------------------------------------------
# ITS report pages (rendered last, after every per-chromosome zoom; see main())
# ---------------------------------------------------------------------------

def _its_atlas_chroms(df, arm_blocks, chrom_sizes):
    """Scaffolds carrying a telomere or ITS row, longest first (ideogram row order)."""
    telomere_chroms = {c for c, blist in arm_blocks.items() if blist}
    its_chroms = set(df["chr"].unique()) if not df.empty else set()
    return sorted(telomere_chroms | its_chroms, key=lambda c: (-chrom_sizes.get(c, 0), c))


def _atlas_tiers(atlas_chroms, chrom_sizes):
    """Size tiers of whole homolog groups: the fewest spanning ATLAS_TIER_RATIO each, cut at the widest size gaps."""
    groups = _homolog_groups(atlas_chroms, chrom_sizes)
    tops = np.log([max(max(chrom_sizes.get(c, 0) for c in g), 1) for g in groups])
    n = len(groups)
    gap = np.append(tops[:-1] - tops[1:], 0.0)
    count, score, cut = np.zeros(n + 1), np.zeros(n + 1), np.full(n + 1, n)
    for i in range(n - 1, -1, -1):  # best split of groups[i:], built from the smallest groups up
        ends = np.arange(i + 1, np.searchsorted(-tops, np.log(ATLAS_TIER_RATIO) - tops[i], side="right") + 1)
        j = ends[np.argmin((count[ends] + 1) * 1e6 - score[ends] - gap[ends - 1])]
        count[i], score[i], cut[i] = count[j] + 1, score[j] + gap[j - 1], j
    tiers, i = [], 0
    while i < n:
        tiers.append([c for g in groups[i:cut[i]] for c in g])
        i = cut[i]
    return tiers


def _atlas_extent(chroms, chrom_sizes, cells):
    """Drawn length per scaffold: its size, or its last ITS end when the size is unknown."""
    return {c: max(int(chrom_sizes.get(c, 0)), int(cells[c][1].max()) if c in cells else 0, 1) for c in chroms}


def _nearest_end(pos, chrom_size):
    """Return ('p'|'q', distance) for whichever scaffold end pos sits closest to."""
    dist_p, dist_q = pos, max(chrom_size - pos, 0)
    return ("p", dist_p) if dist_p <= dist_q else ("q", dist_q)


# ---------------------------------------------------------------------------
# Region-zoom toolkit shared with plot_its.py (single-locus figures)
# ---------------------------------------------------------------------------

def _region_unit(start, end):
    """Locus-page unit: Mbp for spans from 5 Mb or windows past 10 Mb, kbp from 1 kb, else bp."""
    span = end - start
    return (1e6, "Mbp") if span >= 5_000_000 or end >= 10_000_000 else (1e3, "kbp") if span >= 1_000 else (1.0, "bp")


def _region_axis_spec(view_start, view_end, plain=False):
    """Round genomic ticks for an absolute-position window; in Mbp, or with plain=True in the span's unit."""
    if plain:
        div, unit = _region_unit(view_start, view_end)
        ticks = [t * div for t in ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(
            view_start / div, view_end / div) if view_start <= t * div <= view_end]
        return {"axis_start": view_start, "axis_end": view_end, "tick_pos": ticks,
                "tick_labels": [_plain_number(t / div) for t in ticks], "xlabel": f"Position ({unit})"}
    span_mb = max(view_end - view_start, 1) / 1e6
    dec = _mb_tick_decimals(span_mb)
    ticks = [t for t in ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(
        view_start, view_end) if view_start <= t <= view_end]
    if len(ticks) < 2:
        ticks = list(np.linspace(view_start, view_end, 5))
    return {
        "axis_start": view_start, "axis_end": view_end, "tick_pos": ticks,
        "tick_labels": [f"{t / 1e6:.{dec}f}" for t in ticks],
        "xlabel": "Position (Mbp)",
    }


def _block_key_color(block, threshold):
    """Double-key colour of one parsed BED block from its four match counts."""
    strand, canon = composition_class(*(block.get(k, 0) for k in ("fwdCan", "revCan", "fwdNonCan", "revNonCan")),
                                      threshold)
    return DOUBLE_KEY_COLORS[int(canon)][int(strand)]


def _draw_locus_ideogram(ax, chrom_size, view_start, view_end, blocks_list, its_blocks_list,
                         contig_blocks_list=None, label="Scaffold", threshold=None):
    """Whole-scaffold overview for a single-locus zoom: terminal caps, all ITS ticks, zoom box.

    threshold: when set, ITS ticks take their double-key colours instead of p/q/b.
    """
    bar_y, bar_h = 0.60, 0.34
    bottom = bar_y - bar_h / 2.0
    ax.add_patch(Rectangle((0, bottom), chrom_size, bar_h, facecolor="#ececec",
                           edgecolor="#b9b9b9", linewidth=0.5, zorder=1))

    cap_min = chrom_size * 0.004
    for b in blocks_list or []:
        if b.get("term") == "contig":
            continue
        w = max(b["end"] - b["start"], cap_min)
        x = min(b["start"], chrom_size - w) if b["start"] > chrom_size / 2.0 else b["start"]
        ax.add_patch(Rectangle((x, bottom), w, bar_h, facecolor=COLORS["terminal"],
                               edgecolor="none", zorder=2))

    for b in contig_blocks_list or []:
        w = max(b["end"] - b["start"], cap_min)
        x = min(b["start"], chrom_size - w) if b["start"] > chrom_size / 2.0 else b["start"]
        ax.add_patch(Rectangle((x, bottom), w, bar_h, facecolor="none",
                               edgecolor=COLORS["terminal"], linewidth=0.6, zorder=3))

    if its_blocks_list:
        # Vectorised: a chromosome can carry thousands of ITS rows (fast-mode genomes),
        # and one ax.plot() call per row was the dominant cost of this figure.
        xs = np.array([(b["start"] + b["end"]) / 2.0 for b in its_blocks_list])
        labels = np.array([b.get("label", "") if threshold is None else _block_key_color(b, threshold)
                           for b in its_blocks_list])
        for k in np.unique(labels):
            mask = labels == k
            segs = np.empty((int(mask.sum()), 2, 2))
            segs[:, 0, 0] = xs[mask]; segs[:, 0, 1] = bottom
            segs[:, 1, 0] = xs[mask]; segs[:, 1, 1] = bottom + bar_h
            ax.add_collection(LineCollection(
                segs, colors=ITS_ORIENT_COLORS.get(k, COLORS["its"]) if threshold is None else k,
                linewidths=0.7, zorder=3, rasterized=True))

    mid = (view_start + view_end) / 2.0
    half = max((view_end - view_start) / 2.0, chrom_size * 0.0035)
    box_l, box_r = max(0, mid - half), min(chrom_size, mid + half)
    box_bottom = bottom - 0.12
    ax.add_patch(Rectangle((box_l, box_bottom), box_r - box_l, bar_h + 0.24,
                           facecolor="none", edgecolor=ZOOM_COLOR, linewidth=1.0, zorder=5))

    ax.set_xlim(-chrom_size * 0.012, chrom_size * 1.012)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ticks = [t for t in ticker.MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(
        0, chrom_size) if 0 <= t <= chrom_size]
    dec = _mb_tick_decimals(chrom_size / 1e6)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t / 1e6:.{dec}f}" for t in ticks])
    ax.tick_params(axis="x", length=2, width=0.4, pad=1.0,
                   labelsize=AXIS_TICK_SIZE, colors="black")
    ax.set_xlabel("Scaffold (Mbp)", fontsize=AXIS_TICK_SIZE, color="black", labelpad=1.5)
    _set_track_label(ax, label)
    return box_l, box_r, box_bottom


def _draw_its_blocks_track(ax, view_start, view_end, blocks_list, its_blocks_list,
                           gap_blocks_list, label="Blocks", contig_blocks_list=None, threshold=None):
    """Single-panel block track in genomic coordinates; ITS colored by orientation.

    threshold: when set, ITS take their double-key colours instead of p/q/b.

    Rectangles are batched into one PatchCollection per group (not one add_patch per
    block), since a locus window in a dense fast-mode genome can hold thousands of rows.
    """
    backbone_y = 0.5
    view_span = max(view_end - view_start, 1)
    ax.plot([view_start, view_end], [backbone_y, backbone_y],
            color="#d6d6d6", linewidth=0.9, solid_capstyle="round", zorder=1)

    def _clip(seq):
        clipped = []
        for b in seq or []:
            if b["end"] <= view_start or b["start"] >= view_end:
                continue
            ds, de = max(b["start"], view_start), min(b["end"], view_end)
            clipped.append((b, ds, de))
        return clipped

    def _draw(seq, color_fn, zorder):
        clipped = _clip(seq)
        if not clipped:
            return
        rects = [Rectangle((ds, backbone_y - 0.07), de - ds, 0.14) for _, ds, de in clipped]
        ax.add_collection(PatchCollection(
            rects, facecolor=[color_fn(b) for b, _, _ in clipped], edgecolor="none",
            alpha=0.98, zorder=zorder, rasterized=len(clipped) > 200))
        for b, ds, de in clipped:
            if _block_symbol_fits(de - ds, view_span, b.get("label", "")):
                _draw_block_symbol(ax, (ds + de) / 2, backbone_y, b.get("label", ""), zorder + 1)

    def _draw_outline(seq, color, zorder):
        clipped = _clip(seq)
        if not clipped:
            return
        rects = [Rectangle((ds, backbone_y - 0.07), de - ds, 0.14) for _, ds, de in clipped]
        ax.add_collection(PatchCollection(
            rects, facecolor="none", edgecolor=color, linewidth=0.6, zorder=zorder,
            rasterized=len(clipped) > 200))

    _draw([b for b in blocks_list if b.get("term") != "contig"] if blocks_list else [],
          lambda b: COLORS["terminal"], 2)
    _draw_outline([b for b in blocks_list if b.get("term") == "contig"] if blocks_list else [],
                  COLORS["terminal"], 3)
    _draw(its_blocks_list, (lambda b: ITS_ORIENT_COLORS.get(b.get("label", ""), COLORS["its"]))
          if threshold is None else (lambda b: _block_key_color(b, threshold)), 4)
    _draw(gap_blocks_list, lambda b: COLORS["gap"], 6)
    _draw_outline(contig_blocks_list, COLORS["terminal"], 5)

    ax.set_ylim(0.28, 0.72)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    _set_track_label(ax, label)
    _hide_x_axis(ax)
    return True


# ---------------------------------------------------------------------------
# ITS-1 "Interstitial telomeres: genome view"
# ---------------------------------------------------------------------------

def _its_labels(df, chrom_sizes, arm_blocks=None, names=None):
    names = _its_atlas_chroms(df, arm_blocks or {}, chrom_sizes) if names is None else names
    # Stable aliases preserve identity even when a long suffix is also shared.
    labels = {}
    reserved = set(names)
    used = set()
    for i, name in enumerate(names):
        label = name if len(name) <= 18 else f"S{i + 1}:…{name[-10:]}"
        while label != name and (label in reserved or label in used):
            label = "_" + label
        labels[name] = label
        used.add(label)
    return labels


def its_scaffold_summary(df, arm_blocks, chrom_sizes):
    """Complete observed-row accounting, including zero-ITS terminal scaffolds."""
    names = _its_atlas_chroms(df, arm_blocks, chrom_sizes)
    labels = _its_labels(df, chrom_sizes, arm_blocks, names=names)
    grouped = df.groupby("chr").agg(rows=("chr", "size"), its_bp=("teloLen", "sum"),
                                    canonical_bp=("canonical_bp", "sum"), _size=("chrSize", "max"))
    out = grouped.reindex(names, fill_value=0).rename_axis("chr").reset_index()
    out.insert(1, "display_label", [labels[c] for c in names])
    out.insert(2, "plot_extent_bp", [chrom_sizes.get(c, 0) for c in names])
    out.insert(3, "its_size_known", out.pop("_size") > 0)
    return out


def _its_atlas_cells(df, threshold=LABEL_THRESHOLD):
    """Per scaffold: ITS starts, ends and double-key colours, widest first so narrow ITS stay on top."""
    if df.empty:
        return {}
    colors = np.array(double_key_colors(df, threshold), dtype=object)
    starts, ends = df["start"].to_numpy(), df["end"].to_numpy()
    out = {}
    for chrom, idx in df.groupby("chr", sort=False).indices.items():
        order = idx[np.argsort(-(ends[idx] - starts[idx]), kind="stable")]
        out[chrom] = (starts[order], ends[order], colors[order])
    return out


def _its_top_hits(df, pairs, clusters):
    """Scaffold -> [(ranked tag, midpoint)] for up to five loci per category."""
    top = {}
    for prefix, frame in (("C", clusters.head(5)), ("F", pairs.head(5)), ("L", rank_long_its(df, 5))):
        for rank, row in enumerate(frame.itertuples(index=False), 1):
            top.setdefault(row.chr, []).append((f"{prefix}{rank}", (row.start + row.end) / 2.0))
    return top


_HAPLOTYPE_TOKEN = re.compile(
    r"(?i)(?:^|[_.\-|])(?:mat(?:ernal)?|pat(?:ernal)?|hap(?:lotype)?[_.\-]?[12]|h[12])(?=$|[_.\-|])")


def _homolog_key(name):
    """Scaffold name without its haplotype token (mat/pat, hap1/hap2, h1/h2, PanSN #1#/#2#)."""
    return _HAPLOTYPE_TOKEN.sub("", re.sub(r"#[12]#", "#", name)).strip("_.-|") or name


def _homolog_groups(chroms, chrom_sizes):
    """Homolog groups, largest group first; members in name order so haplotypes keep one order."""
    groups = OrderedDict()
    for chrom in chroms:
        groups.setdefault(_homolog_key(chrom), []).append(chrom)
    return [sorted(g) for g in sorted(groups.values(),
                                      key=lambda g: (-max(chrom_sizes.get(c, 0) for c in g), min(g)))]


def _set_mb_ticks(ax, longest, unit=1e6):
    """Plain ticks from 0 to a scaffold length in Mbp (or another unit), whole steps from 2 units up."""
    ticks = ticker.MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10], integer=longest >= 2 * unit).tick_values(0, longest / unit)
    ax.set_xticks([t * unit for t in ticks if 0 <= t * unit <= longest])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _plain_number(x / unit)))


def _atlas_rows(names):
    """Row centres in inches below a segment top, homolog groups spaced apart, and the segment height."""
    ys, y = [], -ATLAS_PITCH_IN / 2.0
    for i, c in enumerate(names):
        y += ATLAS_PITCH_IN + (ATLAS_GROUP_GAP_IN if i and _homolog_key(c) != _homolog_key(names[i - 1]) else 0.0)
        ys.append(y)
    return ys, y + ATLAS_PITCH_IN / 2.0


def _atlas_columns(chroms, sizes, col_h, tiers=None):
    """Flow homolog groups down columns of col_h inches; each column holds one axis segment per tier it touches."""
    tiers = _atlas_tiers(chroms, sizes) if tiers is None else tiers
    page = set(chroms)
    head = ATLAS_AXIS_IN + (ATLAS_TITLE_IN if len(tiers) > 1 else 0.0)
    columns, used = [[]], 0.0
    for k, tier in enumerate(tiers):
        for group in _homolog_groups([c for c in tier if c in page], sizes):
            same = bool(columns[-1]) and columns[-1][-1][0] == k
            cost = len(group) * ATLAS_PITCH_IN + (ATLAS_GROUP_GAP_IN if same else head)
            if columns[-1] and used + cost > col_h + 1e-9:
                columns.append([])
                used, same = 0.0, False
                cost = len(group) * ATLAS_PITCH_IN + head
            if not same:
                columns[-1].append((k, []))
            columns[-1][-1][1].extend(group)
            used += cost
    return columns, len(tiers) > 1


def _atlas_layout(chroms, sizes, tiers=None):
    """Page height, balanced columns and split flag: one full-width column, else two, on one slide, else a double slide."""
    tiers = _atlas_tiers(chroms, sizes) if tiers is None else tiers
    for page_h in ATLAS_HEIGHTS:
        col_h = page_h - ATLAS_HEADER_IN - ATLAS_KEY_IN - ATLAS_FOOT_IN
        for n in range(1, ATLAS_COLUMNS + 1):
            columns, split = _atlas_columns(chroms, sizes, col_h, tiers)
            if len(columns) > n:
                continue
            lo, hi = 0.0, col_h  # shortest column height that still needs only n columns
            for _ in range(24):
                mid = (lo + hi) / 2.0
                lo, hi = (lo, mid) if len(_atlas_columns(chroms, sizes, mid, tiers)[0]) <= n else (mid, hi)
            return page_h, _atlas_columns(chroms, sizes, hi, tiers)[0], split
    return page_h, columns, split


def _paginate_atlas(atlas_chroms, chrom_sizes):
    """Atlas pages holding as many scaffolds as fit a double slide at the fixed row pitch; homolog groups stay whole."""
    col_h = ATLAS_HEIGHTS[-1] - ATLAS_HEADER_IN - ATLAS_KEY_IN - ATLAS_FOOT_IN
    tiers = _atlas_tiers(atlas_chroms, chrom_sizes)
    pages, page = [], []
    for group in _homolog_groups(atlas_chroms, chrom_sizes):
        if page and len(_atlas_columns(page + group, chrom_sizes, col_h, tiers)[0]) > ATLAS_COLUMNS:
            pages.append(page)
            page = []
        page.extend(group)
    return pages + [page] if page or not pages else pages


def plot_its_overview_page(df, pairs, clusters, arm_blocks, chrom_sizes, params,
                           chroms=None, page_number=1, page_count=1, cells=None,
                           all_chroms=None, top_hits=None, known_sizes=None):
    """ITS atlas: every ITS at its true position on a shared Mb axis, coloured by its double-key class."""
    all_chroms = _its_atlas_chroms(df, arm_blocks, chrom_sizes) if all_chroms is None else all_chroms
    chroms = _paginate_atlas(all_chroms, chrom_sizes)[0] if chroms is None else list(chroms)
    threshold = params.get("label_threshold", LABEL_THRESHOLD)
    cells = _its_atlas_cells(df, threshold) if cells is None else cells
    labels = df.attrs.get("display_labels") or _its_labels(df, chrom_sizes, arm_blocks)
    top = _its_top_hits(df, pairs, clusters) if top_hits is None else top_hits
    known = (df.groupby("chr")["chrSize"].max().to_dict() if known_sizes is None else known_sizes)
    limit = params.get("terminal_limit") or 0
    end_scan = params.get("ultra_fast") is True and limit > 0
    extent = _atlas_extent(list(dict.fromkeys([*all_chroms, *chroms])), chrom_sizes, cells)
    tiers = _atlas_tiers(list(extent), extent)  # dataset-wide, so a tier keeps one scale on every page
    height, columns, split = _atlas_layout(chroms, extent, tiers)
    width = FIG_WIDTH_DOUBLE
    fig = plt.figure(figsize=(width, height))
    page_suffix = f" ({page_number}/{page_count})" if page_count > 1 else ""
    _page_title(fig, "ITS atlas" + page_suffix)
    if not chroms:
        fig.text(0.5, 0.45, "No ITS", ha="center", va="center",
                 fontsize=PLACEHOLDER_TEXT_SIZE, color="#bbbbbb")
        return fig

    # One-line row labels; the widest sets the label column so every bar starts at the same x.
    row_label = {c: f"{labels.get(c, c)}  {'≥' if known.get(c) == 0 else ''}{_fmt_bp_fixed(extent[c])}"
                 for c in chroms}
    font = FontProperties(size=AXIS_TICK_SIZE)
    label_w = max(TextToPath().get_text_width_height_descent(t, font, ismath=False)[0]
                  for t in row_label.values()) / 72 + 0.10
    margin = 0.15
    col_w = (width - 2 * margin - ATLAS_COL_GAP_IN * (len(columns) - 1)) / len(columns)
    axes_w = col_w - label_w
    longest = [max(extent[c] for c in tier) for tier in tiers]  # one scale per tier, shared across columns and pages
    lift = transforms.ScaledTranslation(0, 0.6 / 72, fig.dpi_scale_trans)  # triangle tip just clears the bar
    drew_caps, drew_marks = False, set()
    used = max(sum(_atlas_rows(names)[1] + ATLAS_AXIS_IN + (ATLAS_TITLE_IN if split else 0.0) for _, names in column)
               for column in columns)
    drop = max(height - ATLAS_HEADER_IN - ATLAS_KEY_IN - ATLAS_FOOT_IN - used, 0.0) / 2.0  # key and rows centred as one block
    for n, column in enumerate(columns):
        x_in = margin + n * (col_w + ATLAS_COL_GAP_IN) + label_w
        top_in = height - ATLAS_HEADER_IN - ATLAS_KEY_IN - drop
        for k, names in column:
            top_in -= ATLAS_TITLE_IN if split else 0.0
            ys, seg_h = _atlas_rows(names)
            ax = _inch_axes(fig, x_in, top_in - seg_h, axes_w, seg_h)
            top_in -= seg_h + ATLAS_AXIS_IN
            span = longest[k]
            min_w = span * (1.2 / 72) / axes_w  # 1.2 pt floor keeps single ITS visible
            backbone, masks, caps, its_rects, its_colors = [], [], [], [], []
            for y, c in zip(ys, names):
                size = extent[c]
                y0 = y - ATLAS_BAR_IN / 2.0
                backbone.append(Rectangle((0, y0), size, ATLAS_BAR_IN))
                if end_scan and known.get(c, size) > 0 and size > 2 * limit:
                    masks.append(Rectangle((limit, y0), size - 2 * limit, ATLAS_BAR_IN))
                if c in cells:
                    for s, e, color in zip(*cells[c]):
                        w = max(e - s, min_w)
                        its_rects.append(Rectangle((min(s, size - w), y0), w, ATLAS_BAR_IN))
                        its_colors.append(color)
                for block in arm_blocks.get(c, []):
                    w = max(block["end"] - block["start"], min_w)
                    x = min(block["start"], size - w) if block["start"] > size / 2.0 else block["start"]
                    caps.append(Rectangle((max(x, 0), y0), w, ATLAS_BAR_IN))
                marks = []
                for tag, x in sorted(top.get(c, []), key=lambda m: m[1]):
                    if marks and (x - marks[-1][1]) / span * axes_w < max(0.22, 0.05 * len(marks[-1][0]) + 0.08):
                        marks[-1] = (marks[-1][0] + "/" + tag, marks[-1][1])
                    else:
                        marks.append((tag, x))
                for tag, x in marks:
                    tag = "/".join(sorted(tag.split("/"), key=lambda t: ("CFL".find(t[0]), int(t[1:]))))
                    drew_marks.update(tag.split("/"))
                    ax.plot([x], [y0], marker="v", markersize=3, color="#222222", markeredgewidth=0,
                            linestyle="none", transform=ax.transData + lift, clip_on=False, zorder=6)
                    ax.annotate(tag, (x, y0), xytext=(2.0, 1.6), textcoords="offset points", ha="left",
                                va="center", fontsize=MIN_TEXT_SIZE, fontweight="bold", color="#222222",
                                annotation_clip=False, zorder=6,
                                path_effects=[patheffects.withStroke(linewidth=1.2, foreground="white")])
            ax.add_collection(PatchCollection(backbone, facecolor=ATLAS_BACKBONE, edgecolor="none", zorder=1))
            if masks:
                ax.add_collection(PatchCollection(masks, facecolor="white", edgecolor=COLORS["gap"],
                                                  linewidth=0.4, zorder=2))
            if its_rects:
                ax.add_collection(PatchCollection(its_rects, facecolor=its_colors, edgecolor="none",
                                                  zorder=3, rasterized=len(its_rects) > 2000))
            if caps:
                drew_caps = True
                ax.add_collection(PatchCollection(caps, facecolor=COLORS["terminal"], edgecolor="none",
                                                  zorder=4))
            # Thin left bracket joins the members of a homolog group.
            bracket = transforms.blended_transform_factory(ax.transAxes, ax.transData)
            i = 0
            while i < len(names):
                j = i
                while j + 1 < len(names) and _homolog_key(names[j + 1]) == _homolog_key(names[i]):
                    j += 1
                if j > i:
                    ax.plot([-0.008, -0.008], [ys[i] - ATLAS_BAR_IN / 2, ys[j] + ATLAS_BAR_IN / 2],
                            color="#8C8C8C", lw=0.5, transform=bracket, clip_on=False, solid_capstyle="butt")
                i = j + 1
            ax.set_yticks(ys)
            ax.set_yticklabels([row_label[c] for c in names], fontsize=AXIS_TICK_SIZE)
            ax.tick_params(axis="y", length=0, pad=4)
            ax.set_ylim(seg_h, 0)
            ax.set_xlim(0, span)
            unit, name = _pick_bp_unit(max(span, 1e3))
            _set_mb_ticks(ax, span, unit)
            ax.tick_params(axis="x", length=2, width=0.4, pad=1.0, labelsize=AXIS_TICK_SIZE)
            ax.set_xlabel(f"Position ({name}p)", fontsize=AXIS_TICK_SIZE, labelpad=1.5)
            for side in ("left", "right", "top"):
                ax.spines[side].set_visible(False)
            ax.spines["bottom"].set_linewidth(0.45)
            if split:  # tier range in the row-label format, the unit written once when both ends share it
                lo, hi = _fmt_bp_fixed(min(extent[c] for c in tiers[k])), _fmt_bp_fixed(span)
                _panel_title(ax, f"Scaffolds {lo.rsplit(' ', 1)[0] if lo[-2:] == hi[-2:] else lo}–{hi}")

    # Centre the rows, from labels to the last tick, on the page.
    renderer = fig.canvas.get_renderer()
    boxes = [ax.get_tightbbox(renderer) for ax in fig.axes]
    dx = (width - (max(b.x1 for b in boxes) + min(b.x0 for b in boxes)) / fig.dpi) / 2.0 / width
    for ax in fig.axes:
        pos = ax.get_position()
        ax.set_position([pos.x0 + dx, pos.y0, pos.width, pos.height])

    handles = [Patch(facecolor=ATLAS_BACKBONE, edgecolor="none", label="Scaffold")]
    if drew_caps:
        handles.append(Patch(facecolor=COLORS["terminal"], edgecolor="none", label="Terminal telomere"))
    if end_scan:
        handles.append(Patch(facecolor="white", edgecolor=COLORS["gap"], linewidth=0.4,
                             label=f"Not scanned (end scan, {_fmt_bp(limit)} per end)"))
    handles += [Line2D([], [], marker="v", markersize=3, color="#222222", markeredgewidth=0, linestyle="none",
                       label=f"({tag})  {name}") for tag, name in (("C", "Top clusters"), ("F", "Top fusions"),
                                                                  ("L", "Longest canonical ITS"))
                if any(mark.startswith(tag) for mark in drew_marks)]
    # Key and legend share one strip under the title, centred as a pair.
    key_side = LOCUS_KEY_SIDE
    key_y = height - ATLAS_HEADER_IN - key_side - 0.04 - drop
    key = _inch_axes(fig, 0, key_y, key_side, key_side)
    _draw_double_key(key)
    legend = fig.legend(handles=handles, loc="center left", frameon=False, fontsize=LEGEND_TEXT_SIZE,
                        ncol=(len(handles) + 1) // 2, bbox_to_anchor=(0, (key_y + key_side / 2) / height),
                        handlelength=1.4, handleheight=0.8, borderaxespad=0, labelspacing=0.6, columnspacing=1.6)
    key_box = key.get_tightbbox(renderer)
    key_w = (key_box.x1 - key_box.x0) / fig.dpi
    legend_w = legend.get_window_extent(renderer).width / fig.dpi
    x0 = (width - key_w - 0.35 - legend_w) / 2.0
    key.set_position([(x0 + key.get_position().x0 * width - key_box.x0 / fig.dpi) / width,
                      key_y / height, key_side / width, key_side / height])
    legend.set_bbox_to_anchor(((x0 + key_w + 0.35) / width, (key_y + key_side / 2) / height), fig.transFigure)
    return fig


def _inch_axes(fig, x, y, w, h, **kwargs):
    """Axes placed in inches from the lower-left corner, so panels keep a fixed physical size."""
    width, height = fig.get_size_inches()
    return fig.add_axes([x / width, y / height, w / width, h / height], **kwargs)


def _its_axis_style(ax, left=True):
    """Terminal-page chrome: light 0.45-pt bottom/left spines, short ticks, grey tick labels."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_linewidth(0.45)
        ax.spines[side].set_color("#bcbcbc")
    ax.spines["left"].set_visible(left)
    ax.tick_params(length=2.0, width=0.45, pad=1.5, labelsize=AXIS_TICK_SIZE,
                   color="#bcbcbc", labelcolor="#555555")
    ax.minorticks_off()


def _its_shares(df):
    """Forward and canonical share of each ITS's repeat counts (NaN when it has none)."""
    total = (df["fwdCan"] + df["revCan"] + df["fwdNonCan"] + df["revNonCan"]).to_numpy(dtype=float)
    total[total == 0] = np.nan
    return ((df["fwdCan"] + df["fwdNonCan"]).to_numpy(dtype=float) / total,
            (df["fwdCan"] + df["revCan"]).to_numpy(dtype=float) / total)


def _its_marker_area(lengths, longest, dense):
    """Marker area in pt^2 proportional to ITS length."""
    return np.maximum((24.0 if dense else 60.0) * lengths / max(longest, 1), 1.0)


def _plain_log_axis(ax, which, lo, hi, margin=1.25):
    """Log x or y axis spanning at least a decade, with plain-number ticks at decades or 1-2-5 steps."""
    lo, hi = lo / margin, hi * margin
    if hi / lo < 10:
        mid = np.sqrt(lo * hi)
        lo, hi = mid / np.sqrt(10), mid * np.sqrt(10)
    k0, k1 = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))
    ticks = [10.0 ** k for k in range(k0, k1 + 1) if lo <= 10.0 ** k <= hi]
    if len(ticks) < 2:
        ticks = [m * 10.0 ** k for k in range(k0, k1 + 1) for m in (1, 2, 5) if lo <= m * 10.0 ** k <= hi]
    getattr(ax, f"set_{which}scale")("log")
    getattr(ax, f"set_{which}lim")(lo, hi)
    axis = getattr(ax, f"{which}axis")
    axis.set_ticks(ticks, [f"{t:.0f}" if t >= 1 else f"{t:g}" for t in ticks])
    axis.set_minor_locator(ticker.NullLocator())


def _draw_its_length_panel(ax, lengths):
    """Log-binned ITS length histogram with the median marked."""
    lengths = lengths[lengths > 0]
    lo, hi = np.log10(lengths.min()), np.log10(lengths.max())
    n_bins = int(np.clip(round((hi - lo) * 10), 1, 40))
    edges = np.logspace(lo - 0.05, hi + 0.05, n_bins + 1)
    counts = np.histogram(lengths, edges)[0]
    ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge", color=ITS_LENGTH_COLOR, edgecolor="white", linewidth=0.3)
    median = _median(lengths)
    ax.axvline(median, color="#222222", lw=0.6, zorder=3)
    right_half = np.log10(median) > (lo + hi) / 2
    ax.text(median, 0.94, f" median\n {_fmt_bp(median)} ", transform=ax.get_xaxis_transform(),
            ha="right" if right_half else "left", va="top", fontsize=LEGEND_TEXT_SIZE, color="#222222")
    ax.set_ylim(0, max(1, counts.max()) * 1.15)
    _plain_log_axis(ax, "x", edges[0], edges[-1], margin=1.0)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, integer=True))
    _its_axis_style(ax)
    ax.set_xlabel("ITS length (bp)", fontsize=AXIS_LABEL_SIZE, labelpad=2)
    ax.set_ylabel("# ITSs", fontsize=AXIS_LABEL_SIZE, labelpad=2)


def _draw_its_scaffold_panel(ax, df, chrom_sizes):
    """ITS count against density, with a descriptive log-space fit and sparse labels."""
    per = df.groupby("chr").agg(n=("teloLen", "size"), bp=("teloLen", "sum"), size=("chrSize", "max"))
    per["size"] = [max(chrom_sizes.get(c, 0), s) for c, s in zip(per.index, per["size"])]
    per = per[(per["size"] > 0) & (per["bp"] > 0)]
    if per.empty:
        _hide_panel(ax, "No scaffold size")
        return
    x, y = per["n"].astype(float), per["bp"] / per["size"] * 1e6  # density already carries the size
    if len(per) >= 3 and x.nunique() >= 2:
        slope, intercept = np.polyfit(np.log10(x), np.log10(y), 1)
        ends = np.array([x.min(), x.max()])
        ax.plot(ends, 10 ** (intercept + slope * np.log10(ends)), color="#777777", lw=0.6, zorder=1)
    # Dot diameter is linear in log10 scaffold size, so dot chromosomes stand apart; largest drawn first.
    lo, hi = np.log10(per["size"].min()), np.log10(per["size"].max())
    diameter = lambda v: 2.5 + 5.0 * ((np.log10(v) - lo) / (hi - lo) if hi > lo else 0.5)
    order = np.argsort(-per["size"].to_numpy(), kind="stable")
    ax.scatter(x.iloc[order], y.iloc[order], s=diameter(per["size"].iloc[order].to_numpy(float)) ** 2,
               facecolor="#9ECAE1", alpha=0.8, edgecolor="#333333", linewidths=0.2, zorder=2)
    nice = [m * 10.0 ** k for k in range(int(np.floor(lo)), int(np.ceil(hi)) + 1) for m in (1, 2, 5)
            if lo - 1e-9 <= np.log10(m) + k <= hi + 1e-9]
    refs = sorted({nice[0], nice[len(nice) // 2], nice[-1]}) if nice else [10 ** hi]
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="none", markersize=diameter(v), markerfacecolor="#9ECAE1",
                              alpha=0.8, markeredgecolor="#333333", markeredgewidth=0.2, label=_fmt_bp(v)) for v in refs],
              title="Scaffold size", title_fontsize=LEGEND_TEXT_SIZE, fontsize=LEGEND_TEXT_SIZE, loc="best",
              frameon=False, handletextpad=0.3, labelspacing=0.5, borderaxespad=0.2, borderpad=0.2)
    _plain_log_axis(ax, "x", x.min(), x.max())
    _plain_log_axis(ax, "y", y.min(), y.max(), margin=1.6)
    _its_axis_style(ax)
    ax.set_xlabel("ITS per scaffold", fontsize=AXIS_LABEL_SIZE, labelpad=2)
    ax.set_ylabel("ITS density (bp/Mbp)", fontsize=AXIS_LABEL_SIZE, labelpad=2)
    # Count all known scaffolds, including those without ITS, before allowing labels.
    if len(set(chrom_sizes) | set(df["chr"])) > 10:
        return
    # Place sparse labels in display coordinates, keeping them apart on log axes.
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    placed = []
    for chrom in per.index:
        text = ax.annotate(chrom, (x[chrom], y[chrom]), xytext=(4, 0),
                           textcoords="offset points", ha="left", va="center",
                           fontsize=ANNOTATION_TEXT_SIZE, color="#222222")
        for dy in (0, 9, -9, 18, -18, 27, -27):
            for dx, align in ((4, "left"), (-4, "right")):
                text.set_position((dx, dy))
                text.set_ha(align)
                text.update_positions(renderer)
                box = matplotlib.text.Text.get_window_extent(text, renderer).expanded(1.04, 1.15)
                if (ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
                        and not any(box.overlaps(other) for other in placed)):
                    break
            else:
                continue
            break
        else:
            text.remove()
            continue
        placed.append(box)


def plot_its_summary_page(df, pairs, clusters, chrom_sizes, params):
    """ITS summary: headline tiles, length (a), scaffold density (b), and composition (c)."""
    fig = plt.figure(figsize=SLIDE_SIZE)
    _page_title(fig, "Assembly ITS summary")
    if df.empty:
        fig.text(0.5, 0.5, "No ITS", ha="center", va="center",
                 fontsize=PLACEHOLDER_TEXT_SIZE, color="#bbbbbb")
        return fig
    lengths = df["teloLen"].to_numpy(dtype=float)
    _draw_stat_tiles(_add_axes_in(fig, (FIG_WIDTH_DOUBLE - 6.0) / 2, 0.36, 6.0, 0.34), [
        (f"{_median(lengths):,.0f}", "Median ITS length (bp)"),
        (f"{len(df):,}", "ITS"),
        (f"{lengths.sum() / 1e3:,.1f}", "ITS content (kbp)"),
        (f"{df['chr'].nunique():,}", "Scaffolds with ITS"),
        (f"{len(clusters):,}", "Clusters"),
        (f"{len(pairs):,}", "Candidate fusions"),
    ])
    hist = _inch_axes(fig, 0.66, 0.63, 1.66, 1.66)  # level with panel c's joint square
    per = _inch_axes(fig, 2.98, 0.63, 1.66, 1.66)
    if (lengths > 0).any():
        _draw_its_length_panel(hist, lengths)
    else:
        _hide_panel(hist, "No ITS length")
    _panel_title(hist, "ITS length distribution", "a", x_in=PANEL_LETTER_X_IN, y=1 + 0.32 / 1.66)
    _draw_its_scaffold_panel(per, df, chrom_sizes)
    _panel_title(per, "ITS scaffold distribution", "b", y=1 + 0.32 / 1.66)  # titles level with c's marginal
    _draw_its_composition(fig, df, params)
    return fig


def _draw_its_joint_panel(ax, ax_top, ax_right, df, threshold):
    """Forward proportion x canonical proportion per ITS, on the double-key class regions, with marginals."""
    fwd, can = _its_shares(df)
    ok = np.isfinite(fwd)
    cuts = (-0.04, 1 - threshold, threshold, 1.04)
    for c in range(3):
        for s in range(3):
            ax.add_patch(Rectangle((cuts[s], cuts[c]), cuts[s + 1] - cuts[s], cuts[c + 1] - cuts[c],
                                   facecolor=DOUBLE_KEY_COLORS[c][s], alpha=0.25, linewidth=0, zorder=0))
    for cut in cuts[1:3]:
        ax.axvline(cut, color="#bcbcbc", lw=0.45, zorder=1)
        ax.axhline(cut, color="#bcbcbc", lw=0.45, zorder=1)
    lengths = df["teloLen"].to_numpy(dtype=float)
    dense = len(df) > 500
    order = np.argsort(-lengths[ok], kind="mergesort")
    colors = np.asarray(double_key_colors(df, threshold), dtype=object)[ok][order]
    ax.scatter(fwd[ok][order], can[ok][order],
               s=_its_marker_area(lengths[ok][order], lengths.max(), dense), c=list(colors),
               alpha=0.8 if dense else 1.0, edgecolors="#333333",
               linewidths=0.2, rasterized=dense, zorder=2)
    ticks = [0, 1 - threshold, threshold, 1]
    labels = [f"{v:.2f}".rstrip("0").rstrip(".") for v in ticks]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_ticks(ticks, labels)
    ax.set_xlim(cuts[-1], cuts[0])  # mirrored: fwd (5' end) on the left, as in the key
    ax.set_ylim(cuts[0], cuts[-1])
    _its_axis_style(ax)
    ax.set_xlabel("← Forward proportion", fontsize=AXIS_LABEL_SIZE, labelpad=2)
    ax.set_ylabel("Canonical proportion →", fontsize=AXIS_LABEL_SIZE, labelpad=2)

    bins = np.linspace(0, 1, 21)
    for marginal, values, orientation in ((ax_top, fwd[ok], "vertical"), (ax_right, can[ok], "horizontal")):
        counts = np.histogram(values, bins)[0]
        marginal.stairs(counts, bins, orientation=orientation, fill=True, facecolor="#9A9A9A", alpha=0.35, linewidth=0)
        marginal.stairs(counts, bins, orientation=orientation, color="#6F6F6F", linewidth=0.6)
    for marginal, line in ((ax_top, ax_top.axvline), (ax_right, ax_right.axhline)):
        for cut in cuts[1:3]:
            line(cut, color="#bcbcbc", lw=0.45)
        _its_axis_style(marginal, left=False)
        marginal.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    ax_right.spines["bottom"].set_visible(False)
    ax_right.spines["left"].set_visible(True)
    ax_top.set_xlim(ax.get_xlim())
    ax_right.set_ylim(ax.get_ylim())


def _draw_its_size_key(ax, lengths, dense):
    """Three reference markers for the length-to-area scale of the joint panel."""
    longest = lengths.max()
    top = 10 ** np.floor(np.log10(max(longest, 1)))
    refs = [v for v in (top / 100, top / 10, top) if v >= max(lengths.min(), 1)] or [longest]
    xs = np.arange(len(refs)) - (len(refs) - 1) / 2.0  # centred on the axes whatever the count
    ax.scatter(xs, np.zeros(len(refs)), s=_its_marker_area(np.array(refs), longest, dense),
               facecolors="none", edgecolors="#555555", linewidths=0.5)
    for x, ref in zip(xs, refs):
        ax.text(x, -0.9, _fmt_bp(ref), ha="center", va="top", fontsize=LEGEND_TEXT_SIZE, color="#222222")
    ax.text(0, 1.0, "Marker area scales with length", ha="center", va="bottom",
            fontsize=LEGEND_TEXT_SIZE, color="#222222")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-2.2, 1.6)
    ax.axis("off")


def _draw_its_composition(fig, df, params):
    """Summary panel c: joint composition and class totals in the same nine regions."""
    threshold = params.get("label_threshold", LABEL_THRESHOLD)
    lengths = df["teloLen"].to_numpy(dtype=float)
    joint = _inch_axes(fig, 5.21, 0.63, 1.66, 1.66)
    top = _inch_axes(fig, 5.21, 2.33, 1.66, 0.28)
    right = _inch_axes(fig, 6.91, 0.63, 0.16, 1.66)
    _draw_its_joint_panel(joint, top, right, df, threshold)
    _panel_title(top, "ITS composition", "c")
    strand, canon = composition_class(df["fwdCan"], df["revCan"], df["fwdNonCan"], df["revNonCan"], threshold)
    counts, bp = np.zeros((3, 3)), np.zeros((3, 3))
    np.add.at(counts, (canon, strand), 1)
    np.add.at(bp, (canon, strand), lengths)
    cuts = (-0.04, 1 - threshold, threshold, 1.04)
    for c in range(3):
        for s in range(3):
            share = 100 * bp[c, s] / lengths.sum() if lengths.sum() > 0 else 0
            joint.text((cuts[s] + cuts[s + 1]) / 2, (cuts[c] + cuts[c + 1]) / 2,
                       f"{int(counts[c, s]):,}\n{share:.0f}%", ha="center", va="center",
                       fontsize=MIN_TEXT_SIZE, color="black", zorder=4)


_ITS_SEGMENT_KEYS = (("fwdCan", 2, 2), ("fwdNonCan", 0, 2), ("revNonCan", 0, 0), ("revCan", 2, 0))
ITS_GLYPH_ROW_IN, ITS_GLYPH_BAR = 0.30, 0.19  # candidate glyph row pitch (in) and bar height (rows): label line over a thin bar


def _its_glyph_axes(ax, lo, hi, n_rows, label):
    """Glyph rows on a bp axis with plain numeric ticks and the unit in the axis label."""
    divisor, unit = _pick_bp_unit(max(-lo, hi, 1))
    ax.set_xlim(lo * 1.02, max(hi, 1) * 1.02)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{abs(v) / divisor:g}"))
    ax.set_xlabel(label.format(unit if unit == "bp" else unit + "p"), fontsize=AXIS_LABEL_SIZE, labelpad=1.5)
    _its_axis_style(ax, left=False)


def _its_row_labels(ax, i, left, right):
    """Identity label and one self-describing value on a line above a glyph row."""
    row_y = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    for x, text, ha, color in ((0, left, "left", "#222222"), (1, right, "right", "#555555")):
        ax.text(x, i - ITS_GLYPH_BAR / 2 - 0.08, text, transform=row_y, ha=ha, va="bottom", fontsize=AXIS_TICK_SIZE, color=color)


def _its_locus_name(df, chrom):
    """Display scaffold name, trimmed to fit the label column."""
    name = df.attrs.get("display_labels", {}).get(chrom, chrom)
    return name if len(name) <= 18 else name[:17] + "…"


def _draw_its_long_glyphs(ax, df):
    """Panel c: top ITS as to-scale bars split into their four strand/canonical segments, fwd first."""
    top = rank_long_its(df, 5)
    keys = ["chr", "start", "end", "teloLen", "teloLabel", "teloType"]
    counts = [k for k, _, _ in _ITS_SEGMENT_KEYS]
    top = top.merge(df.drop_duplicates(keys)[keys + counts], on=keys, how="left")
    _its_glyph_axes(ax, 0, top["teloLen"].max(), len(top), "ITS length ({})")
    colors = np.array([DOUBLE_KEY_COLORS[c][s] for _, c, s in _ITS_SEGMENT_KEYS])
    for i, row in enumerate(top.itertuples(index=False)):
        values = np.array([getattr(row, k) for k in counts], dtype=float)
        total = values.sum()
        if total > 0:
            widths = row.teloLen * values / total
            keep = widths > 0
            ax.barh(np.full(keep.sum(), i), widths[keep], left=np.r_[0, np.cumsum(widths)[:-1]][keep],
                    height=ITS_GLYPH_BAR, color=colors[keep], linewidth=0)
        else:
            ax.barh(i, row.teloLen, height=ITS_GLYPH_BAR, color=DOUBLE_KEY_COLORS[1][1], linewidth=0)
        # Count-based share, the composition page's y axis: canonical over all matches.
        _its_row_labels(ax, i, f"L{i + 1}  {_its_locus_name(df, row.chr)}",
                        f"{100 * (row.fwdCan + row.revCan) / total:.0f}% canonical" if total > 0 else "NA")


def _draw_its_cluster_glyphs(ax, df, clusters, threshold):
    """Panel a: each cluster as a to-scale track of its span with ITS cells in class colours."""
    top = clusters.head(5)
    _its_glyph_axes(ax, 0, float(top["span"].max()), len(top), "Position in cluster ({})")
    width_in = ax.get_position().width * ax.figure.get_size_inches()[0]
    min_width = ax.get_xlim()[1] * (1.2 / 72) / width_in
    for i, row in enumerate(top.itertuples(index=False)):
        ax.add_patch(Rectangle((0, i - 0.02), row.span, 0.04, facecolor=COLORS["gap"], linewidth=0))
        members = df[(df["chr"] == row.chr) & (df["start"] >= row.start) & (df["end"] <= row.end)]
        members = members.assign(dk_color=double_key_colors(members, threshold)).sort_values("teloLen", ascending=False)
        for m in members.itertuples(index=False):
            ax.add_patch(Rectangle((m.start - row.start, i - ITS_GLYPH_BAR / 2), max(m.end - m.start, min_width), ITS_GLYPH_BAR,
                                   facecolor=m.dk_color, linewidth=0))
        _its_row_labels(ax, i, f"C{i + 1}  {_its_locus_name(df, row.chr)}", f"{row.rows:,} ITS")


def _overlap_color(frame, start, end):
    """Class colour of the ITS overlapping [start, end) most, or the neutral ITS grey when none does."""
    overlap = np.minimum(frame["end"].to_numpy(), end) - np.maximum(frame["start"].to_numpy(), start)
    if not len(overlap) or overlap.max() <= 0:
        return COLORS["its"]
    return frame["dk_color"].iloc[int(np.argmax(overlap))]


def _draw_its_fusion_glyphs(ax, df, pairs, threshold):
    """Panel b: q block left and p block right of their junction, to scale."""
    top = pairs.head(5)
    half = top["spacer_bp"] / 2
    reach = float(np.maximum(top["q_bp"] + half, top["p_bp"] + half).max())
    _its_glyph_axes(ax, -reach, reach, len(top), "q  ←  Distance to junction ({})  →  p")
    ax.axvline(0, color="#dddddd", lw=0.45, zorder=0)
    colored = df.assign(dk_color=double_key_colors(df, threshold))
    for i, (row, h) in enumerate(zip(top.itertuples(index=False), half)):
        scaffold = colored[colored["chr"] == row.chr]
        ax.plot([-h, h], [i, i], color="#555555", lw=0.6, solid_capstyle="butt")
        for lo, x0, width in ((row.start, -h - row.q_bp, row.q_bp), (row.end - row.p_bp, h, row.p_bp)):
            ax.add_patch(Rectangle((x0, i - ITS_GLYPH_BAR / 2), width, ITS_GLYPH_BAR, facecolor=_overlap_color(scaffold, lo, lo + width),
                                   linewidth=0))
        _its_row_labels(ax, i, f"F{i + 1}  {_its_locus_name(df, row.chr)}", f"{row.spacer_bp:,} bp spacer")


def _its_segment_legend(ax):
    """One centred line naming the segment colours: hue is strand, lighter is less canonical."""
    handles = [Patch(facecolor=DOUBLE_KEY_COLORS[c][s], label=label) for c, s, label in
               ((2, 2, "fwd canonical"), (0, 2, "fwd non-canonical"), (2, 1, "both strands"),
                (0, 0, "rev non-canonical"), (2, 0, "rev canonical"))]
    ax.legend(handles=handles, loc="center", ncol=len(handles), frameon=False, fontsize=LEGEND_TEXT_SIZE,
              handlelength=0.9, handleheight=0.9, handletextpad=0.4, columnspacing=1.4, borderaxespad=0)
    ax.axis("off")


def plot_its_composition_page(df, pairs, clusters, params):
    """ITS candidates: clusters, candidate fusions and long ITS as to-scale glyph rows, in three columns."""
    fig = plt.figure(figsize=SLIDE_SIZE)
    _page_title(fig, "ITS candidates")
    threshold = params.get("label_threshold", LABEL_THRESHOLD)
    panels = (("Top ITS clusters", clusters, _draw_its_cluster_glyphs, (df, clusters, threshold)),
              ("Candidate fusions", pairs, _draw_its_fusion_glyphs, (df, pairs, threshold)),
              ("Longest canonical ITS", df, _draw_its_long_glyphs, (df,)))
    width, gap, x = 1.86, 0.44, 0.60
    for letter, (title, frame, draw, args) in zip("abc", panels):
        h = ITS_GLYPH_ROW_IN * max(1, min(5, len(frame)))
        ax = _inch_axes(fig, x, 2.60 - h, width, h)
        x += width + gap
        if frame.empty:
            _hide_panel(ax, "None")
        else:
            draw(ax, *args)
        _panel_title(ax, title, letter, x_in=PANEL_LETTER_X_IN if letter == "a" else None)
    if not df.empty:
        _its_segment_legend(_inch_axes(fig, 0, 3.02, SLIDE_SIZE[0], 0.14))
    return fig

# ---------------------------------------------------------------------------
# ITS-3 "Top loci"
# ---------------------------------------------------------------------------

def resolve_its_loci(clusters, df, pairs, chrom_sizes):
    """Resolve up to five padded locus windows per category, preserving candidate ranks."""
    loci = []
    for prefix, frame in (("C", clusters.head(5)), ("F", pairs.head(5)), ("L", rank_long_its(df, 5))):
        for rank, row in enumerate(frame.itertuples(index=False), 1):
            s, e = pad_window(int(row.start), int(row.end), chrom_sizes.get(row.chr, row.end))
            loci.append((f"{prefix}{rank}", row.chr, s, e))
    return loci


def _fmt_region(start, end):
    """Window as '117–171 kbp' in the locus axis unit, with just enough decimals to tell the ends apart."""
    divisor, unit = _region_unit(start, end)
    span = max(end - start, 1) / divisor
    dec = 0 if unit == "bp" else min(3, max(0, int(np.ceil(np.log10(10.0 / span)))))
    return f"{start / divisor:.{dec}f}–{end / divisor:.{dec}f} {unit}"


LOCUS_KEY_SIDE, LOCUS_KEY_TOP = 0.45, 0.14  # inches; a tiny key in the header's top-right corner


def _place_locus_key(fig, key_ax, side=LOCUS_KEY_SIDE, top=LOCUS_KEY_TOP):
    """Park the drawn 3x3 key in the header's top-right corner, its labels right-aligned with the tracks."""
    w, h = fig.get_size_inches()
    key_ax.set_position([0.97 - side / w, 1 - (top + side) / h, side / w, side / h])
    overhang = key_ax.get_tightbbox(fig.canvas.get_renderer()).x1 / fig.bbox.width - 0.97
    key_ax.set_position([0.97 - side / w - overhang, 1 - (top + side) / h, side / w, side / h])


def plot_its_loci_page(loci, chrom_sizes, arm_blocks, its_blocks, gap_blocks,
                       density_data, canonical_data, strand_data, gc_data, entropy_data,
                       unknown_extents=(), display_labels=None, threshold=LABEL_THRESHOLD):
    """Selected ITS locus: one terminal-style track-stack column per ranked locus."""
    if not loci:
        return _placeholder_figure("ITS loci", "No ITS loci to show.")

    track_specs = [("blocks", None)]
    if density_data:
        track_specs.append(("density", None))
    if canonical_data:
        track_specs.append(("canonical", None))
    if strand_data:
        track_specs.append(("strand", None))
    if gc_data:
        track_specs.append(("gc", None))
    if entropy_data:
        track_specs.append(("entropy", None))

    height_map = {"blocks": 0.20, "density": 0.34, "canonical": 0.34, "strand": 0.34,
                 "gc": 0.34, "entropy": 0.34}
    height_ratios = [0.34] + [height_map[name] for name, _ in track_specs]
    n_cols = len(loci)
    fig, axes = plt.subplots(
        len(height_ratios), n_cols, figsize=SLIDE_SIZE,
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.35,
                    "wspace": 0.30 if n_cols > 1 else 0.0},
        squeeze=False)
    fig.subplots_adjust(left=0.145, right=0.97, top=0.775, bottom=0.16)  # terminal page margins

    for col, (label, chrom, start, end) in enumerate(loci):
        size = int(chrom_sizes.get(chrom, end))
        blist = arm_blocks.get(chrom, [])
        itslist = its_blocks.get(chrom, []) if its_blocks else []
        gaplist = gap_blocks.get(chrom, []) if gap_blocks else []
        show_label = col == 0

        ideo_ax = axes[0][col]
        box_l, box_r, box_bottom = _draw_locus_ideogram(
            ideo_ax, size, start, end, blist, itslist,
            label=("Observed\nextent (Mbp)" if chrom in unknown_extents else "Scaffold\n(Mbp)")
            if show_label else None, threshold=threshold)
        ideo_ax.set_xlabel("")  # the blocks row would hide it; the unit sits in the track label
        _set_mb_ticks(ideo_ax, size)  # same ticks as the atlas
        # Trim the empty band under the bar so its tick labels clear the blocks row.
        cut = box_bottom - 0.04
        pos = ideo_ax.get_position()
        ideo_ax.set_position([pos.x0, pos.y0 + cut * pos.height, pos.width, (1 - cut) * pos.height])
        ideo_ax.set_ylim(cut, 1)
        ideo_ax.tick_params(axis="x", pad=2.0)

        axis_spec = _region_axis_spec(start, end, plain=True)
        visible_rows = []
        for i, (name, _) in enumerate(track_specs):
            row = i + 1
            ax = axes[row][col]
            if name == "blocks":
                visible = _draw_its_blocks_track(ax, start, end, blist, itslist, gaplist,
                                                 label="Blocks" if show_label else None,
                                                 threshold=threshold)
            elif name == "density":
                visible = _draw_fraction_track(ax, density_data.get(chrom), start, end, size,
                                               "full", COLORS["density"],
                                               label="Repeat\ndensity" if show_label else None)
            elif name == "canonical":
                visible = _draw_fraction_track(ax, canonical_data.get(chrom), start, end, size,
                                               "full", COLORS["canonical"],
                                               label="Canonical\nratio" if show_label else None)
            elif name == "gc":
                bg = gc_data.get(chrom)
                if bg is None or len(bg[0]) == 0:
                    ax.set_visible(False)
                    visible = False
                else:
                    starts, ends, values = bg
                    visible = _draw_fraction_track(ax, (starts, ends, values / 100.0), start, end,
                                                   size, "full", COLORS["gc"],
                                                   label="GC\ncontent" if show_label else None,
                                                   y_max=1.0, y_min=0.0)
            elif name == "entropy":
                visible = _draw_fraction_track(ax, entropy_data.get(chrom), start, end, size,
                                               "full", COLORS["entropy"],
                                               label="Shannon\nentropy" if show_label else None,
                                               y_max=2.0, y_min=1.0)
            else:
                visible = _draw_strand_track(ax, strand_data.get(chrom), start, end, size,
                                             "full", label="Strand\nbias" if show_label else None)
            if visible:
                ax.set_xlim(axis_spec["axis_start"], axis_spec["axis_end"])
                visible_rows.append(row)

        if visible_rows:
            for row in visible_rows[:-1]:
                _hide_x_axis(axes[row][col])
            _apply_terminal_x_axis(axes[visible_rows[-1]][col], axis_spec, "full")

        _panel_title(ideo_ax, f"{label}  {_fmt_region(start, end)}",
                     "abc"[col] if n_cols > 1 else None, x_in=PANEL_LETTER_X_IN if col == 0 else None)

    chroms = {chrom for _, chrom, _, _ in loci}
    if len(chroms) == 1:
        chrom = loci[0][1]
        size = int(chrom_sizes.get(chrom, loci[0][3]))
        _page_title(fig, f"{chrom}  ({'≥' if chrom in unknown_extents else ''}{_fmt_bp(size)})")
    else:
        _page_title(fig, "ITS loci")
    key_ax = fig.add_axes([0, 0, 1, 1])
    _draw_double_key(key_ax)
    _place_locus_key(fig, key_ax)
    return fig

# ---------------------------------------------------------------------------
# Tier 2: Terminal zoom figures
# ---------------------------------------------------------------------------

def compute_view_windows(blocks_list, chrom_size):
    """Compute (p_window, q_window) for terminal zoom panels.

    Each window is (start, end) or None if no blocks at that end.
    Target: ~50% telomeric content per panel (PADDING_FACTOR = 2.0).
    When both arms exist and extents are within MAX_RATIO, a common limit
    is used so arm lengths can be compared directly; otherwise independent
    per-arm limits preserve detail for the smaller arm.
    Only considers arm telomere blocks (term=="scaffold"), not contig rows.
    """
    PADDING_FACTOR = 2.0
    MIN_WINDOW = 1_000
    MAX_RATIO = 5.0
    ROUND_TO_BP = 1_000

    def _normalize_terminal_limit_bp(limit_bp):
        if limit_bp is None:
            return None
        chrom_limit = int(chrom_size)
        if chrom_limit <= 0:
            return None
        limit_bp = min(chrom_limit, max(int(limit_bp), min(MIN_WINDOW, chrom_limit)))
        if chrom_limit <= ROUND_TO_BP:
            return chrom_limit
        rounded = int(np.ceil(limit_bp / float(ROUND_TO_BP)) * ROUND_TO_BP)
        return min(chrom_limit, rounded)

    def _arm_limit(arm_blocks, arm):
        if not arm_blocks:
            return None
        intervals = [
            _terminal_distance_interval(
                int(b["start"]), int(b["end"]), int(chrom_size), arm)
            for b in arm_blocks
        ]
        furthest = max(max(s, e) for s, e in intervals)
        padded = min(int(chrom_size), max(int(round(furthest * PADDING_FACTOR)), MIN_WINDOW))
        return _normalize_terminal_limit_bp(padded)

    p_blocks = [b for b in blocks_list if b.get("term") != "contig" and b["closestEnd"] == "p"]
    q_blocks = [b for b in blocks_list if b.get("term") != "contig" and b["closestEnd"] == "q"]
    p_limit = _arm_limit(p_blocks, "p")
    q_limit = _arm_limit(q_blocks, "q")

    if p_limit is not None and q_limit is not None:
        bigger = max(p_limit, q_limit)
        smaller = max(min(p_limit, q_limit), 1)
        if bigger / smaller <= MAX_RATIO:
            # Similar scale: common limit for direct visual comparison
            p_window = (0, bigger)
            q_window = (max(0, int(chrom_size) - bigger), chrom_size)
        else:
            # Very different scales: independent limits to preserve detail
            p_window = (0, p_limit)
            q_window = (max(0, int(chrom_size) - q_limit), chrom_size)
    elif p_limit is not None:
        p_window = (0, p_limit)
        q_window = (max(0, int(chrom_size) - p_limit), chrom_size)
    elif q_limit is not None:
        p_window = (0, q_limit)
        q_window = (max(0, int(chrom_size) - q_limit), chrom_size)
    else:
        p_window = None
        q_window = None

    return p_window, q_window


def clip_bedgraph(bg_tuple, view_start, view_end):
    """Return (starts, ends, values) arrays clipped to [view_start, view_end)."""
    starts, ends, values = bg_tuple
    mask = (ends > view_start) & (starts < view_end)
    cs = np.maximum(starts[mask], view_start)
    ce = np.minimum(ends[mask], view_end)
    return cs, ce, values[mask]


def _draw_blocks_track(ax, blocks_list, view_start, view_end, chrom_size, arm, label=None,
                       telomere_present=True, its_blocks_list=None, gap_blocks_list=None,
                       contig_blocks_list=None):
    """Draw telomere blocks as a thin track on a backbone line.

    Nature-style monochrome gradient:
      terminal blocks → COLORS["terminal"] (dark, filled)
      contig rows     → COLORS["terminal"] (dark, outline-only)
      ITS blocks      → COLORS["its"]      (medium grey)
      gap blocks      → COLORS["gap"]      (light grey)
    p/q arm letters are centered on each block when the block is wide enough.
    """
    backbone_y = 0.5
    display_start, display_end = _project_terminal_interval(
        view_start, view_end, chrom_size, arm)
    view_span = abs(display_end - display_start) if display_end != display_start else 1
    ax.plot([display_start, display_end], [backbone_y, backbone_y],
            color="#d6d6d6", linewidth=0.9, solid_capstyle="round",
            zorder=1)

    # ---- Terminal blocks (p / q / balanced) ----
    for b in blocks_list:
        if b["end"] <= view_start or b["start"] >= view_end:
            continue
        cs = max(b["start"], view_start)
        ce = min(b["end"], view_end)
        ds, de = _project_terminal_interval(cs, ce, chrom_size, arm)
        if de < ds:
            ds, de = de, ds
        rect = Rectangle(
            (ds, backbone_y - 0.07), de - ds, 0.14,
            facecolor=COLORS["terminal"], edgecolor="none", alpha=0.98, zorder=2,
        )
        ax.add_patch(rect)
        if _block_symbol_fits(de - ds, view_span, b["label"]):
            _draw_block_symbol(ax, (ds + de) / 2, backbone_y, b["label"], zorder=4)

    # ---- Contig rows (outline-only) ----
    if contig_blocks_list:
        for b in contig_blocks_list:
            if b["end"] <= view_start or b["start"] >= view_end:
                continue
            cs = max(b["start"], view_start)
            ce = min(b["end"], view_end)
            ds, de = _project_terminal_interval(cs, ce, chrom_size, arm)
            if de < ds:
                ds, de = de, ds
            rect = Rectangle(
                (ds, backbone_y - 0.07), de - ds, 0.14,
                facecolor="none", edgecolor=COLORS["terminal"], linewidth=0.6, zorder=3,
            )
            ax.add_patch(rect)

    # ---- ITS blocks ----
    if its_blocks_list:
        for b in its_blocks_list:
            if b["end"] <= view_start or b["start"] >= view_end:
                continue
            cs = max(b["start"], view_start)
            ce = min(b["end"], view_end)
            ds, de = _project_terminal_interval(cs, ce, chrom_size, arm)
            if de < ds:
                ds, de = de, ds
            rect = Rectangle(
                (ds, backbone_y - 0.07), de - ds, 0.14,
                facecolor=COLORS["its"], edgecolor="none",
                alpha=0.98, zorder=3,
            )
            ax.add_patch(rect)
            if _block_symbol_fits(de - ds, view_span, b.get("label", "")):
                _draw_block_symbol(ax, (ds + de) / 2, backbone_y, b.get("label", ""), zorder=5)

    # ---- Gap blocks ----
    if gap_blocks_list:
        for b in gap_blocks_list:
            if b["end"] <= view_start or b["start"] >= view_end:
                continue
            cs = max(b["start"], view_start)
            ce = min(b["end"], view_end)
            ds, de = _project_terminal_interval(cs, ce, chrom_size, arm)
            if de < ds:
                ds, de = de, ds
            rect = Rectangle(
                (ds, backbone_y - 0.07), de - ds, 0.14,
                facecolor=COLORS["gap"], edgecolor="none",
                alpha=0.98, zorder=4,
            )
            ax.add_patch(rect)

    ax.set_ylim(0.28, 0.72)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    _set_track_label(ax, label)
    if not telomere_present:
        ax.text(0.5, 0.84, "No telomere", transform=ax.transAxes,
                ha="center", va="center", fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999")
    _hide_x_axis(ax)
    return True


def _draw_fraction_track(ax, track_data, view_start, view_end, chrom_size, arm, color, label=None, y_max=1.0, y_min=0.0):
    """Draw a quantitative track as a thin step plot with light fill."""
    starts, ends, values = clip_bedgraph(track_data, view_start, view_end)
    if len(starts) == 0:
        ax.set_visible(False)
        return False

    starts, ends, values = _project_terminal_series(
        starts, ends, values, chrom_size, arm)

    valid = np.isfinite(values) & (values >= 0)
    no_data = ~valid

    if np.any(no_data):
        # One collection for all no-data spans, full axes height like axvspan.
        spans = [((a, 0), (b, 0), (b, 1), (a, 1)) for a, b in zip(starts[no_data], ends[no_data])]
        ax.add_collection(PolyCollection(spans, facecolor=COLORS["no_data"], alpha=0.22, linewidth=0, zorder=0,
                                         transform=ax.get_xaxis_transform()), autolim=False)

    # One fill (runs split by NaN) and one line collection per track, however fragmented the data.
    pieces = [_bedgraph_to_step(starts[a:b], ends[a:b], values[a:b]) for a, b in _iter_true_runs(valid)]
    if pieces:
        ax.fill_between(np.concatenate([np.append(xs, xs[-1]) for xs, _ in pieces]),
                        np.concatenate([np.append(ys, np.nan) for _, ys in pieces]),
                        step="post", color=color, alpha=0.16, linewidth=0)
        ax.add_collection(LineCollection([np.column_stack(cbook.pts_to_poststep(xs, ys)) for xs, ys in pieces],
                                         colors=color, linewidths=0.6, capstyle=plt.rcParams["lines.solid_capstyle"],
                                         joinstyle=plt.rcParams["lines.solid_joinstyle"]), autolim=False)

    _style_fraction_axis(ax, label, y_max=y_max, y_min=y_min)
    return True


def _draw_strand_track(ax, strand_data, view_start, view_end, chrom_size, arm, label=None):
    """Draw strand bias = |fwdRatio - revRatio| as a 0..1 track."""
    starts, ends, values = clip_bedgraph(strand_data, view_start, view_end)
    if len(starts) == 0:
        ax.set_visible(False)
        return False

    bias_values = values.copy()
    valid = np.isfinite(values) & (values >= 0)
    bias_values[valid] = np.abs((2.0 * values[valid]) - 1.0)
    return _draw_fraction_track(
        ax,
        (starts, ends, bias_values),
        view_start,
        view_end,
        chrom_size,
        arm,
        COLORS["strand_bias"],
        label,
    )


def _hide_panel(ax, message="No telomere"):
    """Hide axes and show a centered gray label."""
    ax.set_visible(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.5, 0.5, message, transform=ax.transAxes,
            ha="center", va="center", fontsize=PLACEHOLDER_TEXT_SIZE, color="#bbbbbb")


def plot_terminal_zoom(chrom, chrom_size, blocks_list,
                       density_data=None, canonical_data=None, strand_data=None,
                       its_blocks_list=None, gc_data=None, entropy_data=None,
                       gap_blocks_list=None, contig_blocks_list=None):
    """Two-column terminal zoom: p-end (left) and q-end (right).

    Each column shows tracks (blocks, density, canonical ratio, strand bias, GC, entropy).
    If windows overlap on a short chromosome, a single merged panel is used.
    contig_blocks_list: optional contig-terminal rows to draw as outline-only.
    """
    has_density = density_data is not None and len(density_data[0]) > 0
    has_canonical = canonical_data is not None and len(canonical_data[0]) > 0
    has_strand = strand_data is not None and len(strand_data[0]) > 0
    has_its = its_blocks_list is not None and len(its_blocks_list) > 0
    has_gc = gc_data is not None and len(gc_data[0]) > 0
    has_entropy = entropy_data is not None and len(entropy_data[0]) > 0
    p_has_telomere = any(b["closestEnd"] == "p" for b in blocks_list)
    q_has_telomere = any(b["closestEnd"] == "q" for b in blocks_list)

    p_window, q_window = compute_view_windows(blocks_list, chrom_size)
    fallback_full_chrom = p_window is None and q_window is None
    if fallback_full_chrom:
        p_window = (0, chrom_size)
    overlap = p_window is not None and q_window is not None and p_window[1] >= q_window[0]
    merged = fallback_full_chrom or (overlap and p_has_telomere and q_has_telomere)

    track_specs = [("blocks", None)]
    if has_density:
        track_specs.append(("density", density_data))
    if has_canonical:
        track_specs.append(("canonical", canonical_data))
    if has_strand:
        track_specs.append(("strand", strand_data))
    if has_gc:
        track_specs.append(("gc", gc_data))
    if has_entropy:
        track_specs.append(("entropy", entropy_data))

    n_tracks = len(track_specs)

    n_cols = 1 if merged else 2
    height_map = {
        "blocks":    0.18,
        "density":   0.30,
        "canonical": 0.30,
        "strand":    0.30,
        "gc":        0.30,
        "entropy":   0.30,
    }
    height_ratios = [height_map[name] for name, _ in track_specs]

    fig_height = 2.2 + 0.25 * n_tracks
    fig, axes = plt.subplots(
        n_tracks, n_cols,
        figsize=(FIG_WIDTH_DOUBLE, fig_height),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.32,
                     "wspace": 0.18},
        squeeze=False,
        sharey="row",
    )

    fig.subplots_adjust(left=0.145, right=0.97, top=0.775, bottom=0.16)

    # Determine which columns to draw
    if merged:
        # Single merged panel
        columns = [(0, (0, chrom_size), "full")]
    else:
        columns = []
        if p_window:
            columns.append((0, p_window, "p"))
        if q_window:
            columns.append((1 if n_cols == 2 else 0, q_window, "q"))

    label_col = columns[0][0] if columns else 0

    # Track which column indices are actually drawn
    drawn_cols = set()
    for col_idx, window, arm in columns:
        drawn_cols.add(col_idx)
        view_start, view_end = window
        display_start, display_end = _project_terminal_interval(
            view_start, view_end, chrom_size, arm)
        axis_spec = _get_terminal_axis_spec(display_start, display_end, arm)
        visible_rows = []

        for row, (track_name, track_data) in enumerate(track_specs):
            ax = axes[row][col_idx]
            show_label = col_idx == label_col

            if track_name == "blocks":
                telomere_present = (
                    True if arm == "full" else
                    p_has_telomere if arm == "p" else
                    q_has_telomere
                )
                visible = _draw_blocks_track(
                    ax, blocks_list, view_start, view_end, chrom_size, arm,
                    label="Blocks" if show_label else None,
                    telomere_present=telomere_present,
                    its_blocks_list=its_blocks_list if has_its else None,
                    gap_blocks_list=gap_blocks_list,
                    contig_blocks_list=contig_blocks_list,
                )
            elif track_name == "density":
                visible = _draw_fraction_track(
                    ax, track_data, view_start, view_end,
                    chrom_size, arm,
                    COLORS["density"],
                    label="Repeat\ndensity" if show_label else None,
                )
            elif track_name == "canonical":
                visible = _draw_fraction_track(
                    ax, track_data, view_start, view_end,
                    chrom_size, arm,
                    COLORS["canonical"],
                    label="Canonical\nratio" if show_label else None,
                )
            elif track_name == "gc":
                starts, ends, values = track_data
                gc_frac = values / 100.0
                visible = _draw_fraction_track(
                    ax, (starts, ends, gc_frac), view_start, view_end,
                    chrom_size, arm,
                    COLORS["gc"],
                    label="GC\ncontent" if show_label else None,
                    y_max=1.0,
                    y_min=0.0,
                )
            elif track_name == "entropy":
                visible = _draw_fraction_track(
                    ax, track_data, view_start, view_end,
                    chrom_size, arm,
                    COLORS["entropy"],
                    label="Shannon\nentropy" if show_label else None,
                    y_max=2.0,
                    y_min=1.0,
                )
            else:
                visible = _draw_strand_track(
                    ax, track_data, view_start, view_end,
                    chrom_size, arm,
                    label="Strand\nbias" if show_label else None,
                )

            if visible:
                if axis_spec is not None:
                    if arm == "q":
                        ax.set_xlim(axis_spec["axis_end"], axis_spec["axis_start"])
                    else:
                        ax.set_xlim(axis_spec["axis_start"], axis_spec["axis_end"])
                visible_rows.append(row)

        if visible_rows:
            for row in visible_rows[:-1]:
                _hide_x_axis(axes[row][col_idx])
            _apply_terminal_x_axis(
                axes[visible_rows[-1]][col_idx],
                axis_spec,
                arm,
            )

    # Hide columns with no drawable data
    if n_cols == 2 and not merged:
        for col_idx in range(n_cols):
            if col_idx not in drawn_cols:
                for row in range(n_tracks):
                    _hide_panel(axes[row][col_idx], message="")

    # Column titles
    if merged:
        axes[0][0].set_title("Full scaffold", fontsize=PANEL_TITLE_SIZE, pad=3)
    else:
        if n_cols == 2:
            axes[0][0].set_title("p-arm", fontsize=PANEL_TITLE_SIZE, pad=3)
            axes[0][1].set_title("q-arm", fontsize=PANEL_TITLE_SIZE, pad=3)

    if not merged and n_cols == 2:
        fig.text(0.033, 0.825, "a", fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="bottom", ha="left")
        fig.text(0.533, 0.825, "b", fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="bottom", ha="left")

    fig.suptitle(f"{chrom}  ({_fmt_bp(chrom_size)})",
                 fontsize=PANEL_LABEL_SIZE, fontweight="bold", y=0.958, x=0.5, ha="center")
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate publication-ready figures from Teloscope output.",
        epilog="Example: python teloscope_report.py output/ -o report.pdf",
    )
    parser.add_argument("directory", help="Teloscope output directory")
    parser.add_argument("-o", "--output", default=None,
                        help="PDF stem in split mode, exact path for one PDF, or PNG directory")
    parser.add_argument("--png", action="store_true",
                        help="Save individual PNG files instead of a single PDF")
    parser.add_argument("--section", choices=("split", "all", "terminal", "its"), default="split",
                        help="Default: separate PDFs; all: combined PDF; terminal/its: one section")
    parser.add_argument("--dpi", type=int, default=450,
                        help="DPI for raster output (default: 450)")
    parser.add_argument("--draft", action="store_true",
                        help="Draft mode: render at 150 DPI for fast iteration")
    args = parser.parse_args()

    if args.draft:
        args.dpi = 150
    if args.dpi <= 0:
        parser.error("--dpi must be positive")

    if not os.path.isdir(args.directory):
        sys.exit(f"Error: '{args.directory}' is not a directory.")

    files = find_files(args.directory)
    if (args.section == "terminal" and "terminal" not in files or
            args.section == "its" and "interstitial" not in files or
            not ({"terminal", "interstitial"} & set(files))):
        sys.exit(f"Error: Missing teloscope output files in '{args.directory}'.\n"
                 f"Run teloscope first to generate output files.")
    if "report" not in files and args.section != "its":
        _warn(f"No '*_report.tsv' file found in '{args.directory}'; the overview classification panel will show no data.")

    blocks = parse_terminal_bed(files["terminal"]) if "terminal" in files else {}
    include_terminal = args.section != "its" and "terminal" in files

    # Split terminal blocks into arm telomeres (scaffold) and contig-terminal rows (contig)
    arm_blocks = {}
    contig_blocks = {}
    for chrom, blist in blocks.items():
        arm_blist = [b for b in blist if b.get("term") == "scaffold"]
        contig_blist = [b for b in blist if b.get("term") == "contig"]
        if arm_blist:
            arm_blocks[chrom] = arm_blist
        if contig_blist:
            contig_blocks[chrom] = contig_blist

    # The C parser releases the GIL while tokenising, so the five BEDgraph tracks load concurrently.
    tracks = ("density", "canonical_ratio", "strand_ratio", "gc", "entropy")
    with ThreadPoolExecutor(max_workers=len(tracks)) as pool:
        loading = {k: pool.submit(parse_bedgraph, files[k]) for k in tracks if k in files}
        its_blocks = parse_terminal_bed(files["interstitial"]) if "interstitial" in files else None
        gap_blocks = parse_interval_bed(files["gaps"]) if "gaps" in files else None
    density_data, canonical_data, strand_data, gc_data, entropy_data = (
        loading[k].result() if k in loading else None for k in tracks)

    chrom_sizes = get_chrom_sizes(blocks, density_data, canonical_data, strand_data)

    classifications = parse_report(files["report"]) if "report" in files else OrderedDict()
    # Only generate figures for chromosomes that have telomere blocks
    profile_chroms = [chrom for chrom in sorted(arm_blocks.keys(),
                                                key=lambda c: chrom_sizes.get(c, 0), reverse=True)
                      if chrom_sizes.get(chrom, 0) > 0]
    skipped_no_size = sorted(set(arm_blocks) - set(profile_chroms))
    if skipped_no_size:
        _warn(f"Skipping {len(skipped_no_size)} chromosome(s) with no usable size: {', '.join(skipped_no_size[:5])}"
              f"{' ...' if len(skipped_no_size) > 5 else ''}.")

    # ITS pages: skipped entirely when the interstitial BED is missing or empty.
    its_page = None
    if "interstitial" in files and args.section != "terminal":
        params = read_params(files.get("report"))
        its_frame = load_its_frame(files["interstitial"], motif_len=params["motif_len"], blocks=its_blocks)
        if len(its_frame) or args.section in ("its", "split") or not include_terminal:
            extents = its_frame.groupby("chr").agg(size=("chrSize", "max"), end=("end", "max"))
            for chrom, row in extents.iterrows():
                chrom_sizes[chrom] = max(chrom_sizes.get(chrom, 0), int(row["size"]), int(row["end"]))
            gaps_frame = load_gaps_frame(files["gaps"]) if "gaps" in files else pd.DataFrame(columns=["chr", "start", "end"])
            pairs = pair_fusions(its_frame, gaps_frame, params["max_block_dist"])
            clusters = summarize_its_clusters(its_frame)
            its_page = (its_frame, pairs, clusters, params)
            its_frame.attrs["display_labels"] = _its_labels(its_frame, chrom_sizes, arm_blocks)

    fallback_pages = []

    if args.png:
        out_dir = args.output or args.directory
    else:
        default_name = ("teloscope_report.pdf" if args.section in ("all", "split")
                        else f"teloscope_{args.section}_report.pdf")
        out_path = args.output or os.path.join(args.directory, default_name)
        out_dir = os.path.dirname(out_path) or "."
    os.makedirs(out_dir, exist_ok=True)

    if its_page is not None:
        interstitial_base = os.path.basename(files["interstitial"])
        suffix = "_interstitial_telomeres.bed"
        prefix = interstitial_base[:-len(suffix)] if interstitial_base.endswith(suffix) else os.path.splitext(interstitial_base)[0]
        its_frame, pairs, clusters, params = its_page
        os.makedirs(out_dir, exist_ok=True)
        its_top_hits(os.path.join(out_dir, f"{prefix}_its_top_hits.tsv"), its_frame, pairs, clusters)
        its_frame.to_csv(os.path.join(out_dir, f"{prefix}_its_rows.tsv"), sep="\t", index=False, na_rep="NA")
        its_scaffold_summary(its_frame, arm_blocks, chrom_sizes).to_csv(
            os.path.join(out_dir, f"{prefix}_its_scaffolds.tsv"), sep="\t", index=False)

    # --- One page list shared by the PDF and PNG branches ---
    pages = [
        ("overview-1", "teloscope_overview_1.png",
         lambda: plot_overview_page1(classifications, arm_blocks, chrom_sizes),
         "Assembly overview (page 1)",
         "Failed to render overview page 1. A placeholder page was written instead."),
        ("overview-2", "teloscope_overview_2.png",
         lambda: plot_overview_page2(arm_blocks, chrom_sizes),
         "Assembly overview (page 2)",
         "Failed to render overview page 2. A placeholder page was written instead."),
    ]
    if not include_terminal:
        pages = []
    for chrom in (profile_chroms if include_terminal else []):
        csize = chrom_sizes.get(chrom, 0)
        if csize == 0:
            continue
        arm_blist = arm_blocks.get(chrom, [])
        contig_blist = contig_blocks.get(chrom, [])
        den = density_data.get(chrom) if density_data else None
        can = canonical_data.get(chrom) if canonical_data else None
        strand = strand_data.get(chrom) if strand_data else None
        its = its_blocks.get(chrom, []) if its_blocks else None
        gaps = gap_blocks.get(chrom, []) if gap_blocks else None
        gc = gc_data.get(chrom) if gc_data else None
        ent = entropy_data.get(chrom) if entropy_data else None
        pages.append((
            chrom, f"teloscope_{_sanitize_filename(chrom)}.png",
            lambda chrom=chrom, csize=csize, arm_blist=arm_blist, contig_blist=contig_blist,
                   den=den, can=can, strand=strand, its=its, gc=gc, ent=ent, gaps=gaps:
                plot_terminal_zoom(chrom, csize, arm_blist, den, can, strand, its, gc, ent,
                                   gaps, contig_blocks_list=contig_blist),
            chrom,
            f"Failed to render the terminal zoom for {chrom}. A placeholder page was written instead.",
        ))

    terminal_page_count = len(pages)
    if its_page is not None:
        its_frame, pairs, clusters, params = its_page
        loci = resolve_its_loci(clusters, its_frame, pairs, chrom_sizes)
        unknown_extents = set(its_frame.loc[its_frame["chrSize"] == 0, "chr"])
        pages.append((
            "its-summary", "teloscope_its_summary.png",
            lambda: plot_its_summary_page(its_frame, pairs, clusters, chrom_sizes, params),
            "ITS summary", "Failed to render ITS summary.",
        ))
        # Atlas geometry and rankings are dataset-wide; compute once here rather than per page.
        atlas = _its_atlas_chroms(its_frame, arm_blocks, chrom_sizes)
        cells = _its_atlas_cells(its_frame, params["label_threshold"])
        top_hits = _its_top_hits(its_frame, pairs, clusters)
        known_sizes = its_frame.groupby("chr")["chrSize"].max().to_dict()
        atlas_pages = _paginate_atlas(atlas, _atlas_extent(atlas, chrom_sizes, cells))
        atlas_count = len(atlas_pages)
        for page_idx, chroms in enumerate(atlas_pages):
            pages.append((
                f"its-atlas-{page_idx + 1}", f"teloscope_its_atlas_{page_idx + 1:03d}.png",
                lambda chroms=chroms, page_idx=page_idx: plot_its_overview_page(
                    its_frame, pairs, clusters, arm_blocks, chrom_sizes, params,
                    chroms, page_idx + 1, atlas_count, cells, atlas, top_hits, known_sizes),
                f"ITS atlas ({page_idx + 1}/{atlas_count})", "Failed to render ITS atlas.",
            ))
        pages.append((
            "its-candidates", "teloscope_its_candidates.png",
            lambda: plot_its_composition_page(its_frame, pairs, clusters, params),
            "ITS candidates", "Failed to render ITS candidates.",
        ))
        for locus in loci:
            pages.append((
                f"its-locus-{locus[0]}", f"teloscope_its_locus_{locus[0]}.png",
                lambda locus=locus: plot_its_loci_page(
                    [locus], chrom_sizes, arm_blocks, its_blocks, gap_blocks,
                    density_data, canonical_data, strand_data, gc_data, entropy_data,
                    unknown_extents, its_frame.attrs["display_labels"], params["label_threshold"]),
                f"ITS locus {locus[0]}", "Failed to render ITS locus.",
            ))

    n_figures = len(pages)
    print("\n+++ Plot Report +++", file=sys.stderr)

    def emit_page(save_fig, builder, title, message, name, index):
        """Build, save via save_fig, and report progress/fallback for one page."""
        ok, error_text = _save_figure_with_fallback(save_fig, builder, title, message)
        if not ok:
            fallback_pages.append((name, error_text))
        print(f"[{index}/{n_figures}] {title}{' [warning]' if not ok else ''}", file=sys.stderr)

    # --- Write-and-close pattern: one figure in memory at a time ---
    if args.png:
        os.makedirs(out_dir, exist_ok=True)
        used_names = set()
        for i, (name, png_name, builder, title, message) in enumerate(pages, start=1):
            if png_name in used_names:
                png_name = f"{i:04d}_{png_name}"
            used_names.add(png_name)
            path = os.path.join(out_dir, png_name)
            emit_page(lambda fig, path=path: fig.savefig(path, dpi=args.dpi), builder, title, message, name, i)
        print(f"Saved: {out_dir}/", file=sys.stderr)
    else:
        if args.section == "split":
            stem = os.path.splitext(out_path)[0]
            groups = [(stem + "_terminal.pdf", pages[:terminal_page_count]),
                      (stem + "_its.pdf", pages[terminal_page_count:])]
        else:
            groups = [(out_path, pages)]
        i = 0
        for pdf_path, section_pages in groups:
            if not section_pages:
                continue
            with PdfPages(pdf_path) as pdf:
                for name, png_name, builder, title, message in section_pages:
                    i += 1
                    emit_page(lambda fig: pdf.savefig(fig, dpi=args.dpi), builder, title, message, name, i)
            print(f"Saved: {pdf_path}", file=sys.stderr)

    if fallback_pages:
        preview = ", ".join(f"{name} ({err})" for name, err in fallback_pages[:5])
        more = " ..." if len(fallback_pages) > 5 else ""
        _warn(f"Report generation completed with {len(fallback_pages)} placeholder figure(s): {preview}{more}")


if __name__ == "__main__":
    main()
