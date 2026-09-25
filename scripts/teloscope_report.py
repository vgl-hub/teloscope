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

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch, ConnectionPatch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection, PatchCollection
import matplotlib.ticker as ticker
import matplotlib.patheffects as patheffects
from matplotlib import transforms
from matplotlib.colors import LinearSegmentedColormap, LogNorm

# ---------------------------------------------------------------------------
# Nature-style configuration
# ---------------------------------------------------------------------------

# Nature figure widths: 89 mm (single column), 183 mm (double column)
FIG_WIDTH_SINGLE = 3.50    # inches (89 mm)
FIG_WIDTH_DOUBLE = 7.20    # inches (183 mm)
REPORT_PAGE_HEIGHT = 2.2 + 0.25 * 6  # terminal page with all six tracks
ITS_ATLAS_ROWS = 20  # two columns of ten; paginate instead of shrinking text
ITS_POSITION_BINS = 100

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

# ITS junction classes, fixed order; fusion is the one accent hue, the rest recede to greys picked for >=3:1 contrast on white.
CLASS_ORDER = ["fusion", "tail_to_tail", "fragmentation", "single"]
CLASS_SHORT = {"fusion": "fusion", "tail_to_tail": "t2t", "fragmentation": "frag", "single": "single"}
CLASS_COLORS = {
    "fusion":        COLORS["q"],
    "tail_to_tail":  COLORS["terminal"],  # #4A4A4A, contrast 8.9:1
    "fragmentation": "#6E6E6E",           # contrast 5.1:1
    "single":        COLORS["its"],       # #8C8C8C, contrast 3.4:1
}

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

FIGURE_TITLE_SIZE = 9.3
FIGURE_SUMMARY_SIZE = 6.0
PANEL_LABEL_SIZE = 8.0
PANEL_TITLE_SIZE = 6.9
AXIS_LABEL_SIZE = 6.2
AXIS_TICK_SIZE = 5.5
PANEL_D_TICK_SIZE = 5.0  # panel d's two narrow log-x sub-panels need smaller decade ticks
LEGEND_TEXT_SIZE = 5.4
TABLE_TEXT_SIZE = 6.5
TABLE_ROW_IN = 0.135  # fixed physical row height (header + data rows alike) for ITS tables
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
                         dtype={"chrom": str, "start": np.int64,
                                "end": np.int64, "value": np.float64},
                         engine="c", on_bad_lines="skip")
    except pd.errors.EmptyDataError:
        _warn(f"BEDgraph file '{path}' is empty; skipping.")
        return OrderedDict()
    except Exception as exc:
        _warn(f"Fast BEDgraph parse failed for '{path}' ({exc}); retrying line-by-line.")
        return _parse_bedgraph_fallback(path)

    valid = ((df["start"] >= 0) & (df["end"] > df["start"]) &
             (df["end"] <= MAX_COORD) & np.isfinite(df["value"]))
    if not valid.all():
        _warn(f"Skipped {int((~valid).sum())} invalid BEDgraph row(s) from '{path}'.")
    df = df.loc[valid].sort_values(["chrom", "start"], kind="mergesort")
    data = OrderedDict()
    for chrom, grp in df.groupby("chrom", sort=False):
        data[chrom] = (grp["start"].values, grp["end"].values, grp["value"].values)
    return data


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
              "motif_len": 6, "min_canonical_count": 4, "known_params": set()}
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

        fh.write(f"# section 2: {top_longest} longest interstitial telomere rows by canonical bp "
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
    """Pick one bp/kb/Mb unit for a value or panel from its largest value."""
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
    decimals = 2 if unit == "Mb" else 1
    text = f"{bp / divisor:.{decimals}f}".rstrip("0").rstrip(".")
    return f"{text} {unit}"


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


def _panel_label(ax, label, x=-0.05, y=1.15):
    """Add Nature-style panel label (bold lowercase letter, top-left)."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="top", ha="left")


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
    left_legend.set_bbox_to_anchor((left_anchor, 0.50), transform=ax.transAxes)
    right_legend.set_bbox_to_anchor((right_anchor, 0.50), transform=ax.transAxes)
    return left_legend, right_legend


def _draw_ranked_bar_panel(ax, rows, x_label, xticks, title=None):
    """Draw a ranked horizontal bar panel with a shared Nature-style axis treatment."""
    label_space = 1.42
    y_pos = np.arange(len(rows), dtype=np.float64)
    values = [row["value"] for row in rows]
    colors = [row["color"] for row in rows]

    ax.barh(
        y_pos,
        values,
        height=0.20,
        color=colors,
        edgecolor="white",
        linewidth=0.45,
        zorder=2,
    )
    for y_pos_i, row in zip(y_pos, rows):
        ax.text(
            -0.08,
            y_pos_i,
            row["label"],
            ha="right",
            multialignment="right",
            va="center",
            fontsize=MIN_TEXT_SIZE,
            color="#222222",
            linespacing=1.0,
        )

    x_max = max(max(xticks), max(values) + 0.12) if values else max(xticks)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_yticks([])
    ax.set_xlim(-label_space, x_max)
    ax.set_xlabel(x_label, fontsize=AXIS_LABEL_SIZE)
    ax.set_xticks([tick for tick in xticks if tick <= x_max + 0.1])
    ax.tick_params(axis="x", length=2.0, width=0.35, pad=0.6, labelsize=AXIS_TICK_SIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(0.35)
    ax.spines["bottom"].set_color("black")
    ax.axvline(0.0, color="black", linewidth=0.35, zorder=1)
    ax.minorticks_off()
    if title:
        ax.set_title(title, loc="left", fontsize=PANEL_TITLE_SIZE, pad=0.0, y=0.985)


def _draw_summary_bar(ax, segments, title, x_label, fmt_value, x_formatter=None):
    """Draw a thin stacked summary bar for overview panel a."""
    total = sum(segment[1] for segment in segments)
    if title:
        ax.set_title(title, loc="left", fontsize=6.8, pad=0.0, y=0.995)

    if total <= 0:
        ax.text(0.5, 0.5, "No classification data",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        return

    left = 0.0
    bar_center = 0.103
    for segment in segments:
        label, value, color = segment[:3]
        label_value = segment[3] if len(segment) > 3 else value
        if value <= 0:
            continue
        ax.barh(
            bar_center, value, left=left, height=0.034,
            color=color, edgecolor="white", linewidth=0.6,
            zorder=2,
        )
        if value / total >= 0.08:
            ax.text(
                left + value / 2.0, bar_center, fmt_value(label_value),
                ha="center", va="center",
                fontsize=ANNOTATION_TEXT_SIZE, color="#222222", zorder=3,
            )
        left += value

    ax.set_xlim(0, total)
    ax.set_ylim(0.055, 0.142)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=2.0, width=0.45, pad=0.6, labelsize=AXIS_TICK_SIZE)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#c7c7c7")
    ax.spines["bottom"].set_linewidth(0.45)
    ax.minorticks_off()
    if x_formatter is not None:
        ax.xaxis.set_major_formatter(x_formatter)
    ax.set_xlabel(x_label, fontsize=AXIS_LABEL_SIZE, labelpad=0.15)


def _draw_dual_summary_panel(ax, all_segments, telo_segments, title=None):
    """Draw a single two-row stacked summary panel with a shared percentage axis."""
    if title:
        ax.set_title(title, loc="left", fontsize=PANEL_TITLE_SIZE, pad=0.0, y=0.985)

    plot_rows = [
        (0.70, "All", all_segments),
        (0.54, "With telomeres", telo_segments),
    ]
    for y_pos, _, segments in plot_rows:
        left = 0.0
        for segment in segments:
            _, value, color = segment[:3]
            cnt = segment[3] if len(segment) > 3 else None
            if value <= 0:
                continue
            ax.barh(
                y_pos, value, left=left, height=0.060,
                color=color, edgecolor="white", linewidth=0.6,
                zorder=2,
            )
            if cnt is not None and value >= 5.0:
                ax.text(
                    left + value / 2.0, y_pos, str(cnt),
                    ha="center", va="center",
                    fontsize=ANNOTATION_TEXT_SIZE, color="#222222", zorder=3,
                )
            left += value

    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(0.43, 0.81)
    ax.set_yticks([row[0] for row in plot_rows])
    ax.set_yticklabels([row[1] for row in plot_rows], fontsize=AXIS_TICK_SIZE)
    ax.tick_params(axis="y", length=0, pad=3.2)
    ax.tick_params(axis="x", length=2.0, width=0.45, pad=0.55, labelsize=AXIS_TICK_SIZE)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}%"))
    ax.set_xlabel("Sequences (%)", fontsize=AXIS_LABEL_SIZE, labelpad=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("black")
    ax.spines["bottom"].set_linewidth(0.35)
    ax.minorticks_off()


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
        ax.text(0.5, 0.5, "No flagged scaffolds\n(fragmented / discordant)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999", linespacing=1.4)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        if title:
            ax.set_title(title, loc="left", fontsize=PANEL_TITLE_SIZE, pad=0.0, y=0.985)
        return

    rows = []
    for entry in flagged_top:
        wrapped = _wrap_scaffold_name(entry["chrom"], width=14)
        rows.append({
            "label": wrapped,
            "value": np.log10(entry["size"] + 1.0),
            "color": COLORS.get(entry["category"], "#aaaaaa"),
        })

    _draw_ranked_bar_panel(
        ax,
        rows,
        x_label="Scaffold size (log10 bp)",
        xticks=[0, 2, 4, 6, 8],
        title=title,
    )


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


def plot_overview_page1(classifications, blocks, chrom_sizes):
    """Page 1: Scaffold classification bars + legend + flagged scaffolds panel."""
    FLAGGED_DISTANCE_BP = 1000

    cat_labels = list(classifications.keys())
    cat_counts = [len(v) for v in classifications.values()]
    cat_colors = [COLORS.get(c, "#aaaaaa") for c in cat_labels]
    total = sum(cat_counts)

    block_rows = _compute_block_rows(blocks, chrom_sizes)
    distance_flagged_blocks = [r for r in block_rows if r["distance"] > FLAGGED_DISTANCE_BP]
    flagged_scaffolds = {
        chrom
        for category in FLAGGED_SCAFFOLD_CATEGORIES
        for chrom in classifications.get(category, [])
    }
    flagged_total_scaffolds = len(flagged_scaffolds)
    distance_flagged_block_total = len(distance_flagged_blocks)

    all_len = [r["length"] for r in block_rows if r["label"] in ("p", "q", "b")]
    total_blocks = len(block_rows)
    total_paths = total if total > 0 else max(len(chrom_sizes), len(blocks))
    median_bp = int(round(_median(all_len))) if all_len else None

    # ---- Figure layout: aligned top-row panels with legends below ----
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, 3.55))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.22],
                          width_ratios=[1.0, 0.85], hspace=0.18, wspace=0.28)

    ax_summary = fig.add_subplot(gs[0, 0])
    ax_flagged_scaffolds = fig.add_subplot(gs[0, 1])
    ax_legends = fig.add_subplot(gs[1, :])
    fig.subplots_adjust(left=0.096, right=0.962, top=0.81, bottom=0.08)

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
    _draw_dual_summary_panel(ax_summary, absolute_pct_segments, relative_segments, title=None)

    # ---- Quality legend (5 categories) ----
    quality_handles = [
        Patch(facecolor=COLORS["T2T"],         label="T2T"),
        Patch(facecolor=COLORS["Incomplete"],   label="Incomplete"),
        Patch(facecolor=COLORS["Fragmented"],   label="Fragmented"),
        Patch(facecolor=COLORS["Discordant"],   label="Discordant"),
        Patch(facecolor=COLORS["No telomeres"], label="No telomeres"),
    ]

    # ---- Shared legend strip: gap key + quality legend ----
    gap_handles = [
        Patch(facecolor="#555555", label="Gapless"),
        Patch(facecolor="#cccccc", label="Gapped"),
    ]
    _draw_centered_legend_pair(
        ax_legends,
        gap_handles,
        quality_handles,
        gap=0.016,
        left_legend_kwargs={
            "ncol": 1,
            "fontsize": LEGEND_TEXT_SIZE,
        },
        right_legend_kwargs={
            "ncol": 3,
            "fontsize": LEGEND_TEXT_SIZE,
        },
    )

    # ---- Panel b: Flagged scaffolds (discordant / fragmented by size) ----
    _draw_flagged_scaffolds_panel(ax_flagged_scaffolds, classifications, chrom_sizes, title=None)

    # ---- Titles and panel labels ----
    fig.canvas.draw()
    label_style = dict(fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="bottom", ha="left")
    title_style = dict(fontsize=PANEL_TITLE_SIZE, va="bottom", ha="center")

    panel_label_y = 0.825
    panel_title_y = 0.827

    bbox_a = ax_summary.get_position()
    fig.text(max(0.005, bbox_a.x0 - 0.032), panel_label_y, "a", **label_style)
    fig.text((bbox_a.x0 + bbox_a.x1) / 2.0, panel_title_y, "Scaffold classification", **title_style)
    bbox_b = ax_flagged_scaffolds.get_position()
    fig.text(max(0.005, bbox_b.x0 - 0.032), panel_label_y, "b", **label_style)
    fig.text((bbox_b.x0 + bbox_b.x1) / 2.0, panel_title_y, "Flagged scaffolds", **title_style)

    fig.suptitle("Assembly summary", fontsize=FIGURE_TITLE_SIZE, fontweight="bold",
                 x=0.5, y=0.968, ha="center")
    summary_line = (
        f"Scaffolds (n={total_paths:,}, flagged={flagged_total_scaffolds:,}), "
        f"telomere blocks (n={total_blocks:,}, distance-flagged={distance_flagged_block_total:,}), "
        f"median {median_bp:,} bp"
        if median_bp is not None else
        f"Scaffolds (n={total_paths:,}, flagged={flagged_total_scaffolds:,}), "
        f"telomere blocks (n={total_blocks:,}, distance-flagged={distance_flagged_block_total:,}), "
        "median NA"
    )
    fig.text(0.5, 0.918, summary_line,
             ha="center", va="top", fontsize=FIGURE_SUMMARY_SIZE, color="#444444")

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

    # ---- Figure layout ----
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, 3.15))
    gs = fig.add_gridspec(1, 3, wspace=0.30, width_ratios=[1.0, 1.0, 0.80])
    ax_rain = fig.add_subplot(gs[0, 0])
    ax_scatter = fig.add_subplot(gs[0, 1])
    ax_flagged = fig.add_subplot(gs[0, 2])
    fig.subplots_adjust(left=0.096, right=0.962, top=0.81, bottom=0.14)

    # ---- Panel c: Vertical raincloud ----
    groups = [
        ("p-arm", p_len, COLORS["p"]),
        ("q-arm", q_len, COLORS["q"]),
        ("balanced", b_len, COLORS["b"]),
    ]
    groups = [(lbl, vals, col) for lbl, vals, col in groups if vals]

    if groups:
        rng = np.random.default_rng(0)
        spacing = 0.36
        offsets = np.arange(len(groups), dtype=np.float64)
        offsets -= (len(groups) - 1) / 2.0
        positions = 0.75 + (spacing * offsets)
        for position, (lbl, vals, col) in zip(positions, groups):
            log_vals = np.log10(np.asarray(vals, dtype=np.float64) + 1.0)
            _draw_raincloud_group(ax_rain, log_vals, position, col, rng, vert=True)

        x_labels = [f"{lbl}\n(n={len(vals)})" for lbl, vals, _ in groups]
        ax_rain.set_xticks(positions)
        ax_rain.set_xticklabels(x_labels, fontsize=AXIS_TICK_SIZE)
        ax_rain.tick_params(axis="x", length=2, pad=2.5)
        ax_rain.tick_params(axis="y", length=2, labelsize=AXIS_TICK_SIZE)
        ax_rain.set_xlim(positions.min() - 0.34, positions.max() + 0.30)
        ax_rain.set_box_aspect(1.0)
        ax_rain.set_anchor("S")
        ax_rain.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))

        median_log = np.log10(_median(all_len) + 1.0)
        ax_rain.axhline(median_log, color="black",
                        linestyle=OVERVIEW_DASH_STYLE, linewidth=OVERVIEW_DASH_WIDTH)
        median_trans = transforms.blended_transform_factory(ax_rain.transAxes, ax_rain.transData)
        ax_rain.text(1.01, median_log, "median",
                     transform=median_trans, ha="left", va="center",
                     rotation=90, fontsize=LEGEND_TEXT_SIZE, color="#555555")
        ax_rain.set_ylabel("Telomere length (log10bp)", fontsize=AXIS_LABEL_SIZE)
    else:
        lbl_text = "No labeled telomere blocks" if block_rows else "No telomere blocks"
        ax_rain.text(0.5, 0.5, lbl_text, ha="center", va="center",
                     transform=ax_rain.transAxes, fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999")
        ax_rain.set_ylabel("Telomere length (log10bp)", fontsize=AXIS_LABEL_SIZE)
        ax_rain.set_xticks([])
        ax_rain.set_box_aspect(1.0)
        ax_rain.set_anchor("S")

    # ---- Panel d: Scatter (terminal offset vs length) ----
    scatter_groups = [
        ("p-arm", "p", COLORS["p"]),
        ("q-arm", "q", COLORS["q"]),
        ("balanced", "b", COLORS["b"]),
    ]
    scatter_rows = [row for row in block_rows if row["label"] in {"p", "q", "b"}]
    if scatter_rows:
        max_log_x = 0.0
        for leg_lbl, arm_key, col in scatter_groups:
            field = "label" if arm_key == "b" else "closestEnd"
            arm_rows = [row for row in scatter_rows if row[field] == arm_key]
            if not arm_rows:
                continue
            x = np.log10(np.asarray([row["distance"] for row in arm_rows], dtype=np.float64) + 1.0)
            y = np.asarray([row["length"] / 1e3 for row in arm_rows], dtype=np.float64)
            max_log_x = max(max_log_x, float(np.max(x)))
            ax_scatter.scatter(x, y, s=13, color=col, alpha=0.55,
                               edgecolors="white", linewidths=0.28,
                               rasterized=True, label=leg_lbl)

        ax_scatter.set_xlabel("Distance to end (log10bp)", fontsize=AXIS_LABEL_SIZE)
        ax_scatter.set_ylabel("Telomere length (kbp)", fontsize=AXIS_LABEL_SIZE)
        ax_scatter.set_xticks([0, 1, 2, 3, 4])
        x_right = max(4.12, max_log_x + 0.18)
        ax_scatter.axvspan(3.0, x_right, facecolor="#ededed", linewidth=0, zorder=0)
        ax_scatter.set_xlim(-0.10, x_right)
        y_max_sc = max(float(np.max(np.asarray([row["length"] / 1e3 for row in scatter_rows],
                                                dtype=np.float64))), 0.0)
        ax_scatter.set_ylim(-0.45, max(1.2, y_max_sc * 1.10))
        ax_scatter.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
        ax_scatter.tick_params(axis="both", labelsize=AXIS_TICK_SIZE)
        ax_scatter.set_box_aspect(1.0)
        ax_scatter.set_anchor("S")
        leg_d = ax_scatter.legend(loc="upper right", fontsize=LEGEND_TEXT_SIZE,
                                  frameon=True, facecolor="white", edgecolor="black",
                                  framealpha=1.0, fancybox=False,
                                  handletextpad=0.35, borderaxespad=0.2)
        leg_d.get_frame().set_linewidth(0.4)
        ax_scatter.axvline(3.0, color="black", linewidth=OVERVIEW_DASH_WIDTH,
                           linestyle=OVERVIEW_DASH_STYLE, zorder=1)
    else:
        ax_scatter.text(0.5, 0.5, "No telomere blocks",
                        transform=ax_scatter.transAxes,
                        ha="center", va="center", fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999")
        ax_scatter.set_xlabel("Distance to end (log10bp)", fontsize=AXIS_LABEL_SIZE)
        ax_scatter.set_ylabel("Telomere length (kbp)", fontsize=AXIS_LABEL_SIZE)
        ax_scatter.set_box_aspect(1.0)
        ax_scatter.set_anchor("S")

    # ---- Panel e: Flagged telomeres lollipop ----
    if not block_rows:
        ax_flagged.text(0.5, 0.5, "No telomere blocks",
                        transform=ax_flagged.transAxes,
                        ha="center", va="center", fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999")
        ax_flagged.set_xticks([])
        ax_flagged.set_yticks([])
        for spine in ax_flagged.spines.values():
            spine.set_visible(False)
    elif not flagged_top:
        ax_flagged.text(
            0.5, 0.5,
            f"No flagged blocks\n(distance > {FLAGGED_DISTANCE_BP // 1000} kbp)",
            transform=ax_flagged.transAxes,
            ha="center", va="center", fontsize=PLACEHOLDER_TEXT_SIZE, color="#999999", linespacing=1.4,
        )
        ax_flagged.set_xticks([])
        ax_flagged.set_yticks([])
        for spine in ax_flagged.spines.values():
            spine.set_visible(False)
    else:
        rows = []
        for scaffold, info in flagged_top:
            wrapped = _wrap_scaffold_name(scaffold, width=14)
            dominant_arm = info["arms"].most_common(1)[0][0]
            rows.append({
                "label": f"{wrapped}\n(n={info['count']})",
                "value": np.log10(info["length"] + 1.0),
                "color": COLORS.get(dominant_arm, "#aaaaaa"),
            })

        _draw_ranked_bar_panel(
            ax_flagged,
            rows,
            x_label="Flagged length (log10 bp)",
            xticks=[0, 1, 2, 3, 4],
        )

    ax_flagged.set_anchor("S")

    # ---- Panel labels and titles ----
    fig.canvas.draw()
    label_style = dict(fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="bottom", ha="left")
    title_style = dict(fontsize=PANEL_TITLE_SIZE, ha="center", va="bottom")

    panel_meta = [
        (ax_rain,    "c", "Length by arm"),
        (ax_scatter, "d", "Telomere positioning"),
        (ax_flagged, "e", "Flagged telomeres"),
    ]
    title_y = 0.845
    for ax, lbl, title in panel_meta:
        bbox = ax.get_position()
        fig.text(max(0.005, bbox.x0 - 0.032), title_y, lbl, **label_style)
        fig.text((bbox.x0 + bbox.x1) / 2.0, title_y + 0.002, title, **title_style)

    return fig


# ---------------------------------------------------------------------------
# ITS report pages (rendered last, after every per-chromosome zoom; see main())
# ---------------------------------------------------------------------------

def _class_totals(df):
    """Per-class row count and total bp, reindexed to the fixed CLASS_ORDER."""
    counts = df.groupby("teloType")["teloLen"].size().reindex(CLASS_ORDER, fill_value=0)
    totals_bp = df.groupby("teloType")["teloLen"].sum().reindex(CLASS_ORDER, fill_value=0)
    return counts, totals_bp


def _fmt_frac(value):
    """Format a fraction, or 'NA' when it is undefined."""
    return f"{value:.2f}" if np.isfinite(value) else "NA"


def _fmt_log_bp_tick(value):
    """Format a signed log10(bp) decade tick (value=0 at the scaffold end) as clean bp/kb/Mb text."""
    if value == 0:
        return "0"
    bp = 10 ** abs(value)
    divisor, unit = _pick_bp_unit(bp)
    return f"{bp / divisor:g} {unit}"


def _short_scaffold_labels(names, width=18):
    """Strip a shared prefix (common in PanSN-style names) then clip for display.

    Rows need the *distinguishing* suffix, not a sample/haplotype tag every row
    shares, so this trims to the last separator inside the common prefix.
    """
    if len(set(names)) > 1:
        common = os.path.commonprefix(names)
        cut = 0
        for i in range(len(common), 0, -1):
            if common[i - 1] in "#_./:|-":
                cut = i
                break
        names = [n[cut:] or n for n in names]
    return [n if len(n) <= width else n[:width - 1] + "…" for n in names]


def _its_atlas_chroms(df, arm_blocks, chrom_sizes):
    """Scaffolds carrying a telomere or ITS row, longest first (ideogram row order)."""
    telomere_chroms = {c for c, blist in arm_blocks.items() if blist}
    its_chroms = set(df["chr"].unique()) if not df.empty else set()
    return sorted(telomere_chroms | its_chroms, key=lambda c: (-chrom_sizes.get(c, 0), c))


def _split_long_short_chroms(atlas_chroms, chrom_sizes):
    """Split scaffolds into >=20%-of-longest and shorter, preserving length-descending order."""
    if not atlas_chroms:
        return [], []
    cutoff = chrom_sizes.get(atlas_chroms[0], 0) * 0.20
    long_chroms = [c for c in atlas_chroms if chrom_sizes.get(c, 0) >= cutoff]
    short_chroms = [c for c in atlas_chroms if chrom_sizes.get(c, 0) < cutoff]
    return long_chroms, short_chroms


def _nearest_end(pos, chrom_size):
    """Return ('p'|'q', distance) for whichever scaffold end pos sits closest to."""
    dist_p, dist_q = pos, max(chrom_size - pos, 0)
    return ("p", dist_p) if dist_p <= dist_q else ("q", dist_q)


# ---------------------------------------------------------------------------
# Region-zoom toolkit shared with plot_its.py (single-locus figures)
# ---------------------------------------------------------------------------

def _region_axis_spec(view_start, view_end):
    """Round genomic ticks for an absolute-position window; positions always in Mbp."""
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


def _draw_locus_ideogram(ax, chrom_size, view_start, view_end, blocks_list, its_blocks_list,
                         contig_blocks_list=None, label="Scaffold"):
    """Whole-scaffold overview for a single-locus zoom: terminal caps, all ITS ticks, zoom box."""
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
        labels = np.array([b.get("label", "") for b in its_blocks_list])
        for k in np.unique(labels):
            mask = labels == k
            segs = np.empty((int(mask.sum()), 2, 2))
            segs[:, 0, 0] = xs[mask]; segs[:, 0, 1] = bottom
            segs[:, 1, 0] = xs[mask]; segs[:, 1, 1] = bottom + bar_h
            ax.add_collection(LineCollection(
                segs, colors=ITS_ORIENT_COLORS.get(k, COLORS["its"]),
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
                           gap_blocks_list, label="Blocks", contig_blocks_list=None):
    """Single-panel block track in genomic coordinates; ITS colored by orientation.

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
    _draw(its_blocks_list, lambda b: ITS_ORIENT_COLORS.get(b.get("label", ""), COLORS["its"]), 4)
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

def _its_orientation(df):
    return df["teloLabel"].where(df["teloLabel"].isin(("p", "q", "b")), "unknown")


def _scan_note(params):
    mode = params.get("ultra_fast")
    if mode is False:
        return "Full scan; observed ITS rows"
    if mode is True:
        limit = params.get("terminal_limit")
        return ("End scan (initial " + _fmt_bp(limit) + " per end); observed ITS rows"
                if limit else "End scan; scan limit unknown; observed ITS rows")
    return "Scan scope unknown; observed ITS rows"


def _its_header(fig, title, subtitle):
    fig.suptitle(title, fontsize=PANEL_LABEL_SIZE, fontweight="bold", y=0.96)
    fig.text(0.5, 0.895, subtitle, ha="center", va="top",
             fontsize=FIGURE_SUMMARY_SIZE, color="black")


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


def _its_position_counts(df, chrom_sizes):
    """One count per row midpoint in 100 equal-width bins per scaffold."""
    out = {}
    for chrom, grp in df.groupby("chr", sort=False):
        extent = max(int(chrom_sizes.get(chrom, 0)), int(grp["end"].max()), 1)
        edges = np.linspace(0, extent, ITS_POSITION_BINS + 1)
        mid = grp["start"].to_numpy(dtype=float) + (
            grp["end"].to_numpy(dtype=float) - grp["start"].to_numpy(dtype=float)) / 2
        out[chrom] = (edges, np.histogram(mid, bins=edges)[0])
    return out


def _its_top_hits(df, pairs, clusters):
    """Scaffold -> which top-ranked ITS list(s) (L1/C1/F1) it heads, for atlas star labels."""
    top = {}
    for prefix, frame in (("L1", rank_long_its(df, 1)), ("C1", clusters.head(1)), ("F1", pairs.head(1))):
        if not frame.empty:
            top.setdefault(frame.iloc[0]["chr"], []).append(prefix)
    return top


def plot_its_overview_page(df, pairs, clusters, arm_blocks, chrom_sizes, params,
                           chroms=None, page_number=1, page_count=1, position_counts=None,
                           all_chroms=None, top_hits=None, known_sizes=None):
    """Fixed-size atlas page; count bins include every observed midpoint once."""
    all_chroms = _its_atlas_chroms(df, arm_blocks, chrom_sizes) if all_chroms is None else all_chroms
    chroms = all_chroms[:ITS_ATLAS_ROWS] if chroms is None else list(chroms)
    counts = _its_position_counts(df, chrom_sizes) if position_counts is None else position_counts
    labels = df.attrs.get("display_labels") or _its_labels(df, chrom_sizes, arm_blocks)
    vmax = max((int(v.max()) for _, v in counts.values()), default=1)
    norm = LogNorm(vmin=1, vmax=max(2, vmax))
    cmap = LinearSegmentedColormap.from_list("its_counts", ["#DCECF5", COLORS["p"]])
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, REPORT_PAGE_HEIGHT))
    gs = fig.add_gridspec(1, 2, left=0.18, right=0.96, bottom=0.22, top=0.78, wspace=0.85)
    top = _its_top_hits(df, pairs, clusters) if top_hits is None else top_hits
    known = (df.groupby("chr")["chrSize"].max().to_dict() if known_sizes is None else known_sizes)
    chunks = [chroms[:10], chroms[10:]]
    for col, names in enumerate(chunks):
        ax = fig.add_subplot(gs[0, col])
        if not names:
            ax.set_axis_off()
            continue
        # Each row uses a normalized 0-1 axis (own physical scale, shown in the label) so microchromosomes don't disappear beside large scaffolds.
        for y, chrom in enumerate(names):
            extent = max(int(chrom_sizes.get(chrom, 0)), 1)
            ax.plot([0, 1], [y, y], color=COLORS["gap"], lw=7, solid_capstyle="butt", zorder=0)
            if (params.get("ultra_fast") is True and params.get("terminal_limit") and
                    known.get(chrom, extent) > 0):
                fraction = min(params["terminal_limit"] / extent, 0.5)
                ax.add_patch(Rectangle((fraction, y - 0.32), 1 - 2 * fraction, 0.64,
                                       facecolor="white", edgecolor="none", zorder=1))
            if chrom in counts:
                edges, values = counts[chrom]
                occupied = np.flatnonzero(values)
                patches = [Rectangle((edges[i] / extent, y - 0.30),
                                     (edges[i + 1] - edges[i]) / extent, 0.60) for i in occupied]
                if patches:
                    collection = PatchCollection(patches, cmap=cmap, norm=norm,
                                                 edgecolors="none", rasterized=True, zorder=2)
                    collection.set_array(values[occupied])
                    ax.add_collection(collection)
            for block in arm_blocks.get(chrom, []):
                side = block.get("closestEnd")
                if side in ("p", "q"):
                    ax.plot([0 if side == "p" else 1], [y], marker="|",
                            color=COLORS["terminal"], ms=8, mew=1, zorder=3)
        tick_labels = []
        for c in names:
            suffix = " " + "/".join(top[c]) if c in top else ""
            extent_label = _fmt_bp(chrom_sizes.get(c, 0))
            if c in known and known[c] == 0:
                extent_label += "*"
            tick_labels.append(f"{labels[c]}{suffix}\n{extent_label}")
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(tick_labels, fontsize=MIN_TEXT_SIZE)
        ax.set_ylim(9.6, -0.7)
        ax.set_xlim(-0.01, 1.01)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "50", "100"])
        ax.set_xlabel("Position within plotted extent (%)")
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    cax = fig.add_axes([0.40, 0.08, 0.25, 0.025])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation="horizontal")
    cb.set_ticks(sorted(set([1, max(1, vmax)])))
    cb.set_ticklabels([str(v) for v in sorted(set([1, max(1, vmax)]))])
    cb.ax.minorticks_off()
    fig.text(0.39, 0.095, "Rows per 1% bin", ha="right", va="center", fontsize=MIN_TEXT_SIZE)
    fig.text(0.98, 0.022, "Midpoint counts; white: outside initial end windows; *extent inferred; caps: terminal calls",
             ha="right", fontsize=MIN_TEXT_SIZE)
    _its_header(fig, f"Interstitial telomeres: atlas ({page_number}/{page_count})",
                f"{_scan_note(params)}; {len(chroms)} of {len(all_chroms)} scaffolds on this page")
    return fig


def _draw_its_length_hist_panel(ax, df):
    """Exact ECDFs, including ties, singleton groups and unknown orientations."""
    styles = ["-", "--", ":", "-."]
    orientation = _its_orientation(df)
    for (key, name), style in zip(ITS_ORIENT_LABELS, styles):
        values = df.loc[orientation == key, "teloLen"].to_numpy(dtype=float)
        if not len(values):
            continue
        x, counts = np.unique(values, return_counts=True)
        cumulative = counts.cumsum() / len(values)
        ax.step(np.r_[x[0], x], np.r_[0, cumulative], where="post",
                color=ITS_ORIENT_COLORS[key], linestyle=style, linewidth=1,
                label=f"{name} (n={len(values):,})")
        if len(x) == 1:
            ax.plot(x, [1], marker="o", ms=2, color=ITS_ORIENT_COLORS[key])
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("ITS row length (bp)")
    ax.set_ylabel("Cumulative fraction")
    if len(df):
        ax.legend(loc="lower right", fontsize=MIN_TEXT_SIZE)
    ax.xaxis.set_minor_locator(ticker.NullLocator())


def _draw_its_class_totals_panel(ax_count, ax_mb, df):
    classes = df["teloType"].where(df["teloType"].isin(CLASS_ORDER), "other")
    order = CLASS_ORDER + (["other"] if (classes == "other").any() else [])
    counts = df.groupby(classes)["teloLen"].size().reindex(order, fill_value=0)
    totals = df.groupby(classes)["teloLen"].sum().reindex(order, fill_value=0)
    for ax, series, label, fmt in ((ax_count, counts, "Rows", lambda x: f"{int(x):,}"),
                                   (ax_mb, totals, "Summed row length (bp)", _fmt_bp)):
        vals = series.to_numpy(dtype=float)
        ax.set_xscale("log")
        positive = vals[vals > 0]
        low = min(positive.min() / 2, 0.5) if len(positive) else 0.5
        high = max(vals.max(), 1) * 20
        ax.set_xlim(low, high)
        for i, (key, val) in enumerate(zip(order, vals)):
            if val > 0:
                ax.plot(val, i, marker="o", ms=3, color=CLASS_COLORS.get(key, COLORS["its"]))
            ax.text(max(val, low) * 1.3, i, fmt(val), va="center", fontsize=MIN_TEXT_SIZE)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order if ax is ax_count else [])
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.set_xlabel(label)
        ax.xaxis.set_major_locator(ticker.LogLocator(numticks=3))
        ax.xaxis.set_minor_locator(ticker.NullLocator())
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)


def _draw_its_composition_panel(ax, df, motif_len, min_canonical_count):
    """Positive canonical estimates at their true coordinates; excluded zeros counted."""
    valid = (df["teloLen"] > 0) & (df["canonical_bp"] > 0)
    positive = df.loc[valid]
    excluded = len(df) - len(positive)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("ITS row length (bp)")
    ax.set_ylabel("Canonical match bp")
    if positive.empty:
        ax.set_xlim(1, 10)
        ax.set_ylim(1, 10)
        ax.text(0.5, 0.5, "No positive canonical counts", transform=ax.transAxes,
                ha="center", fontsize=PLACEHOLDER_TEXT_SIZE)
    else:
        x = positive["teloLen"].to_numpy(dtype=float)
        y = positive["canonical_bp"].to_numpy(dtype=float)
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        ax.set_xlim(x.min() / 1.5, x.max() * 1.5)
        ax.set_ylim(y.min() / 1.5, y.max() * 1.5)
        if len(positive) > 500:
            hb = ax.hexbin(x, y, xscale="log", yscale="log", gridsize=35, mincnt=1,
                           cmap=LinearSegmentedColormap.from_list("composition", ["#DCECF5", COLORS["p"]]),
                           norm=LogNorm(vmin=1, vmax=max(2, len(positive))), linewidths=0,
                           rasterized=True)
            # Normalize to actual occupied-bin counts (not the number of input rows).
            hb.set_clim(1, max(2, float(hb.get_array().max())))
            cb = ax.figure.colorbar(hb, ax=ax, fraction=0.045, pad=0.025)
            cb.set_label("Rows per hexagon (log)", fontsize=MIN_TEXT_SIZE)
            cb.ax.tick_params(labelsize=MIN_TEXT_SIZE)
        else:
            orientation = _its_orientation(positive)
            for (key, name), marker in zip(ITS_ORIENT_LABELS, ("<", ">", "D", "o")):
                mask = orientation == key
                if mask.any():
                    ax.scatter(x[mask], y[mask], s=9, color=ITS_ORIENT_COLORS[key],
                               edgecolors="black", linewidths=0.2, marker=marker, label=name)
            ax.legend(loc="upper left", fontsize=MIN_TEXT_SIZE, handletextpad=0.2)
        ax.plot([lo / 2, hi * 2], [lo / 2, hi * 2], color="black", lw=0.5,
                linestyle="--", zorder=0)
        if min_canonical_count is not None and min_canonical_count > 0:
            floor = min_canonical_count * motif_len
            ax.axhline(floor, color="black", lw=0.5, linestyle=":", zorder=0)
    ax.text(0.02, -0.29, f"Positive n={len(positive):,}; zero canonical n={excluded:,}",
            transform=ax.transAxes, fontsize=MIN_TEXT_SIZE)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())


def plot_its_statistics_page(df, pairs, clusters, params):
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, REPORT_PAGE_HEIGHT))
    gs = fig.add_gridspec(2, 2, left=0.085, right=0.97, top=0.77, bottom=0.22,
                         width_ratios=[1.05, 1], hspace=1.0, wspace=0.55)
    comp = fig.add_subplot(gs[:, 0])
    ecdf = fig.add_subplot(gs[0, 1])
    sub = gs[1, 1].subgridspec(1, 2, wspace=0.55)
    count, length = fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])
    known = params.get("known_params", set())
    floor = params.get("min_canonical_count") if {"motif_len", "min_canonical_count"} <= known else None
    _draw_its_composition_panel(comp, df, params.get("motif_len", 6), floor)
    _draw_its_length_hist_panel(ecdf, df)
    _draw_its_class_totals_panel(count, length, df)
    _its_header(fig, "Interstitial telomeres: distributions",
                f"{_scan_note(params)}; n={len(df):,} rows; {len(pairs):,} candidate pairs; {len(clusters):,} clusters")
    motif_note = "" if "motif_len" in known else " (assumed)"
    note = f"Canonical match bp = count × {params.get('motif_len', 6)} bp{motif_note}; dashed: equal bp"
    if floor is not None:
        note += "; dotted: configured count floor"
    fig.text(0.5, 0.027, note, ha="center", fontsize=MIN_TEXT_SIZE)
    return fig


def _table_bbox(n_data_rows, ax_height_in):
    """Top-aligned bbox giving every table row a fixed physical height instead of stretching."""
    frac = min(1.0, (n_data_rows + 1) * TABLE_ROW_IN / ax_height_in)
    return [0.0, 1.0 - frac, 1.0, frac]


def _style_table(tab):
    """Compact table style: small text, thin horizontal-only rules, no vertical lines."""
    tab.auto_set_font_size(False)
    tab.set_fontsize(TABLE_TEXT_SIZE)
    for cell in tab.get_celld().values():
        cell.set_linewidth(0.3)
        cell.set_edgecolor("#cccccc")
        cell.visible_edges = "horizontal"


def _draw_its_summary_tables(ax_long, ax_cluster, ax_fusion, df, clusters, pairs, ax_height_in):
    """Panels b-d: Long ITS / Clusters / Candidate fusions, compact ax.table renders."""
    for ax in (ax_long, ax_cluster, ax_fusion):
        ax.axis("off")

    ax_long.set_title("Long ITS, ranked by canonical bp", loc="left", fontsize=PANEL_TITLE_SIZE, pad=2)
    top_long = rank_long_its(df, 5)
    if top_long.empty:
        ax_long.text(0.02, 0.9, "none", fontsize=TABLE_TEXT_SIZE, va="top", color="black")
    else:
        short_names = [df.attrs.get("display_labels", {}).get(c, c) for c in top_long["chr"]]
        rows = [[f"L{i + 1}", name, f"{(row.start + row.end) / 2e6:.3f}", _fmt_bp(int(row.teloLen)),
                f"{int(row.canonical_bp):,}", _fmt_frac(row.can_prop),
                row.teloLabel if row.teloLabel in ("p", "q", "b") else "other",
                CLASS_SHORT.get(row.teloType, "other")]
               for i, (row, name) in enumerate(zip(top_long.itertuples(index=False), short_names))]
        tab = ax_long.table(
            cellText=rows,
            colLabels=["#", "scaffold", "pos (Mb)", "length", "can. bp", "ratio", "strand", "class"],
            colWidths=[0.08, 0.23, 0.13, 0.13, 0.12, 0.11, 0.10, 0.10],
            loc="upper left", cellLoc="left", colLoc="left",
            bbox=_table_bbox(len(rows), ax_height_in),
        )
        _style_table(tab)

    ax_cluster.set_title("ITS clusters (≥3 rows within 50 kb)", loc="left",
                         fontsize=PANEL_TITLE_SIZE, pad=2)
    top_clusters = clusters.head(5)
    if top_clusters.empty:
        ax_cluster.text(0.02, 0.9, "none", fontsize=TABLE_TEXT_SIZE, va="top", color="black")
    else:
        short_names = [df.attrs.get("display_labels", {}).get(c, c) for c in top_clusters["chr"]]
        rows = [[f"C{i + 1}", name, f"{row.start / 1e6:.3f}-{row.end / 1e6:.3f}", f"{row.rows}",
                _fmt_bp(int(row.span)), _fmt_bp(int(row.its_bp)), _fmt_bp(int(row.canonical_bp))]
               for i, (row, name) in enumerate(zip(top_clusters.itertuples(index=False), short_names))]
        tab = ax_cluster.table(
            cellText=rows,
            colLabels=["#", "scaffold", "start-end (Mb)", "rows", "span", "ITS bp", "can. bp"],
            colWidths=[0.07, 0.19, 0.25, 0.09, 0.14, 0.13, 0.13],
            loc="upper left", cellLoc="left", colLoc="left",
            bbox=_table_bbox(len(rows), ax_height_in),
        )
        _style_table(tab)

    ax_fusion.set_title("Candidate fusions (q→p), ranked by shorter arm", loc="left",
                        fontsize=PANEL_TITLE_SIZE, pad=2)
    top_pairs = pairs.head(5)
    if top_pairs.empty:
        ax_fusion.text(0.02, 0.9, "none", fontsize=TABLE_TEXT_SIZE, va="top", color="black")
    else:
        short_names = [df.attrs.get("display_labels", {}).get(c, c) for c in top_pairs["chr"]]
        rows = [[f"F{i + 1}", name, f"{(row.start + row.end) / 2e6:.3f}",
                f"{row.q_bp:,}", f"{row.p_bp:,}", f"{row.spacer_bp:,}",
                f"{_fmt_frac(row.q_can_prop)}/{_fmt_frac(row.p_can_prop)}"]
               for i, (row, name) in enumerate(zip(top_pairs.itertuples(index=False), short_names))]
        tab = ax_fusion.table(
            cellText=rows,
            colLabels=["#", "scaffold", "pos (Mb)", "q bp", "p bp", "spacer", "ratio q/p"],
            colWidths=[0.06, 0.25, 0.13, 0.13, 0.13, 0.11, 0.19],
            loc="upper left", cellLoc="left", colLoc="left",
            bbox=_table_bbox(len(rows), ax_height_in),
        )
        _style_table(tab)


def plot_its_composition_page(df, pairs, clusters, params):
    """Compact ranked tables; population statistics live on their own page."""
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, REPORT_PAGE_HEIGHT))
    gs = fig.add_gridspec(3, 1, left=0.06, right=0.98, top=0.77, bottom=0.13, hspace=0.42)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(3)]
    height = axes[0].get_position().height * REPORT_PAGE_HEIGHT
    _draw_its_summary_tables(*axes, df, clusters, pairs, height)
    _its_header(fig, "Interstitial telomeres: selected candidates",
                f"Top 5 per table; {len(df):,} observed rows; {len(clusters):,} clusters; {len(pairs):,} candidate pairs")
    distance = params.get("max_block_dist", 1000)
    assumed = "" if "max_block_dist" in params.get("known_params", set()) else " (assumed)"
    fig.text(0.5, 0.028, f"Pairs: q→p, distance ≤{distance:,} bp{assumed}, engine fusion label; candidates are not confirmed fusions.",
             ha="center", fontsize=MIN_TEXT_SIZE)
    return fig

# ---------------------------------------------------------------------------
# ITS-3 "Top loci"
# ---------------------------------------------------------------------------

def resolve_its_loci(clusters, df, pairs, chrom_sizes):
    """Resolve up to 3 (label, chrom, start, end) windows for C1, L1, F1; skip missing ones."""
    loci = []
    if not clusters.empty:
        c0 = clusters.iloc[0]
        s, e = pad_window(int(c0["start"]), int(c0["end"]), chrom_sizes.get(c0["chr"], c0["end"]))
        loci.append(("C1", c0["chr"], s, e))
    long_its = rank_long_its(df, 1)
    if not long_its.empty:
        l0 = long_its.iloc[0]
        s, e = pad_window(int(l0["start"]), int(l0["end"]), chrom_sizes.get(l0["chr"], l0["end"]))
        loci.append(("L1", l0["chr"], s, e))
    if not pairs.empty:
        f0 = pairs.iloc[0]
        s, e = pad_window(int(f0["start"]), int(f0["end"]), chrom_sizes.get(f0["chr"], f0["end"]))
        loci.append(("F1", f0["chr"], s, e))
    return loci


def plot_its_loci_page(loci, chrom_sizes, arm_blocks, its_blocks, gap_blocks,
                       density_data, canonical_data, strand_data, gc_data, entropy_data,
                       unknown_extents=(), display_labels=None):
    """ITS-3 'Top loci': one terminal-style track-stack column per locus (C1, L1, F1)."""
    if not loci:
        return _placeholder_figure("Interstitial telomeres: top loci", "No ITS loci to show.")

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
    fig_h = REPORT_PAGE_HEIGHT
    fig, axes = plt.subplots(
        len(height_ratios), n_cols, figsize=(FIG_WIDTH_DOUBLE, fig_h),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.35,
                    "wspace": 0.30 if n_cols > 1 else 0.0},
        squeeze=False)
    fig.subplots_adjust(left=0.145, right=0.97, top=0.76, bottom=0.16)

    for col, (label, chrom, start, end) in enumerate(loci):
        size = int(chrom_sizes.get(chrom, end))
        blist = arm_blocks.get(chrom, [])
        itslist = its_blocks.get(chrom, []) if its_blocks else []
        gaplist = gap_blocks.get(chrom, []) if gap_blocks else []
        show_label = col == 0

        ideo_ax = axes[0][col]
        box_l, box_r, box_bottom = _draw_locus_ideogram(
            ideo_ax, size, start, end, blist, itslist,
            label=("Observed\nextent*" if chrom in unknown_extents else "Scaffold") if show_label else None)
        if chrom in unknown_extents:
            ideo_ax.set_xlabel("Observed extent (Mbp); scaffold length unknown", fontsize=AXIS_TICK_SIZE)

        axis_spec = _region_axis_spec(start, end)
        visible_rows = []
        for i, (name, _) in enumerate(track_specs):
            row = i + 1
            ax = axes[row][col]
            if name == "blocks":
                visible = _draw_its_blocks_track(ax, start, end, blist, itslist, gaplist,
                                                 label="Blocks" if show_label else None)
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
                                               y_max=2.0, y_min=1.5)
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

        target_row = visible_rows[0] if visible_rows else 1
        target_ax = axes[target_row][col]
        for x_ideo, x_ax in ((box_l, 0.0), (box_r, 1.0)):
            fig.add_artist(ConnectionPatch(
                xyA=(x_ideo, box_bottom), coordsA=ideo_ax.transData,
                xyB=(x_ax, 1.0), coordsB=target_ax.transAxes,
                color=ZOOM_COLOR, linewidth=0.5, linestyle=(0, (3, 2)), alpha=0.55, zorder=0))

        short = (display_labels or {}).get(chrom, _short_scaffold_labels([chrom], width=20)[0])
        ideo_ax.set_title(f"{label}   {short}", fontsize=PANEL_TITLE_SIZE, pad=3, loc="center")

    fig.suptitle("Interstitial telomeres: top loci", fontsize=FIGURE_TITLE_SIZE,
                fontweight="bold", y=1 - 0.20 / fig_h, x=0.5, ha="center")
    fig.text(0.5, 1 - 0.44 / fig_h,
            "   |   ".join(f"{lbl}: {s:,}-{e:,} bp" for lbl, chrom, s, e in loci),
            ha="center", va="top", fontsize=FIGURE_SUMMARY_SIZE, color="black")

    label_style = dict(fontsize=PANEL_LABEL_SIZE, fontweight="bold", va="bottom", ha="left")
    for col in range(n_cols):
        bbox = axes[0][col].get_position()
        fig.text(max(0.005, bbox.x0 - 0.02), bbox.y1 + 0.030, "abc"[col], **label_style)
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
        for start, end in zip(starts[no_data], ends[no_data]):
            ax.axvspan(start, end, ymin=0, ymax=1,
                       facecolor=COLORS["no_data"], alpha=0.22,
                       linewidth=0, zorder=0)

    for run_start, run_end in _iter_true_runs(valid):
        xs, ys = _bedgraph_to_step(
            starts[run_start:run_end],
            ends[run_start:run_end],
            values[run_start:run_end],
        )
        ax.fill_between(xs, ys, step="post", color=color,
                        alpha=0.16, linewidth=0)
        ax.plot(xs, ys, drawstyle="steps-post", color=color,
                linewidth=0.6)

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

    print(f"Found files: {', '.join(files.keys())}", file=sys.stderr)

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

    density_data = parse_bedgraph(files["density"]) if "density" in files else None
    canonical_data = parse_bedgraph(files["canonical_ratio"]) if "canonical_ratio" in files else None
    strand_data  = parse_bedgraph(files["strand_ratio"]) if "strand_ratio" in files else None
    its_blocks = parse_terminal_bed(files["interstitial"]) if "interstitial" in files else None
    gap_blocks = parse_interval_bed(files["gaps"]) if "gaps" in files else None
    gc_data = parse_bedgraph(files["gc"]) if "gc" in files else None
    entropy_data = parse_bedgraph(files["entropy"]) if "entropy" in files else None

    chrom_sizes = get_chrom_sizes(blocks, density_data, canonical_data, strand_data)

    classifications = parse_report(files["report"]) if "report" in files else OrderedDict()
    # a flagged scaffold appears under its completeness class and its anomaly, so count it once
    total_chroms = len({chrom for v in classifications.values() for chrom in v})
    total_telo = sum(len(blist) for blist in arm_blocks.values())
    cat_summary = ", ".join(f"{k}={len(v)}" for k, v in classifications.items()) or "none"
    print(f"Chromosomes: {total_chroms}  |  Telomere blocks: {total_telo}  |  "
          f"Categories: {cat_summary}",
          file=sys.stderr)

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
            print(f"ITS rows: {len(its_frame)}  |  Candidate fusions: {len(pairs)}  |  "
                  f"Clusters: {len(clusters)}", file=sys.stderr)
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
            "its-statistics", "teloscope_its_statistics.png",
            lambda: plot_its_statistics_page(its_frame, pairs, clusters, params),
            "ITS distributions", "Failed to render ITS distributions.",
        ))
        # Atlas geometry and rankings are dataset-wide; compute once here rather than per page.
        atlas = _its_atlas_chroms(its_frame, arm_blocks, chrom_sizes)
        position_counts = _its_position_counts(its_frame, chrom_sizes)
        top_hits = _its_top_hits(its_frame, pairs, clusters)
        known_sizes = its_frame.groupby("chr")["chrSize"].max().to_dict()
        atlas_count = max(1, (len(atlas) + ITS_ATLAS_ROWS - 1) // ITS_ATLAS_ROWS)
        for page_idx in range(atlas_count):
            chroms = atlas[page_idx * ITS_ATLAS_ROWS:(page_idx + 1) * ITS_ATLAS_ROWS]
            pages.append((
                f"its-atlas-{page_idx + 1}", f"teloscope_its_atlas_{page_idx + 1:03d}.png",
                lambda chroms=chroms, page_idx=page_idx: plot_its_overview_page(
                    its_frame, pairs, clusters, arm_blocks, chrom_sizes, params,
                    chroms, page_idx + 1, atlas_count, position_counts, atlas, top_hits, known_sizes),
                f"ITS atlas {page_idx + 1}/{atlas_count}", "Failed to render ITS atlas.",
            ))
        pages.append((
            "its-candidates", "teloscope_its_candidates.png",
            lambda: plot_its_composition_page(its_frame, pairs, clusters, params),
            "ITS selected candidates", "Failed to render ITS candidates.",
        ))
        for locus in loci:
            pages.append((
                f"its-locus-{locus[0]}", f"teloscope_its_locus_{locus[0]}.png",
                lambda locus=locus: plot_its_loci_page(
                    [locus], chrom_sizes, arm_blocks, its_blocks, gap_blocks,
                    density_data, canonical_data, strand_data, gc_data, entropy_data,
                    unknown_extents, its_frame.attrs["display_labels"]),
                f"ITS selected locus {locus[0]}", "Failed to render ITS locus.",
            ))

    n_figures = len(pages)

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
        print(f"Figures saved to {out_dir}/", file=sys.stderr)
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
            print(f"Report saved to {pdf_path}", file=sys.stderr)

    if fallback_pages:
        preview = ", ".join(f"{name} ({err})" for name, err in fallback_pages[:5])
        more = " ..." if len(fallback_pages) > 5 else ""
        _warn(f"Report generation completed with {len(fallback_pages)} placeholder figure(s): {preview}{more}")


if __name__ == "__main__":
    main()
