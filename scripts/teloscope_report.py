#!/usr/bin/env python3
"""
teloscope_report.py — Publication-ready figures from Teloscope output.

Usage:
    python teloscope_report.py <output_directory> [-o report.pdf]
    python teloscope_report.py <output_directory> --png

Reads Teloscope output files from the given directory and generates:
  Page 1: Assembly overview (classification summary + telomere length distributions)
  Page 2+: Per-chromosome terminal zoom figures (blocks, density, canonical ratio, strand bias)

Requires: Python 3.6+, matplotlib, numpy, pandas
"""

import sys
import os
import glob
import re
import argparse
import textwrap
from collections import Counter, defaultdict, OrderedDict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from matplotlib import transforms

# ---------------------------------------------------------------------------
# Nature-style configuration
# ---------------------------------------------------------------------------

# Nature figure widths: 89 mm (single column), 183 mm (double column)
FIG_WIDTH_SINGLE = 3.50    # inches (89 mm)
FIG_WIDTH_DOUBLE = 7.20    # inches (183 mm)

COLORS = {
    # Classification palette (Nature-compatible, colorblind-friendly)
    "T2T":                 "#1A9641",
    "Gapped T2T":          "#9CCF60",
    "Incomplete":          "#FFC754",
    "Gapped Incomplete":   "#FFE885",
    "Misassembly":         "#D6594C",
    "Gapped Misassembly":  "#F58B6D",
    "Discordant":          "#8278F4",
    "Gapped Discordant":   "#B395EB",
    "No telomeres":        "#C8C8C8",
    "Gapped No telomeres": "#F0F0F0",
    # Arm / track colors
    "p":            "#0072B2",
    "q":            "#D55E00",
    "b":            "#FFC754",
    "density":      "#2A9D59",
    "canonical":    "#1B9ECA",
    "strand_bias":  "#8D80C6",
    "no_data":      "#e0e0e0",
}

OVERVIEW_DASH_STYLE = (0, (2.2, 2.2))
OVERVIEW_DASH_WIDTH = 0.35

def _apply_nature_style():
    """Apply Nature journal rcParams globally."""
    plt.rcParams.update({
        # Fonts — Nature requires sans-serif (Arial / Helvetica)
        "font.family":        "sans-serif",
        "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":          7,
        # Axes
        "axes.titlesize":     8,
        "axes.labelsize":     7,
        "axes.linewidth":     0.5,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        # Ticks
        "xtick.labelsize":    6,
        "ytick.labelsize":    6,
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
        "legend.fontsize":    6,
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

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def find_files(directory):
    """Auto-detect Teloscope output files in *directory*."""
    files = {}
    patterns = {
        "terminal":        "*_terminal_telomeres.bed",
        "interstitial":    "*_interstitial_telomeres.bed",
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
    Parse *_terminal_telomeres.bed -> dict[chrom -> list of block dicts].
    Columns: chrom start end length label fwdCount revCount canCount nonCanCount pathSize terminality
    """
    blocks = defaultdict(list)
    malformed = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: expected at least 10 BED columns, found {len(parts)}; skipping.")
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                length = int(parts[3])
                fwd = int(parts[5])
                rev = int(parts[6])
                can = int(parts[7])
                noncan = int(parts[8])
                path_size = int(parts[9]) if parts[9] else 0
            except ValueError as exc:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid numeric field ({exc}); skipping.")
                continue

            if end <= start or length <= 0:
                malformed += 1
                if malformed <= 3:
                    _warn(f"{path}:{lineno}: invalid interval start={start} end={end} length={length}; skipping.")
                continue

            chrom = parts[0]
            blocks[chrom].append({
                "start":    start,
                "end":      end,
                "length":   length,
                "label":    parts[4],
                "fwd":      fwd,
                "rev":      rev,
                "can":      can,
                "noncan":   noncan,
                "pathSize": path_size,
                "term":     parts[10] if len(parts) > 10 else "",
            })
    if malformed:
        suffix = " (first 3 shown above)" if malformed > 3 else ""
        _warn(f"Skipped {malformed} malformed terminal BED line(s) from '{path}'{suffix}.")
    return dict(blocks)


def _parse_bedgraph_fallback(path):
    """Robust line-by-line BEDgraph parser used if pandas parsing fails."""
    per_chrom = OrderedDict()
    malformed = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
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

            if end <= start or not np.isfinite(value):
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
        data[chrom] = (
            np.asarray(starts, dtype=np.int64),
            np.asarray(ends, dtype=np.int64),
            np.asarray(values, dtype=np.float64),
        )
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
            if line.startswith("#") or line.startswith("track"):
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

    data = OrderedDict()
    for chrom, grp in df.groupby("chrom", sort=False):
        data[chrom] = (grp["start"].values, grp["end"].values, grp["value"].values)
    return data


_TYPE_MAP = OrderedDict([
    ("t2t",                "T2T"),
    ("gapped_t2t",         "Gapped T2T"),
    ("incomplete",         "Incomplete"),
    ("gapped_incomplete",  "Gapped Incomplete"),
    ("misassembly",        "Misassembly"),
    ("gapped_misassembly", "Gapped Misassembly"),
    ("gapped_missassembly","Gapped Misassembly"),
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
        ("Misassembly",         []),
        ("Gapped Misassembly",  []),
        ("Discordant",          []),
        ("Gapped Discordant",   []),
        ("No telomeres",        []),
        ("Gapped No telomeres", []),
    ])
    header_idx = None
    type_col = None
    header_col = None
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
                    except ValueError:
                        _warn(f"{path}:{lineno}: report header is missing required 'header'/'type' columns; skipping section.")
                        header_col = None
                        type_col = None
                        continue
                    header_idx = 0
                continue
            if type_col is None or header_col is None or len(parts) <= max(header_col, type_col):
                continue
            chrom = parts[header_col]
            cat = _TYPE_MAP.get(parts[type_col].lower())
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


def _fmt_bp(bp):
    """Format base pairs for display."""
    if bp >= 1_000_000:
        return f"{bp/1e6:.2f} Mb"
    elif bp >= 1_000:
        return f"{bp/1e3:.1f} kb"
    return f"{bp:,} bp"


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
    ys[:n] = np.clip(values, 0, 1)
    xs[n] = ends[-1]
    ys[n] = 0
    return xs, ys


def _panel_label(ax, label, x=-0.05, y=1.15):
    """Add Nature-style panel label (bold lowercase letter, top-left)."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top", ha="left")


def _sanitize_filename(name):
    """Convert a chromosome name into a filesystem-safe stem."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return safe.strip("._") or "unnamed"


def _hide_x_axis(ax):
    """Hide x-axis ticks and labels for non-bottom tracks."""
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines["bottom"].set_visible(False)


def _style_fraction_axis(ax, label=None):
    """Style 0..1 quantitative tracks with minimal ticks."""
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.tick_params(axis="y", length=2.0, width=0.45, pad=2.2,
                   labelsize=5.5, colors="#555555")
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.45)
    ax.spines["left"].set_color("#bcbcbc")
    if label:
        ax.set_ylabel(label, fontsize=6.3, labelpad=4.5)
    else:
        ax.set_ylabel("")
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


def _apply_terminal_x_axis(ax, view_start, view_end, arm):
    """Apply end-relative x-axis with adaptive units (bp / kbp / Mbp)."""
    if view_end <= view_start:
        return

    span = view_end - view_start

    # Adaptive unit selection
    if span >= 2_000_000:
        fmt = lambda x: f"{x / 1e6:.1f}"
        unit = "Mbp"
    elif span >= 2_000:
        fmt = lambda x: _fmt_kbp(max(x, 0))
        unit = "kbp"
    else:
        fmt = lambda x: f"{max(x, 0):.0f}"
        unit = "bp"

    n_ticks = 4 if span >= 50_000 else 3

    if arm == "q":
        tick_pos = np.linspace(view_end, view_start, n_ticks)
    else:
        tick_pos = np.linspace(view_start, view_end, n_ticks)

    tick_labels = [fmt(pos) for pos in tick_pos]

    if arm in ("p", "q"):
        xlabel = f"Distance to end ({unit})"
    else:
        xlabel = f"Position along scaffold ({unit})"

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis="x", bottom=True, labelbottom=True,
                   length=2.3, width=0.45, pad=1.5)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(0.45)
    ax.spines["bottom"].set_color("#bcbcbc")
    ax.set_xlabel(xlabel, fontsize=6.3, labelpad=2.5)
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
    if block["label"] in {"p", "q"}:
        start_dist, end_dist = _terminal_distance_interval(
            int(block["start"]), int(block["end"]), int(chrom_size), block["label"])
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
        violin = ax.violinplot(
            [core_values], positions=[position], vert=vert,
            widths=0.36, showmeans=False, showmedians=False,
            showextrema=False,
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

    jitter = rng.uniform(-0.15, -0.08, size=len(values))
    filled_mask = ~outliers
    if np.any(filled_mask):
        xs = (np.full(np.sum(filled_mask), position) + jitter[filled_mask]
              if not vert else values[filled_mask])
        ys = (values[filled_mask]
              if vert else np.full(np.sum(filled_mask), position) + jitter[filled_mask])
        ax.scatter(xs, ys, s=7, color=color, alpha=0.42,
                   edgecolors="white", linewidths=0.22, rasterized=True, zorder=3)
    if np.any(outliers):
        xs = (np.full(np.sum(outliers), position) + jitter[outliers]
              if not vert else values[outliers])
        ys = (values[outliers]
              if vert else np.full(np.sum(outliers), position) + jitter[outliers])
        ax.scatter(xs, ys, s=7, facecolors="none", edgecolors=color,
                   linewidths=0.55, rasterized=True, zorder=4)

    box = ax.boxplot(
        [core_values],
        positions=[position],
        vert=vert,
        widths=0.065,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        manage_ticks=False,
    )
    for patch in box["boxes"]:
        patch.set_facecolor("white")
        patch.set_edgecolor(color)
        patch.set_linewidth(0.55)
    for key in ("whiskers", "caps", "medians"):
        for artist in box[key]:
            artist.set_color(color)
            artist.set_linewidth(0.55)


def _draw_balanced_legend(ax, handles, y_anchor=0.50):
    """Render a compact two-row legend within a dedicated legend axis."""
    ax.axis("off")
    if not handles:
        return

    ncol = max(1, int(np.ceil(len(handles) / 2.0)))
    leg = ax.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.5, y_anchor),
        ncol=ncol,
        fontsize=4.95,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        handlelength=0.88,
        handleheight=0.7,
        columnspacing=0.95,
        handletextpad=0.28,
        borderaxespad=0.0,
        labelspacing=0.50,
    )
    leg.get_frame().set_linewidth(0.4)


def _draw_summary_bar(ax, segments, title, x_label, fmt_value, x_formatter=None):
    """Draw a thin stacked summary bar for overview panel a."""
    total = sum(segment[1] for segment in segments)
    if title:
        ax.set_title(title, loc="center", fontsize=7.0, pad=0.0, y=0.97)

    if total <= 0:
        ax.text(0.5, 0.5, "No classification data",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=6.3, color="#999999")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        return

    left = 0.0
    for segment in segments:
        label, value, color = segment[:3]
        label_value = segment[3] if len(segment) > 3 else value
        if value <= 0:
            continue
        ax.barh(
            0.0, value, left=left, height=0.10,
            color=color, edgecolor="white", linewidth=0.6,
            rasterized=True, zorder=2,
        )
        if value / total >= 0.08:
            ax.text(
                left + value / 2.0, 0.0, fmt_value(label_value),
                ha="center", va="center",
                fontsize=4.9, color="#222222", zorder=3,
            )
        left += value

    ax.set_xlim(0, total)
    ax.set_ylim(-0.08, 0.14)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=2.0, width=0.45, pad=1.2)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#c7c7c7")
    ax.spines["bottom"].set_linewidth(0.45)
    ax.minorticks_off()
    if x_formatter is not None:
        ax.xaxis.set_major_formatter(x_formatter)
    ax.set_xlabel(x_label, fontsize=6.1, labelpad=0.6)


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
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, 1.8))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.01, 0.92, title, transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top", ha="left")
    ax.text(0.01, 0.72, textwrap.fill(message, width=95), transform=ax.transAxes,
            fontsize=6.5, va="top", ha="left", color="#444444")
    fig.subplots_adjust(left=0.03, right=0.98, top=0.95, bottom=0.08)
    return fig


def _save_figure_with_fallback(save_figure, build_figure, title, message_prefix):
    """Build and save a figure, falling back to a placeholder page on render errors."""
    fig = None
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


# ---------------------------------------------------------------------------
# Tier 1: Assembly overview
# ---------------------------------------------------------------------------

def plot_assembly_overview(classifications, blocks, chrom_sizes):
    """Overview page with classification summary and telomere block distributions."""
    FLAGGED_DISTANCE_BP = 1000
    FLAGGED_TOP_N = 10

    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, 4.72))

    # Two-row nested gridspec: row 0 = bars (2 cols), row 1 = c/d/e (3 cols)
    outer = fig.add_gridspec(2, 1, height_ratios=[0.58, 1.42], hspace=0.14)
    row0 = outer[0].subgridspec(1, 2, wspace=0.12)
    row1 = outer[1].subgridspec(1, 3, wspace=0.34, width_ratios=[1.0, 1.0, 0.92])

    abs_spec = row0[0, 0].subgridspec(2, 1, height_ratios=[1.0, 0.28], hspace=0.16)
    rel_spec = row0[0, 1].subgridspec(2, 1, height_ratios=[1.0, 0.28], hspace=0.16)

    ax_abs = fig.add_subplot(abs_spec[0])
    ax_abs_leg = fig.add_subplot(abs_spec[1])
    ax_rel = fig.add_subplot(rel_spec[0])
    ax_rel_leg = fig.add_subplot(rel_spec[1])
    ax_rain = fig.add_subplot(row1[0, 0])
    ax_scatter = fig.add_subplot(row1[0, 1])
    ax_flagged = fig.add_subplot(row1[0, 2])
    fig.subplots_adjust(left=0.068, right=0.958, top=0.86, bottom=0.11)

    # ---- Data setup ----
    cat_labels = list(classifications.keys())
    cat_counts = [len(v) for v in classifications.values()]
    cat_colors = [COLORS.get(c, "#aaaaaa") for c in cat_labels]
    total = sum(cat_counts)

    block_rows = []
    for chrom, blist in blocks.items():
        chrom_size = chrom_sizes.get(chrom, 0)
        if chrom_size <= 0 and blist:
            chrom_size = max(
                max((b.get("pathSize", 0) for b in blist), default=0),
                max(b["end"] for b in blist),
            )
        for block in blist:
            block_rows.append({
                "chrom": chrom,
                "label": block["label"],
                "length": block["length"],
                "distance": _block_end_distance(block, chrom_size),
            })

    # Flagged scaffolds: blocks far from scaffold ends (ITS / fragmentation suspects)
    flagged_rows = [r for r in block_rows if r["distance"] > FLAGGED_DISTANCE_BP]
    flagged_by_scaffold = defaultdict(lambda: {"length": 0, "count": 0, "arms": Counter()})
    for fr in flagged_rows:
        entry = flagged_by_scaffold[fr["chrom"]]
        entry["length"] += fr["length"]
        entry["count"] += 1
        entry["arms"][fr["label"]] += 1
    flagged_ranked = sorted(flagged_by_scaffold.items(),
                            key=lambda kv: kv[1]["length"], reverse=True)
    flagged_total_scaffolds = len(flagged_ranked)
    flagged_total_blocks = len(flagged_rows)
    flagged_top = flagged_ranked[:FLAGGED_TOP_N]

    p_len = [row["length"] for row in block_rows if row["label"] == "p"]
    q_len = [row["length"] for row in block_rows if row["label"] == "q"]
    b_len = [row["length"] for row in block_rows if row["label"] == "b"]
    all_len = p_len + q_len + b_len
    total_blocks = len(block_rows)
    total_paths = total if total > 0 else max(len(chrom_sizes), len(blocks))
    median_kbp_text = f"{_fmt_kbp(_median(all_len))} kbp" if all_len else "NA"

    # ---- Panel a/b: Classification summary ----
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

    _draw_summary_bar(
        ax_abs,
        absolute_segments,
        "",
        "Sequences (n)",
        fmt_value=lambda value: f"{int(round(value))}",
    )
    ax_abs.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=4))

    _draw_summary_bar(
        ax_rel,
        relative_segments,
        "",
        "Sequences (%)",
        fmt_value=lambda value: f"{int(round(value))}",
        x_formatter=ticker.FuncFormatter(lambda x, _: f"{x:.0f}%"),
    )
    ax_rel.xaxis.set_major_locator(ticker.MultipleLocator(25))
    abs_handles = [Patch(facecolor=color, label=lab) for lab, _, color in absolute_segments]
    rel_handles = [Patch(facecolor=color, label=lab) for lab, _, color in telomeric_segments]
    _draw_balanced_legend(ax_abs_leg, abs_handles, y_anchor=0.16)
    _draw_balanced_legend(ax_rel_leg, rel_handles, y_anchor=0.16)

    # ---- Panel c: Vertical raincloud plot (log10bp on y, arms on x) ----
    groups = [
        ("p-arm", p_len, COLORS["p"]),
        ("q-arm", q_len, COLORS["q"]),
        ("balanced", b_len, COLORS["b"]),
    ]
    groups = [(label, values, color) for label, values, color in groups if values]

    if groups:
        rng = np.random.default_rng(0)
        spacing = 0.36
        offsets = np.arange(len(groups), dtype=np.float64)
        offsets -= (len(groups) - 1) / 2.0
        positions = 0.75 + (spacing * offsets)
        for position, (label, values, color) in zip(positions, groups):
            log_vals = np.log10(np.asarray(values, dtype=np.float64) + 1.0)
            _draw_raincloud_group(ax_rain, log_vals, position, color, rng, vert=True)

        x_labels = [f"{label}\n(n={len(values)})" for label, values, _ in groups]
        ax_rain.set_xticks(positions)
        ax_rain.set_xticklabels(x_labels, fontsize=5.5)
        ax_rain.tick_params(axis="x", length=2, pad=3)
        ax_rain.tick_params(axis="y", length=2)
        ax_rain.set_xlim(positions.min() - 0.34, positions.max() + 0.30)
        ax_rain.set_box_aspect(1.0)
        ax_rain.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))

        median_log = np.log10(_median(all_len) + 1.0)
        ax_rain.axhline(
            median_log,
            color="black",
            linestyle=OVERVIEW_DASH_STYLE,
            linewidth=OVERVIEW_DASH_WIDTH,
        )
        median_trans = transforms.blended_transform_factory(ax_rain.transAxes, ax_rain.transData)
        ax_rain.annotate(
            f"median ({median_kbp_text})",
            xy=(0.03, median_log),
            xycoords=median_trans,
            xytext=(0, 2),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=5.3,
            color="#444444",
        )
        ax_rain.set_ylabel("Telomere length (log10bp)")
    else:
        label = "No labeled telomere blocks" if block_rows else "No telomere blocks"
        ax_rain.text(0.5, 0.5, label, ha="center", va="center",
                     transform=ax_rain.transAxes, fontsize=7, color="#999999")
        ax_rain.set_ylabel("Telomere length (log10bp)")
        ax_rain.set_xticks([])
        ax_rain.set_box_aspect(1.0)

    # ---- Panel d: Terminal offset (log10bp x) vs telomere length (kbp y) ----
    scatter_groups = [
        ("p-arm", "p", COLORS["p"]),
        ("q-arm", "q", COLORS["q"]),
        ("balanced", "b", COLORS["b"]),
    ]
    scatter_rows = [row for row in block_rows if row["label"] in {"p", "q", "b"}]
    if scatter_rows:
        max_log_x = 0.0
        for legend_label, arm, color in scatter_groups:
            arm_rows = [row for row in scatter_rows if row["label"] == arm]
            if not arm_rows:
                continue
            x = np.log10(np.asarray([row["distance"] for row in arm_rows], dtype=np.float64) + 1.0)
            y = np.asarray([row["length"] / 1e3 for row in arm_rows], dtype=np.float64)
            max_log_x = max(max_log_x, float(np.max(x)))
            ax_scatter.scatter(
                x, y,
                s=13,
                color=color,
                alpha=0.55,
                edgecolors="white",
                linewidths=0.28,
                rasterized=True,
                label=legend_label,
            )

        ax_scatter.set_xlabel("Distance to end (log10bp)")
        ax_scatter.set_ylabel("Telomere length (kbp)")
        ax_scatter.set_xticks([0, 1, 2, 3, 4])
        x_right = max(4.05, max_log_x + 0.10)
        ax_scatter.axvspan(3.0, x_right, facecolor="#ededed", linewidth=0, zorder=0)
        ax_scatter.set_xlim(-0.05, x_right)
        y_max = max(float(np.max(np.asarray([row["length"] / 1e3 for row in scatter_rows],
                                             dtype=np.float64))), 0.0)
        ax_scatter.set_ylim(-0.3, max(1.0, y_max * 1.08))
        ax_scatter.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
        ax_scatter.set_box_aspect(1.0)
        leg_d = ax_scatter.legend(
            loc="upper left",
            fontsize=5.6,
            frameon=True,
            facecolor="white",
            edgecolor="black",
            handletextpad=0.35,
            borderaxespad=0.2,
        )
        leg_d.get_frame().set_linewidth(0.4)

        ax_scatter.axvline(
            3.0,
            color="black",
            linewidth=OVERVIEW_DASH_WIDTH,
            linestyle=OVERVIEW_DASH_STYLE,
            zorder=1,
        )
    else:
        ax_scatter.text(0.5, 0.5, "No telomere blocks",
                        transform=ax_scatter.transAxes,
                        ha="center", va="center",
                        fontsize=7, color="#999999")
        ax_scatter.set_xlabel("Distance to end (log10bp)")
        ax_scatter.set_ylabel("Telomere length (kbp)")
        ax_scatter.set_box_aspect(1.0)

    # ---- Panel e: Flagged scaffold lollipop chart ----
    if not block_rows:
        ax_flagged.text(0.5, 0.5, "No telomere blocks",
                        transform=ax_flagged.transAxes,
                        ha="center", va="center",
                        fontsize=7, color="#999999")
        ax_flagged.set_xticks([])
        ax_flagged.set_yticks([])
        for spine in ax_flagged.spines.values():
            spine.set_visible(False)
    elif not flagged_top:
        ax_flagged.text(
            0.5, 0.5,
            f"No flagged blocks\n(distance > {FLAGGED_DISTANCE_BP // 1000} kbp)",
            transform=ax_flagged.transAxes,
            ha="center", va="center",
            fontsize=6, color="#999999", linespacing=1.4,
        )
        ax_flagged.set_xticks([])
        ax_flagged.set_yticks([])
        for spine in ax_flagged.spines.values():
            spine.set_visible(False)
    else:
        names = []
        log_lengths = []
        counts = []
        dot_colors = []
        line_counts = []
        for scaffold, info in flagged_top:
            wrapped = _wrap_scaffold_name(scaffold, width=18)
            names.append(wrapped)
            log_lengths.append(np.log10(info["length"] + 1.0))
            counts.append(info["count"])
            dominant_arm = info["arms"].most_common(1)[0][0]
            dot_colors.append(COLORS.get(dominant_arm, "#aaaaaa"))
            line_counts.append(wrapped.count("\n") + 1)

        y_pos = []
        cursor = 0.0
        for line_count in line_counts:
            y_pos.append(cursor)
            cursor += 1.0 + (0.62 * max(0, line_count - 1))
        y_pos = np.asarray(y_pos, dtype=np.float64)
        max_y = y_pos[-1] if len(y_pos) else 0.0

        for i, (yp, xv, col) in enumerate(zip(y_pos, log_lengths, dot_colors)):
            ax_flagged.hlines(yp, 0, xv, color=col, linewidth=0.8)
            ax_flagged.plot(xv, yp, "o", color=col, markersize=4.5,
                            markeredgecolor="white", markeredgewidth=0.3)
            ax_flagged.annotate(
                f"n={counts[i]}",
                xy=(xv, yp),
                xytext=(4, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=4.5,
                color="#555555",
                annotation_clip=False,
            )

        ax_flagged.set_ylim(max_y + 0.65, -0.55)
        ax_flagged.set_yticks(y_pos)
        ax_flagged.set_yticklabels(names, fontsize=5)
        for tick in ax_flagged.get_yticklabels():
            tick.set_horizontalalignment("right")
            tick.set_multialignment("right")
        x_max = max(log_lengths) if log_lengths else 1.0
        ax_flagged.set_xlim(0.0, x_max + 0.22)
        ax_flagged.set_xlabel("Flagged length (log10bp)", fontsize=6)
        ax_flagged.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax_flagged.spines["top"].set_visible(False)
        ax_flagged.spines["right"].set_visible(False)
        ax_flagged.tick_params(axis="y", pad=1.2)

    ax_flagged.set_box_aspect(2.35)
    ax_flagged.set_anchor("E")

    fig.canvas.draw()
    label_style = dict(fontsize=8, fontweight="bold", va="bottom", ha="left")
    title_style = dict(fontsize=8, ha="center", va="bottom")

    panel_meta = [
        (ax_abs, "a", "Scaffold classification (all)"),
        (ax_rel, "b", "Scaffold classification (with telomeres)"),
        (ax_rain, "c", "Length by arm"),
        (ax_scatter, "d", "Telomere positioning"),
        (ax_flagged, "e", "Flagged telomeres"),
    ]
    for ax, label, title in panel_meta:
        bbox = ax.get_position()
        fig.text(max(0.005, bbox.x0 - 0.032), bbox.y1 + 0.006, label, **label_style)
        fig.text((bbox.x0 + bbox.x1) / 2.0, bbox.y1 + 0.009, title, **title_style)

    summary_title = (
        rf"$\bf{{Assembly\ summary}}$: n={total_paths} | "
        rf"telomeres={total_blocks} | flagged scaffolds={flagged_total_scaffolds} | "
        rf"flagged blocks={flagged_total_blocks}"
    )
    fig.text(0.5, 0.975, summary_title,
             ha="center", va="top", fontsize=6.9, color="#222222")

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
    """
    PADDING_FACTOR = 2.0
    MIN_WINDOW = 1_000
    MAX_RATIO = 5.0

    def _arm_limit(arm_blocks, arm):
        if not arm_blocks:
            return None
        intervals = [
            _terminal_distance_interval(
                int(b["start"]), int(b["end"]), int(chrom_size), arm)
            for b in arm_blocks
        ]
        furthest = max(max(s, e) for s, e in intervals)
        return min(int(chrom_size),
                    max(int(round(furthest * PADDING_FACTOR)), MIN_WINDOW))

    p_blocks = [b for b in blocks_list if b["label"] in ("p", "b")]
    q_blocks = [b for b in blocks_list if b["label"] in ("q", "b")]
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
                       telomere_present=True):
    """Draw telomere blocks as a thin track on a backbone line."""
    if not telomere_present:
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.text(0.5, 0.5, "No telomere", transform=ax.transAxes,
                ha="center", va="center", fontsize=6.4, color="#999999")
        if label:
            ax.set_ylabel(label, fontsize=6.3, labelpad=3)
        else:
            ax.set_ylabel("")
        _hide_x_axis(ax)
        return True

    backbone_y = 0.5
    display_start, display_end = _project_terminal_interval(
        view_start, view_end, chrom_size, arm)
    ax.plot([display_start, display_end], [backbone_y, backbone_y],
            color="#d6d6d6", linewidth=0.9, solid_capstyle="round",
            zorder=1, rasterized=True)

    for b in blocks_list:
        if b["end"] <= view_start or b["start"] >= view_end:
            continue
        cs = max(b["start"], view_start)
        ce = min(b["end"], view_end)
        ds, de = _project_terminal_interval(cs, ce, chrom_size, arm)
        if de < ds:
            ds, de = de, ds
        # Enforce minimum visible width for tiny blocks in large views
        view_span = abs(display_end - display_start)
        min_width = view_span * 0.015
        block_width = de - ds
        expanded = block_width < min_width
        if expanded:
            center = (ds + de) / 2
            ds = center - min_width / 2
            de = center + min_width / 2
            disp_lo = min(display_start, display_end)
            disp_hi = max(display_start, display_end)
            if ds < disp_lo:
                ds, de = disp_lo, disp_lo + min_width
            if de > disp_hi:
                de, ds = disp_hi, disp_hi - min_width
        color = COLORS.get(b["label"], "#999999")
        rect = Rectangle(
            (ds, backbone_y - 0.07), de - ds, 0.14,
            facecolor=color,
            edgecolor=color if expanded else "none",
            linewidth=0.3 if expanded else 0,
            alpha=0.98, zorder=2,
        )
        rect.set_rasterized(True)
        ax.add_patch(rect)

    ax.set_ylim(0.28, 0.72)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    if label:
        ax.set_ylabel(label, fontsize=6.3, labelpad=3)
    else:
        ax.set_ylabel("")
    _hide_x_axis(ax)
    return True


def _draw_fraction_track(ax, track_data, view_start, view_end, chrom_size, arm, color, label=None):
    """Draw a quantitative 0..1 track as a thin step plot with light fill."""
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
                        alpha=0.16, linewidth=0, rasterized=True)
        ax.plot(xs, ys, drawstyle="steps-post", color=color,
                linewidth=0.6, rasterized=True)

    _style_fraction_axis(ax, label)
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
            ha="center", va="center", fontsize=7, color="#bbbbbb")


def plot_terminal_zoom(chrom, chrom_size, blocks_list,
                       density_data=None, canonical_data=None, strand_data=None):
    """Two-column terminal zoom: p-end (left) and q-end (right).

    Each column shows up to 4 tracks (blocks, density, canonical ratio, strand bias).
    If windows overlap on a short chromosome, a single merged panel is used.
    """
    has_density = density_data is not None and len(density_data[0]) > 0
    has_canonical = canonical_data is not None and len(canonical_data[0]) > 0
    has_strand = strand_data is not None and len(strand_data[0]) > 0
    p_has_telomere = any(b["label"] in ("p", "b") for b in blocks_list)
    q_has_telomere = any(b["label"] in ("q", "b") for b in blocks_list)

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

    n_tracks = len(track_specs)

    n_cols = 1 if merged else 2
    height_map = {
        "blocks": 0.18,
        "density": 0.30,
        "canonical": 0.30,
        "strand": 0.30,
    }
    height_ratios = [height_map[name] for name, _ in track_specs]

    fig_height = 3.18
    fig, axes = plt.subplots(
        n_tracks, n_cols,
        figsize=(FIG_WIDTH_DOUBLE, fig_height),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.32,
                     "wspace": 0.18},
        squeeze=False,
        sharey="row",
    )

    fig.subplots_adjust(left=0.095, right=0.97, top=0.785, bottom=0.16)

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
            else:
                visible = _draw_strand_track(
                    ax, track_data, view_start, view_end,
                    chrom_size, arm,
                    label="Strand\nbias" if show_label else None,
                )

            if visible:
                display_start, display_end = _project_terminal_interval(
                    view_start, view_end, chrom_size, arm)
                if arm == "q":
                    ax.set_xlim(display_end, display_start)
                else:
                    ax.set_xlim(display_start, display_end)
                visible_rows.append(row)

        if visible_rows:
            for row in visible_rows[:-1]:
                _hide_x_axis(axes[row][col_idx])
            display_start, display_end = _project_terminal_interval(
                view_start, view_end, chrom_size, arm)
            _apply_terminal_x_axis(
                axes[visible_rows[-1]][col_idx],
                display_start,
                display_end,
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
        axes[0][0].set_title("Full scaffold", fontsize=7, pad=3)
    else:
        if n_cols == 2:
            axes[0][0].set_title("p-arm", fontsize=7, pad=3)
            axes[0][1].set_title("q-arm", fontsize=7, pad=3)
            _panel_label(axes[0][0], "a", x=0.0, y=1.16)
            _panel_label(axes[0][1], "b", x=0.0, y=1.16)

    labels_present = [lab for lab in ("p", "q", "b") if any(b["label"] == lab for b in blocks_list)]
    block_handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=COLORS[lab],
               markersize=4.5, linestyle="none",
               label={"p": "p arm", "q": "q arm", "b": "balanced"}[lab])
        for lab in labels_present
    ]
    if block_handles and (merged or "b" in labels_present):
        fig.legend(
            handles=block_handles, loc="upper left",
            bbox_to_anchor=(0.18, 0.986), fontsize=5.6,
            frameon=False, ncol=len(block_handles),
            handletextpad=0.3, columnspacing=0.8,
        )

    fig.suptitle(f"{chrom}  ({_fmt_bp(chrom_size)})",
                 fontsize=8, fontweight="bold", y=0.985, x=0.09, ha="left")
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
                        help="Output path (default: <directory>/teloscope_report.pdf)")
    parser.add_argument("--png", action="store_true",
                        help="Save individual PNG files instead of a single PDF")
    parser.add_argument("--dpi", type=int, default=450,
                        help="DPI for raster output (default: 450)")
    parser.add_argument("--draft", action="store_true",
                        help="Draft mode: render at 150 DPI for fast iteration")
    args = parser.parse_args()

    if args.draft:
        args.dpi = 150

    if not os.path.isdir(args.directory):
        sys.exit(f"Error: '{args.directory}' is not a directory.")

    files = find_files(args.directory)
    if "terminal" not in files:
        sys.exit(f"Error: Missing teloscope output files in '{args.directory}'.\n"
                 f"Run teloscope first to generate output files.")
    if "report" not in files:
        _warn(f"No '*_report.tsv' file found in '{args.directory}'; the overview classification panel will show no data.")

    print(f"Found files: {', '.join(files.keys())}", file=sys.stderr)

    blocks = parse_terminal_bed(files["terminal"])
    density_data = parse_bedgraph(files["density"]) if "density" in files else None
    canonical_data = parse_bedgraph(files["canonical_ratio"]) if "canonical_ratio" in files else None
    strand_data  = parse_bedgraph(files["strand_ratio"]) if "strand_ratio" in files else None

    chrom_sizes = get_chrom_sizes(blocks, density_data, canonical_data, strand_data)

    classifications = parse_report(files["report"]) if "report" in files else OrderedDict()
    total_chroms = sum(len(v) for v in classifications.values())
    total_telo = sum(len(blist) for blist in blocks.values())
    cat_summary = ", ".join(f"{k}={len(v)}" for k, v in classifications.items()) or "none"
    print(f"Chromosomes: {total_chroms}  |  Telomere blocks: {total_telo}  |  "
          f"Categories: {cat_summary}",
          file=sys.stderr)

    # Only generate figures for chromosomes that have telomere blocks
    profile_chroms = [chrom for chrom in sorted(blocks.keys(),
                                                key=lambda c: chrom_sizes.get(c, 0), reverse=True)
                      if chrom_sizes.get(chrom, 0) > 0]
    skipped_no_size = sorted(set(blocks) - set(profile_chroms))
    if skipped_no_size:
        _warn(f"Skipping {len(skipped_no_size)} chromosome(s) with no usable size: {', '.join(skipped_no_size[:5])}"
              f"{' ...' if len(skipped_no_size) > 5 else ''}.")

    n_figures = len(profile_chroms) + 1  # +1 for overview
    fallback_pages = []

    # --- Write-and-close pattern: one figure in memory at a time ---
    if args.png:
        out_dir = args.output or args.directory
        os.makedirs(out_dir, exist_ok=True)

        # Overview figure
        overview_path = os.path.join(out_dir, "teloscope_overview.png")
        overview_ok, overview_err = _save_figure_with_fallback(
            lambda fig: fig.savefig(overview_path, dpi=args.dpi),
            lambda: plot_assembly_overview(classifications, blocks, chrom_sizes),
            "Assembly overview",
            "Failed to render the assembly overview. A placeholder image was written instead.",
        )
        if not overview_ok:
            fallback_pages.append(("overview", overview_err))
        overview_suffix = " [warning]" if not overview_ok else ""
        print(f"[1/{n_figures}] overview{overview_suffix}", file=sys.stderr)

        # Per-chromosome terminal zoom figures
        for i, chrom in enumerate(profile_chroms, start=2):
            csize = chrom_sizes.get(chrom, 0)
            if csize == 0:
                continue
            blist = blocks.get(chrom, [])
            den = density_data.get(chrom) if density_data else None
            can = canonical_data.get(chrom) if canonical_data else None
            strand = strand_data.get(chrom) if strand_data else None

            safe_name = _sanitize_filename(chrom)
            path = os.path.join(out_dir, f"teloscope_{safe_name}.png")
            ok, error_text = _save_figure_with_fallback(
                lambda fig, path=path: fig.savefig(path, dpi=args.dpi),
                lambda chrom=chrom, csize=csize, blist=blist, den=den, can=can, strand=strand:
                    plot_terminal_zoom(chrom, csize, blist, den, can, strand),
                chrom,
                f"Failed to render the terminal zoom for {chrom}. A placeholder image was written instead.",
            )
            if not ok:
                fallback_pages.append((chrom, error_text))
            suffix = " [warning]" if not ok else ""
            print(f"[{i}/{n_figures}] {chrom}{suffix}", file=sys.stderr)

        print(f"Figures saved to {out_dir}/", file=sys.stderr)

    else:
        out_path = args.output or os.path.join(args.directory, "teloscope_report.pdf")
        with PdfPages(out_path) as pdf:
            # Overview figure
            overview_ok, overview_err = _save_figure_with_fallback(
                lambda fig: pdf.savefig(fig),
                lambda: plot_assembly_overview(classifications, blocks, chrom_sizes),
                "Assembly overview",
                "Failed to render the assembly overview. A placeholder page was written instead.",
            )
            if not overview_ok:
                fallback_pages.append(("overview", overview_err))
            overview_suffix = " [warning]" if not overview_ok else ""
            print(f"[1/{n_figures}] overview{overview_suffix}", file=sys.stderr)

            # Per-chromosome terminal zoom figures
            for i, chrom in enumerate(profile_chroms, start=2):
                csize = chrom_sizes.get(chrom, 0)
                if csize == 0:
                    continue
                blist = blocks.get(chrom, [])
                den = density_data.get(chrom) if density_data else None
                can = canonical_data.get(chrom) if canonical_data else None
                strand = strand_data.get(chrom) if strand_data else None

                ok, error_text = _save_figure_with_fallback(
                    lambda fig: pdf.savefig(fig),
                    lambda chrom=chrom, csize=csize, blist=blist, den=den, can=can, strand=strand:
                        plot_terminal_zoom(chrom, csize, blist, den, can, strand),
                    chrom,
                    f"Failed to render the terminal zoom for {chrom}. A placeholder page was written instead.",
                )
                if not ok:
                    fallback_pages.append((chrom, error_text))
                suffix = " [warning]" if not ok else ""
                print(f"[{i}/{n_figures}] {chrom}{suffix}", file=sys.stderr)

        print(f"Report saved to {out_path}", file=sys.stderr)

    if fallback_pages:
        preview = ", ".join(f"{name} ({err})" for name, err in fallback_pages[:5])
        more = " ..." if len(fallback_pages) > 5 else ""
        _warn(f"Report generation completed with {len(fallback_pages)} placeholder figure(s): {preview}{more}")


if __name__ == "__main__":
    main()
