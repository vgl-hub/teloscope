#!/usr/bin/env python3
"""
teloscope_report.py — Publication-ready figures from Teloscope output.

Usage:
    python teloscope_report.py <output_directory> [-o report.pdf]
    python teloscope_report.py <output_directory> --png

Reads Teloscope output files from the given directory and generates:
  Page 1: Assembly overview (classification donut + telomere length distribution)
  Page 2+: Per-chromosome terminal zoom figures (blocks, density, strand ratio)

Requires: Python 3.6+, matplotlib, numpy, pandas
"""

import sys
import os
import glob
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Nature-style configuration
# ---------------------------------------------------------------------------

# Nature figure widths: 89 mm (single column), 183 mm (double column)
FIG_WIDTH_SINGLE = 3.50    # inches (89 mm)
FIG_WIDTH_DOUBLE = 7.20    # inches (183 mm)

COLORS = {
    # Classification palette (Nature-compatible, colorblind-friendly)
    "T2T":          "#00B945",
    "Incomplete":   "#FF9500",
    "No telomeres": "#9e9e9e",
    "Misassembly":  "#FF2C00",
    "Discordant":   "#845B97",
    # Arm / strand colors
    "p":            "#0C5DA5",
    "q":            "#FF2C00",
    "b":            "#845B97",
    "density":      "#0C5DA5",
    "fwd":          "#0C5DA5",
    "rev":          "#FF2C00",
    "no_data":      "#e0e0e0",
}

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
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.02,
        "savefig.transparent": False,
        # PDF — TrueType embedding (required by Nature)
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
    })

_apply_nature_style()

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
    }
    for key, pat in patterns.items():
        hits = glob.glob(os.path.join(directory, pat))
        if hits:
            files[key] = hits[0]
    return files


def parse_terminal_bed(path):
    """
    Parse *_terminal_telomeres.bed -> dict[chrom -> list of block dicts].
    Columns: chrom start end length label fwdCount revCount canCount nonCanCount pathSize terminality
    """
    blocks = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[0]
            blocks[chrom].append({
                "start":    int(parts[1]),
                "end":      int(parts[2]),
                "length":   int(parts[3]),
                "label":    parts[4],
                "fwd":      int(parts[5]),
                "rev":      int(parts[6]),
                "can":      int(parts[7]),
                "noncan":   int(parts[8]),
                "pathSize": int(parts[9]),
                "term":     parts[10] if len(parts) > 10 else "",
            })
    return dict(blocks)


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

    df = pd.read_csv(path, sep="\t", header=None, skiprows=skip,
                     names=["chrom", "start", "end", "value"],
                     dtype={"chrom": str, "start": np.int64,
                            "end": np.int64, "value": np.float64},
                     engine="c")

    data = OrderedDict()
    for chrom, grp in df.groupby("chrom", sort=False):
        data[chrom] = (grp["start"].values, grp["end"].values, grp["value"].values)
    return data


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------

def classify_chromosomes(blocks, all_chroms=None):
    """Classify chromosomes based on terminal telomere blocks.

    Replicates the C++ labelTerminalBlocks decision tree:
      1. Only scaffold-terminal blocks count for classification.
      2. Find the longest p block and longest q block (by canonical count).
      3. Check orientation validity: p must be closer to the left end,
         q must be closer to the right end. If violated -> Discordant.
      4. Both p and q present: p before q -> T2T, else -> Misassembly.
      5. Single label with duplicates -> Misassembly.
      6. Single block -> Incomplete.
      7. No scaffold-terminal blocks -> No telomeres.
    """
    cats = OrderedDict([
        ("T2T",          []),
        ("Incomplete",   []),
        ("Misassembly",  []),
        ("Discordant",   []),
        ("No telomeres", []),
    ])

    classified = set()
    for chrom, blist in blocks.items():
        classified.add(chrom)

        # Only consider scaffold-terminal blocks
        sblocks = [b for b in blist if b.get("term", "") == "scaffold"]
        if not sblocks:
            cats["No telomeres"].append(chrom)
            continue

        # Find longest p and q blocks by canonical count
        p_blocks = [b for b in sblocks if b["label"] == "p"]
        q_blocks = [b for b in sblocks if b["label"] == "q"]

        best_p = max(p_blocks, key=lambda b: b["can"]) if p_blocks else None
        best_q = max(q_blocks, key=lambda b: b["can"]) if q_blocks else None

        # Check orientation validity (discordant takes priority)
        has_invalid = False
        if best_p:
            left_dist = best_p["start"]
            right_dist = best_p["pathSize"] - (best_p["start"] + best_p["length"])
            if left_dist > right_dist:
                has_invalid = True
        if best_q:
            left_dist = best_q["start"]
            right_dist = best_q["pathSize"] - (best_q["start"] + best_q["length"])
            if left_dist < right_dist:
                has_invalid = True

        if has_invalid:
            cats["Discordant"].append(chrom)
        elif best_p and best_q:
            if best_p["start"] < best_q["start"]:
                cats["T2T"].append(chrom)
            else:
                cats["Misassembly"].append(chrom)
        elif best_p:
            if len(p_blocks) > 1:
                cats["Misassembly"].append(chrom)
            else:
                cats["Incomplete"].append(chrom)
        elif best_q:
            if len(q_blocks) > 1:
                cats["Misassembly"].append(chrom)
            else:
                cats["Incomplete"].append(chrom)
        else:
            cats["No telomeres"].append(chrom)

    if all_chroms:
        for ch in all_chroms:
            if ch not in classified:
                cats["No telomeres"].append(ch)

    return OrderedDict((k, v) for k, v in cats.items() if v)


def get_chrom_sizes(blocks, bedgraph_data=None):
    """Get chromosome sizes from BED pathSize or BEDgraph extents."""
    sizes = {}
    for chrom, blist in blocks.items():
        if blist:
            sizes[chrom] = blist[0]["pathSize"]
    if bedgraph_data:
        for chrom, (starts, ends, values) in bedgraph_data.items():
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
    ys[:n] = np.maximum(values, 0)
    xs[n] = ends[-1]
    ys[n] = 0
    return xs, ys


def _hex_to_rgb(hex_color):
    """Convert '#RRGGBB' to (r, g, b) floats in [0, 1]."""
    return (int(hex_color[1:3], 16) / 255.0,
            int(hex_color[3:5], 16) / 255.0,
            int(hex_color[5:7], 16) / 255.0)

# Pre-compute RGB tuples for strand blending
_FWD_RGB = np.array(_hex_to_rgb(COLORS["fwd"]))
_REV_RGB = np.array(_hex_to_rgb(COLORS["rev"]))
_NODATA_RGBA = np.array((*_hex_to_rgb(COLORS["no_data"]), 0.3))


def _blend_colors_array(values):
    """Vectorized strand color blending -> (N, 4) RGBA array.

    values: numpy array of strand ratios (1=fwd, 0=rev, <0=no data).
    """
    n = len(values)
    rgba = np.empty((n, 4))
    no_data = values < 0
    frac = values.copy()
    frac[no_data] = 0  # placeholder, overwritten below

    # Blend fwd/rev by fraction
    rgba[:, :3] = np.outer(frac, _FWD_RGB) + np.outer(1 - frac, _REV_RGB)
    rgba[:, 3] = 0.7

    # Override no-data entries
    rgba[no_data] = _NODATA_RGBA
    return rgba


def _panel_label(ax, label):
    """Add Nature-style panel label (bold lowercase letter, top-left)."""
    ax.text(-0.05, 1.15, label, transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------------------
# Tier 1: Assembly overview
# ---------------------------------------------------------------------------

def plot_assembly_overview(classifications, blocks, chrom_sizes):
    """Two-panel overview: (a) classification donut, (b) telomere length histogram."""
    fig, (ax_donut, ax_hist) = plt.subplots(
        1, 2, figsize=(FIG_WIDTH_DOUBLE, 2.6),
    )
    fig.subplots_adjust(wspace=0.45, left=0.06, right=0.97, top=0.85, bottom=0.18)

    # ---- Panel a: Classification donut ----
    _panel_label(ax_donut, "a")

    cat_labels = list(classifications.keys())
    cat_counts = [len(v) for v in classifications.values()]
    cat_colors = [COLORS.get(c, "#aaaaaa") for c in cat_labels]
    total = sum(cat_counts)

    if total > 0:
        wedges, _, autotexts = ax_donut.pie(
            cat_counts,
            labels=None,
            colors=cat_colors,
            autopct=lambda p: f"{int(round(p * total / 100))}",
            pctdistance=0.78,
            startangle=90,
            counterclock=False,
            wedgeprops=dict(width=0.45, edgecolor="white", linewidth=0.8),
            textprops=dict(fontsize=7, fontweight="bold", color="white"),
        )
        for w in wedges:
            w.set_rasterized(True)
        for at in autotexts:
            at.set_fontsize(7)
            at.set_fontweight("bold")

        legend_labels = [f"{lab} ({cnt})" for lab, cnt in zip(cat_labels, cat_counts)]
        ax_donut.legend(
            wedges, legend_labels, loc="center left",
            bbox_to_anchor=(-0.12, 0.5), fontsize=6, frameon=False,
            handlelength=0.8, handleheight=0.8,
        )
        ax_donut.set_title(f"Chromosome classification (n={total})", fontsize=8, pad=8)
    else:
        ax_donut.text(0.5, 0.5, "No data", ha="center", va="center",
                      transform=ax_donut.transAxes, fontsize=7, color="#999999")
        ax_donut.set_title("Chromosome classification", fontsize=8, pad=8)

    # ---- Panel b: Telomere length distribution ----
    _panel_label(ax_hist, "b")

    lengths, arm_labels = [], []
    for blist in blocks.values():
        for b in blist:
            lengths.append(b["length"])
            arm_labels.append(b["label"])

    if lengths:
        p_len = [l for l, a in zip(lengths, arm_labels) if a == "p"]
        q_len = [l for l, a in zip(lengths, arm_labels) if a == "q"]
        b_len = [l for l, a in zip(lengths, arm_labels) if a == "b"]
        all_len = p_len + q_len + b_len

        n_bins = min(max(8, len(all_len) // 3), 30)
        lo, hi = min(all_len), max(all_len)
        if lo == hi:
            bins = [lo - 50, lo + 50]
        else:
            step = (hi - lo) / n_bins
            bins = [lo + i * step for i in range(n_bins + 1)]

        data_stack, colors_stack, labels_stack = [], [], []
        for data, color, label in [
            (p_len, COLORS["p"], "p-arm (fwd)"),
            (q_len, COLORS["q"], "q-arm (rev)"),
            (b_len, COLORS["b"], "balanced"),
        ]:
            if data:
                data_stack.append(data)
                colors_stack.append(color)
                labels_stack.append(label)

        if data_stack:
            _, _, patches = ax_hist.hist(
                data_stack, bins=bins, stacked=True,
                color=colors_stack, label=labels_stack,
                edgecolor="white", linewidth=0.4, alpha=0.85,
                rasterized=True,
            )

        median_len = _median(all_len)
        ax_hist.text(
            0.97, 0.93, f"n={len(all_len)}  median={median_len:,.0f} bp",
            transform=ax_hist.transAxes, ha="right", va="top", fontsize=6,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc",
                      alpha=0.9, linewidth=0.4),
        )
        ax_hist.set_ylabel("Count")
        ax_hist.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax_hist.legend(fontsize=6, frameon=False, loc="upper left")

        if hi >= 5000:
            ax_hist.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, _: f"{x/1e3:.1f}"))
            ax_hist.set_xlabel("Telomere block length (kb)")
        else:
            ax_hist.set_xlabel("Telomere block length (bp)")
    else:
        ax_hist.text(0.5, 0.5, "No telomere blocks", ha="center", va="center",
                     transform=ax_hist.transAxes, fontsize=7, color="#999999")
        ax_hist.set_xlabel("Telomere block length")
        ax_hist.set_ylabel("Count")

    ax_hist.set_title("Telomere block lengths", fontsize=8, pad=8)
    return fig


# ---------------------------------------------------------------------------
# Tier 2: Terminal zoom figures
# ---------------------------------------------------------------------------

def compute_view_windows(blocks_list, chrom_size):
    """Compute (p_window, q_window) for terminal zoom panels.

    Each window is (start, end) or None if no blocks at that end.
    If both windows overlap on a short chromosome, returns a single merged window.
    """
    p_blocks = [b for b in blocks_list if b["label"] in ("p", "b")]
    q_blocks = [b for b in blocks_list if b["label"] in ("q", "b")]

    p_window = None
    q_window = None

    if p_blocks:
        p_max_end = max(b["end"] for b in p_blocks)
        extent = p_max_end
        flanking = max(extent * 2, 20_000)
        flanking = min(flanking, 500_000)
        p_window = (0, min(p_max_end + flanking, chrom_size))

    if q_blocks:
        q_min_start = min(b["start"] for b in q_blocks)
        extent = chrom_size - q_min_start
        flanking = max(extent * 2, 20_000)
        flanking = min(flanking, 500_000)
        q_window = (max(q_min_start - flanking, 0), chrom_size)

    # If windows overlap, merge into a single full-width window
    if p_window and q_window and p_window[1] >= q_window[0]:
        merged = (0, chrom_size)
        return merged, None

    return p_window, q_window


def clip_bedgraph(bg_tuple, view_start, view_end):
    """Return (starts, ends, values) arrays clipped to [view_start, view_end)."""
    starts, ends, values = bg_tuple
    mask = (ends > view_start) & (starts < view_end)
    cs = np.maximum(starts[mask], view_start)
    ce = np.minimum(ends[mask], view_end)
    return cs, ce, values[mask]


def _draw_blocks_track(ax, blocks_list, view_start, view_end):
    """Draw telomere blocks as rectangles on a backbone line."""
    backbone_y = 0.5
    ax.plot([view_start, view_end], [backbone_y, backbone_y],
            color="#d0d0d0", linewidth=2.5, solid_capstyle="round",
            zorder=1, rasterized=True)

    for b in blocks_list:
        if b["end"] <= view_start or b["start"] >= view_end:
            continue
        cs = max(b["start"], view_start)
        ce = min(b["end"], view_end)
        color = COLORS.get(b["label"], "#999999")
        rect = Rectangle(
            (cs, backbone_y - 0.3), ce - cs, 0.6,
            facecolor=color, edgecolor="none", alpha=0.85, zorder=2,
        )
        rect.set_rasterized(True)
        ax.add_patch(rect)

    ax.set_ylim(-0.1, 1.1)
    ax.set_ylabel("Blocks", fontsize=7)
    ax.set_yticks([])
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(bottom=False)


def _draw_density_track(ax, density_data, view_start, view_end):
    """Draw repeat density as a filled step plot."""
    starts, ends, values = clip_bedgraph(density_data, view_start, view_end)
    if len(starts) == 0:
        ax.set_visible(False)
        return
    xs, ys = _bedgraph_to_step(starts, ends, values)
    ax.fill_between(xs, ys, step="post", color=COLORS["density"],
                    alpha=0.35, linewidth=0, rasterized=True)
    ax.plot(xs, ys, drawstyle="steps-post", color=COLORS["density"],
            linewidth=0.6, rasterized=True)
    ax.set_ylabel("Repeat\ndensity", fontsize=7)
    ax.set_ylim(0, 1.05)
    ax.tick_params(bottom=False)
    ax.spines["bottom"].set_visible(False)


def _draw_strand_track(ax, strand_data, view_start, view_end):
    """Draw strand ratio as colored bars (vectorized)."""
    starts, ends, values = clip_bedgraph(strand_data, view_start, view_end)
    if len(starts) == 0:
        ax.set_visible(False)
        return

    lefts = starts
    widths = ends - starts
    rgba = _blend_colors_array(values)

    bars = ax.barh(
        np.full(len(lefts), 0.5), widths, left=lefts, height=1.0,
        color=rgba, edgecolor="none", rasterized=True,
    )

    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("Strand\nratio", fontsize=7)


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
                       density_data=None, strand_data=None):
    """Two-column terminal zoom: p-end (left) and q-end (right).

    Each column shows up to 3 tracks (blocks, density, strand ratio).
    If windows overlap on a short chromosome, a single merged panel is used.
    """
    has_density = density_data is not None and len(density_data[0]) > 0
    has_strand = strand_data is not None and len(strand_data[0]) > 0

    p_window, q_window = compute_view_windows(blocks_list, chrom_size)
    merged = (p_window is not None and q_window is None and
              any(b["label"] in ("q", "b") for b in blocks_list) and
              any(b["label"] in ("p", "b") for b in blocks_list))

    n_tracks = 1
    if has_density:
        n_tracks += 1
    if has_strand:
        n_tracks += 1

    n_cols = 1 if merged else 2
    height_ratios = [0.35]
    if has_density:
        height_ratios.append(1.0)
    if has_strand:
        height_ratios.append(0.5)

    fig_height = 0.8 + 0.9 * (n_tracks - 1)
    fig, axes = plt.subplots(
        n_tracks, n_cols,
        figsize=(FIG_WIDTH_DOUBLE, fig_height),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.08,
                     "wspace": 0.25},
        squeeze=False,
    )

    fig.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.18)

    # Determine which columns to draw
    if merged:
        # Single merged panel
        columns = [(0, p_window)]
    else:
        columns = []
        if p_window:
            columns.append((0, p_window))
        if q_window:
            columns.append((1 if n_cols == 2 else 0, q_window))

    # Track which column indices are actually drawn
    drawn_cols = set()
    for col_idx, window in columns:
        drawn_cols.add(col_idx)
        view_start, view_end = window
        row = 0

        # Blocks track
        _draw_blocks_track(axes[row][col_idx], blocks_list, view_start, view_end)
        axes[row][col_idx].set_xlim(view_start, view_end)
        row += 1

        # Density track
        if has_density:
            _draw_density_track(axes[row][col_idx], density_data, view_start, view_end)
            axes[row][col_idx].set_xlim(view_start, view_end)
            row += 1

        # Strand ratio track
        if has_strand:
            _draw_strand_track(axes[row][col_idx], strand_data, view_start, view_end)
            axes[row][col_idx].set_xlim(view_start, view_end)
            row += 1

        # Format x-axis on the bottom track
        view_span = view_end - view_start
        _format_bp_axis(axes[n_tracks - 1][col_idx], view_span)

    # Share Y axes between columns in the same row
    if n_cols == 2 and len(drawn_cols) == 2:
        for row in range(n_tracks):
            axes[row][0].get_shared_y_axes().join(axes[row][0], axes[row][1])

    # Hide columns with no telomere
    if n_cols == 2 and not merged:
        for col_idx in range(n_cols):
            if col_idx not in drawn_cols:
                for row in range(n_tracks):
                    _hide_panel(axes[row][col_idx])

    # Column titles
    if merged:
        axes[0][0].set_title("Full chromosome", fontsize=7, pad=4)
    else:
        if n_cols == 2:
            p_title = "p-end (5')" if p_window or 0 not in drawn_cols else "p-end (5')"
            q_title = "q-end (3')" if q_window or 1 not in drawn_cols else "q-end (3')"
            axes[0][0].set_title(p_title, fontsize=7, pad=4)
            axes[0][1].set_title(q_title, fontsize=7, pad=4)

    # Block legend on top-left panel
    labels_present = set(b["label"] for b in blocks_list)
    legend_elements = []
    for lab, name in [("p", "p-arm (fwd)"), ("q", "q-arm (rev)"), ("b", "balanced")]:
        if lab in labels_present:
            legend_elements.append(
                Line2D([0], [0], marker="s", color="w", markerfacecolor=COLORS[lab],
                       markersize=5, linestyle="none", label=name))
    if legend_elements:
        first_drawn_col = min(drawn_cols)
        axes[0][first_drawn_col].legend(
            handles=legend_elements, loc="upper right", fontsize=6,
            frameon=False, ncol=len(legend_elements), handletextpad=0.2,
        )

    # Strand ratio legend
    if has_strand:
        strand_row = n_tracks - 1
        first_drawn_col = min(drawn_cols)
        axes[strand_row][first_drawn_col].legend(
            handles=[
                Patch(facecolor=COLORS["fwd"], label="Forward"),
                Patch(facecolor="#808080", label="Mixed"),
                Patch(facecolor=COLORS["rev"], label="Reverse"),
                Patch(facecolor=COLORS["no_data"], alpha=0.6, label="No data"),
            ],
            loc="upper right", fontsize=5.5, frameon=False, ncol=4,
            handlelength=0.8, handleheight=0.7, handletextpad=0.2,
        )

    fig.suptitle(f"{chrom}  ({_fmt_bp(chrom_size)})",
                 fontsize=8, fontweight="bold", y=0.95)
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
        sys.exit(f"Error: No *_terminal_telomeres.bed found in '{args.directory}'.\n"
                 f"Run teloscope first to generate output files.")

    print(f"Found files: {', '.join(files.keys())}", file=sys.stderr)

    blocks = parse_terminal_bed(files["terminal"])
    density_data = parse_bedgraph(files["density"]) if "density" in files else None
    strand_data  = parse_bedgraph(files["strand_ratio"]) if "strand_ratio" in files else None

    all_chroms = list(density_data.keys()) if density_data else None
    chrom_sizes = get_chrom_sizes(blocks, density_data)

    classifications = classify_chromosomes(blocks, all_chroms)
    total_chroms = sum(len(v) for v in classifications.values())
    total_telo = sum(len(blist) for blist in blocks.values())
    print(f"Chromosomes: {total_chroms}  |  Telomere blocks: {total_telo}  |  "
          f"Categories: {', '.join(f'{k}={len(v)}' for k, v in classifications.items())}",
          file=sys.stderr)

    # Only generate figures for chromosomes that have telomere blocks
    profile_chroms = sorted(blocks.keys(),
                            key=lambda c: chrom_sizes.get(c, 0), reverse=True)

    n_figures = len(profile_chroms) + 1  # +1 for overview

    # --- Write-and-close pattern: one figure in memory at a time ---
    if args.png:
        out_dir = args.output or args.directory
        os.makedirs(out_dir, exist_ok=True)

        # Overview figure
        fig_overview = plot_assembly_overview(classifications, blocks, chrom_sizes)
        overview_path = os.path.join(out_dir, "teloscope_overview.png")
        fig_overview.savefig(overview_path, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig_overview)
        print(f"[1/{n_figures}] overview", file=sys.stderr)

        # Per-chromosome terminal zoom figures
        for i, chrom in enumerate(profile_chroms, start=2):
            csize = chrom_sizes.get(chrom, 0)
            if csize == 0:
                continue
            blist = blocks.get(chrom, [])
            den = density_data.get(chrom) if density_data else None
            strand = strand_data.get(chrom) if strand_data else None

            fig = plot_terminal_zoom(chrom, csize, blist, den, strand)
            safe_name = chrom.replace("/", "_").replace("\\", "_")
            path = os.path.join(out_dir, f"teloscope_{safe_name}.png")
            fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"[{i}/{n_figures}] {chrom}", file=sys.stderr)

        print(f"Figures saved to {out_dir}/", file=sys.stderr)

    else:
        out_path = args.output or os.path.join(args.directory, "teloscope_report.pdf")
        with PdfPages(out_path) as pdf:
            # Overview figure
            fig_overview = plot_assembly_overview(classifications, blocks, chrom_sizes)
            pdf.savefig(fig_overview, bbox_inches="tight")
            plt.close(fig_overview)
            print(f"[1/{n_figures}] overview", file=sys.stderr)

            # Per-chromosome terminal zoom figures
            for i, chrom in enumerate(profile_chroms, start=2):
                csize = chrom_sizes.get(chrom, 0)
                if csize == 0:
                    continue
                blist = blocks.get(chrom, [])
                den = density_data.get(chrom) if density_data else None
                strand = strand_data.get(chrom) if strand_data else None

                fig = plot_terminal_zoom(chrom, csize, blist, den, strand)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
                print(f"[{i}/{n_figures}] {chrom}", file=sys.stderr)

        print(f"Report saved to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
