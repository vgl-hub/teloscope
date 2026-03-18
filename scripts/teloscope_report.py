#!/usr/bin/env python3
"""
teloscope_report.py — Publication-ready figures from Teloscope output.

Usage:
    python teloscope_report.py <output_directory> [-o report.pdf]
    python teloscope_report.py <output_directory> --png

Reads Teloscope output files from the given directory and generates:
  Page 1: Assembly overview (classification summary + telomere length distribution)
  Page 2+: Per-chromosome terminal zoom figures (blocks, density, canonical ratio, strand ratio)

Requires: Python 3.6+, matplotlib, numpy, pandas
"""

import sys
import os
import glob
import re
import argparse
import textwrap
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
    "canonical":    "#007B5F",
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


_TYPE_MAP = {
    "t2t":                "T2T",
    "gapped_t2t":         "T2T",
    "misassembly":        "Misassembly",
    "gapped_misassembly": "Misassembly",
    "incomplete":         "Incomplete",
    "gapped_incomplete":  "Incomplete",
    "none":               "No telomeres",
    "gapped_none":        "No telomeres",
    "discordant":         "Discordant",
    "gapped_discordant":  "Discordant",
}


def parse_report(path):
    """Parse *_report.tsv -> OrderedDict[category -> list of chrom names].

    Reads the type column from the Path Summary table written by the C++ tool.
    Skips the Assembly Summary sections (lines starting with +++).
    """
    cats = OrderedDict([
        ("T2T",          []),
        ("Incomplete",   []),
        ("Misassembly",  []),
        ("Discordant",   []),
        ("No telomeres", []),
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
        _warn(f"No usable classification rows were parsed from '{path}'; the overview donut will show no data.")

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
    no_data = (~np.isfinite(values)) | (values < 0)
    frac = np.clip(values.copy(), 0, 1)
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
    ax.tick_params(axis="y", length=2.0, width=0.45, pad=1.5,
                   labelsize=5.5, colors="#555555")
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(0.45)
    ax.spines["left"].set_color("#bcbcbc")
    if label:
        ax.set_ylabel(label, fontsize=6.3, labelpad=3)
    else:
        ax.set_ylabel("")
    _hide_x_axis(ax)


def _apply_terminal_x_axis(ax, view_start, view_end, arm):
    """Apply a minimal end-relative x-axis for p/q panels."""
    if view_end <= view_start:
        return

    tick_pos = np.linspace(view_start, view_end, 3)
    if arm == "q":
        tick_labels = [_fmt_kbp(view_end - pos) for pos in tick_pos]
        xlabel = "Distance to end (kbp)"
    elif arm == "p":
        tick_labels = [_fmt_kbp(pos - view_start) for pos in tick_pos]
        xlabel = "Distance to end (kbp)"
    else:
        tick_labels = [_fmt_kbp(pos) for pos in tick_pos]
        xlabel = "Position along scaffold (kbp)"

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
    """Overview page with classification summary and telomere length histogram."""
    fig = plt.figure(figsize=(FIG_WIDTH_DOUBLE, 3.2))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[0.95, 1.85],
        width_ratios=[1.55, 0.95],
        hspace=0.42,
        wspace=0.22,
    )
    ax_class = fig.add_subplot(gs[0, 0])
    ax_stats = fig.add_subplot(gs[0, 1])
    ax_hist = fig.add_subplot(gs[1, :])
    fig.subplots_adjust(left=0.07, right=0.97, top=0.9, bottom=0.14)

    _panel_label(ax_class, "a")
    _panel_label(ax_hist, "b")

    cat_labels = list(classifications.keys())
    cat_counts = [len(v) for v in classifications.values()]
    cat_colors = [COLORS.get(c, "#aaaaaa") for c in cat_labels]
    total = sum(cat_counts)

    # ---- Panel a: Classification composition ----
    ax_class.set_title("Assembly classification", loc="left", fontsize=8, pad=4)
    if total > 0:
        left = 0
        for lab, cnt, color in zip(cat_labels, cat_counts, cat_colors):
            bar = ax_class.barh(
                0, cnt, left=left, height=0.42,
                color=color, edgecolor="white", linewidth=0.8,
                rasterized=True,
            )[0]
            if cnt / total >= 0.1:
                ax_class.text(
                    left + cnt / 2, 0, f"{cnt}",
                    ha="center", va="center", fontsize=6,
                    color="white", fontweight="bold",
                )
            left += cnt

        ax_class.set_xlim(0, total)
        ax_class.set_ylim(-0.7, 0.7)
        ax_class.set_xlabel("Chromosomes", fontsize=6.3, labelpad=2)
        ax_class.set_yticks([])
        ax_class.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=4))
        ax_class.grid(axis="x", color="#ececec", linewidth=0.45)
        ax_class.spines["left"].set_visible(False)
        ax_class.spines["bottom"].set_visible(True)
        ax_class.spines["bottom"].set_color("#bcbcbc")
        ax_class.spines["bottom"].set_linewidth(0.45)
    else:
        ax_class.text(
            0.5, 0.5, "No classification data",
            ha="center", va="center", transform=ax_class.transAxes,
            fontsize=7, color="#999999",
        )
        ax_class.set_xticks([])
        ax_class.set_yticks([])
        ax_class.spines["left"].set_visible(False)
        ax_class.spines["bottom"].set_visible(False)

    # ---- Stats panel ----
    ax_stats.axis("off")

    lengths, arm_labels = [], []
    for blist in blocks.values():
        for b in blist:
            lengths.append(b["length"])
            arm_labels.append(b["label"])

    p_len = [l for l, a in zip(lengths, arm_labels) if a == "p"]
    q_len = [l for l, a in zip(lengths, arm_labels) if a == "q"]
    b_len = [l for l, a in zip(lengths, arm_labels) if a == "b"]
    all_len = p_len + q_len + b_len
    total_blocks = sum(len(blist) for blist in blocks.values())
    total_paths = total if total > 0 else max(len(chrom_sizes), len(blocks))
    median_text = _fmt_bp(int(round(_median(all_len)))) if all_len else "NA"

    summary_rows = [
        ("Chromosomes", f"{total_paths}"),
        ("Terminal blocks", f"{total_blocks}"),
        ("Median block", median_text),
    ]
    y = 0.96
    for label, value in summary_rows:
        ax_stats.text(0.0, y, label, transform=ax_stats.transAxes,
                      ha="left", va="top", fontsize=6, color="#666666")
        ax_stats.text(1.0, y, value, transform=ax_stats.transAxes,
                      ha="right", va="top", fontsize=7.2, fontweight="bold",
                      color="#222222")
        y -= 0.2

    y -= 0.06
    if total > 0:
        for lab, cnt, color in zip(cat_labels, cat_counts, cat_colors):
            pct = 100.0 * cnt / total
            ax_stats.add_patch(Rectangle(
                (0.0, y - 0.028), 0.04, 0.04,
                transform=ax_stats.transAxes,
                facecolor=color, edgecolor="none",
            ))
            ax_stats.text(0.07, y, lab, transform=ax_stats.transAxes,
                          ha="left", va="center", fontsize=6.1, color="#333333")
            ax_stats.text(1.0, y, f"{cnt} ({pct:.0f}%)", transform=ax_stats.transAxes,
                          ha="right", va="center", fontsize=6.1, color="#333333")
            y -= 0.14
    else:
        ax_stats.text(0.0, y, "Classification report unavailable",
                      transform=ax_stats.transAxes, ha="left", va="center",
                      fontsize=6.1, color="#999999")

    # ---- Panel b: Telomere length distribution ----
    ax_hist.set_title("Telomere block length distribution", loc="left", fontsize=8, pad=4)
    if all_len:
        n_bins = min(max(8, len(all_len) // 3), 30)
        lo, hi = min(all_len), max(all_len)
        if lo == hi:
            bins = [lo - 50, lo + 50]
        else:
            step = (hi - lo) / n_bins
            bins = [lo + i * step for i in range(n_bins + 1)]

        data_stack, colors_stack, labels_stack = [], [], []
        for data, color, label in [
            (p_len, COLORS["p"], "p arm"),
            (q_len, COLORS["q"], "q arm"),
            (b_len, COLORS["b"], "balanced"),
        ]:
            if data:
                data_stack.append(data)
                colors_stack.append(color)
                labels_stack.append(label)

        if data_stack:
            ax_hist.hist(
                data_stack, bins=bins, stacked=True,
                color=colors_stack, label=labels_stack,
                edgecolor="white", linewidth=0.4, alpha=0.95,
                rasterized=True,
            )

        median_len = _median(all_len)
        ax_hist.axvline(median_len, color="#555555", linestyle="--", linewidth=0.7)
        ax_hist.text(
            0.99, 0.95, f"n={len(all_len)}   median={_fmt_bp(int(round(median_len)))}",
            transform=ax_hist.transAxes, ha="right", va="top", fontsize=6,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc",
                      alpha=0.9, linewidth=0.4),
        )
        ax_hist.set_ylabel("Count")
        ax_hist.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax_hist.grid(axis="y", color="#ececec", linewidth=0.45)
        ax_hist.legend(fontsize=5.8, frameon=False, loc="upper left",
                       ncol=min(3, len(labels_stack)), handlelength=1.0,
                       columnspacing=0.8, handletextpad=0.4)

        if hi >= 5000:
            ax_hist.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, _: _fmt_kbp(max(x, 0))))
            ax_hist.set_xlabel("Telomere block length (kbp)")
        else:
            ax_hist.set_xlabel("Telomere block length (bp)")
    else:
        label = "No labeled telomere blocks" if lengths else "No telomere blocks"
        ax_hist.text(0.5, 0.5, label, ha="center", va="center",
                     transform=ax_hist.transAxes, fontsize=7, color="#999999")
        ax_hist.set_xlabel("Telomere block length")
        ax_hist.set_ylabel("Count")

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


def _draw_blocks_track(ax, blocks_list, view_start, view_end, label=None):
    """Draw telomere blocks as a thin track on a backbone line."""
    backbone_y = 0.5
    ax.plot([view_start, view_end], [backbone_y, backbone_y],
            color="#d6d6d6", linewidth=1.4, solid_capstyle="round",
            zorder=1, rasterized=True)

    for b in blocks_list:
        if b["end"] <= view_start or b["start"] >= view_end:
            continue
        cs = max(b["start"], view_start)
        ce = min(b["end"], view_end)
        color = COLORS.get(b["label"], "#999999")
        rect = Rectangle(
            (cs, backbone_y - 0.18), ce - cs, 0.36,
            facecolor=color, edgecolor="none", alpha=0.96, zorder=2,
        )
        rect.set_rasterized(True)
        ax.add_patch(rect)

    ax.set_ylim(0.12, 0.88)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    if label:
        ax.set_ylabel(label, fontsize=6.3, labelpad=3)
    else:
        ax.set_ylabel("")
    _hide_x_axis(ax)
    return True


def _draw_fraction_track(ax, track_data, view_start, view_end, color, label=None):
    """Draw a quantitative 0..1 track as a thin step plot with light fill."""
    starts, ends, values = clip_bedgraph(track_data, view_start, view_end)
    if len(starts) == 0:
        ax.set_visible(False)
        return False

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


def _draw_strand_track(ax, strand_data, view_start, view_end, label=None):
    """Draw strand ratio as a thin color band."""
    starts, ends, values = clip_bedgraph(strand_data, view_start, view_end)
    if len(starts) == 0:
        ax.set_visible(False)
        return False

    lefts = starts
    widths = ends - starts
    rgba = _blend_colors_array(values)

    ax.barh(
        np.full(len(lefts), 0.5), widths, left=lefts, height=0.64,
        color=rgba, edgecolor="none", rasterized=True,
    )

    ax.set_ylim(0.1, 0.9)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if label:
        ax.set_ylabel(label, fontsize=6.3, labelpad=3)
    else:
        ax.set_ylabel("")
    _hide_x_axis(ax)
    return True


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

    Each column shows up to 4 tracks (blocks, density, canonical ratio, strand ratio).
    If windows overlap on a short chromosome, a single merged panel is used.
    """
    has_density = density_data is not None and len(density_data[0]) > 0
    has_canonical = canonical_data is not None and len(canonical_data[0]) > 0
    has_strand = strand_data is not None and len(strand_data[0]) > 0

    p_window, q_window = compute_view_windows(blocks_list, chrom_size)
    fallback_full_chrom = p_window is None and q_window is None
    if fallback_full_chrom:
        p_window = (0, chrom_size)
    merged = fallback_full_chrom or (
        p_window is not None and q_window is None and
        any(b["label"] in ("q", "b") for b in blocks_list) and
        any(b["label"] in ("p", "b") for b in blocks_list)
    )

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
        "blocks": 0.28,
        "density": 0.42,
        "canonical": 0.42,
        "strand": 0.18,
    }
    height_ratios = [height_map[name] for name, _ in track_specs]

    fig_height = 1.15 + 1.0 * sum(height_ratios)
    fig, axes = plt.subplots(
        n_tracks, n_cols,
        figsize=(FIG_WIDTH_DOUBLE, fig_height),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.08,
                     "wspace": 0.18},
        squeeze=False,
        sharey="row",
    )

    fig.subplots_adjust(left=0.09, right=0.97, top=0.82, bottom=0.16)

    # Determine which columns to draw
    if merged:
        # Single merged panel
        columns = [(0, p_window, "full")]
    else:
        columns = []
        if p_window:
            columns.append((0, p_window, "p"))
        if q_window:
            columns.append((1 if n_cols == 2 else 0, q_window, "q"))

    label_col = columns[0][0] if columns else 0

    # Track which column indices are actually drawn
    drawn_cols = set()
    strand_has_no_data = False
    for _, window, _ in columns:
        if has_strand:
            _, _, values = clip_bedgraph(strand_data, window[0], window[1])
            if np.any(values < 0):
                strand_has_no_data = True

    for col_idx, window, arm in columns:
        drawn_cols.add(col_idx)
        view_start, view_end = window
        visible_rows = []

        for row, (track_name, track_data) in enumerate(track_specs):
            ax = axes[row][col_idx]
            show_label = col_idx == label_col

            if track_name == "blocks":
                visible = _draw_blocks_track(
                    ax, blocks_list, view_start, view_end,
                    label="Blocks" if show_label else None,
                )
            elif track_name == "density":
                visible = _draw_fraction_track(
                    ax, track_data, view_start, view_end,
                    COLORS["density"],
                    label="Repeat\ndensity" if show_label else None,
                )
            elif track_name == "canonical":
                visible = _draw_fraction_track(
                    ax, track_data, view_start, view_end,
                    COLORS["canonical"],
                    label="Canonical\nratio" if show_label else None,
                )
            else:
                visible = _draw_strand_track(
                    ax, track_data, view_start, view_end,
                    label="Strand\nratio" if show_label else None,
                )

            if visible:
                ax.set_xlim(view_start, view_end)
                visible_rows.append(row)

        if visible_rows:
            for row in visible_rows[:-1]:
                _hide_x_axis(axes[row][col_idx])
            _apply_terminal_x_axis(
                axes[visible_rows[-1]][col_idx],
                view_start,
                view_end,
                arm,
            )

    # Hide columns with no telomere
    if n_cols == 2 and not merged:
        for col_idx in range(n_cols):
            if col_idx not in drawn_cols:
                for row in range(n_tracks):
                    _hide_panel(axes[row][col_idx])

    # Column titles
    if merged:
        axes[0][0].set_title("Full scaffold", fontsize=7, pad=3)
    else:
        if n_cols == 2:
            axes[0][0].set_title("p-end (5')", fontsize=7, pad=3)
            axes[0][1].set_title("q-end (3')", fontsize=7, pad=3)

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
            bbox_to_anchor=(0.09, 0.925), fontsize=5.6,
            frameon=False, ncol=len(block_handles),
            handletextpad=0.3, columnspacing=0.8,
        )

    if has_strand:
        strand_handles = [
            Patch(facecolor=COLORS["rev"], label="Rev"),
            Patch(facecolor="#808080", label="Mix"),
            Patch(facecolor=COLORS["fwd"], label="Fwd"),
        ]
        if strand_has_no_data:
            strand_handles.append(Patch(facecolor=COLORS["no_data"], alpha=0.6, label="No data"))
        fig.legend(
            handles=strand_handles, loc="upper right",
            bbox_to_anchor=(0.97, 0.925), fontsize=5.4,
            frameon=False, ncol=len(strand_handles),
            handlelength=0.8, handleheight=0.6,
            handletextpad=0.3, columnspacing=0.75,
        )

    fig.suptitle(f"{chrom}  ({_fmt_bp(chrom_size)})",
                 fontsize=8, fontweight="bold", y=0.965, x=0.09, ha="left")
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
        _warn(f"No '*_report.tsv' file found in '{args.directory}'; overview classification donut will show no data.")

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
            lambda fig: fig.savefig(overview_path, dpi=args.dpi, bbox_inches="tight"),
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
                lambda fig, path=path: fig.savefig(path, dpi=args.dpi, bbox_inches="tight"),
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
                lambda fig: pdf.savefig(fig, bbox_inches="tight"),
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
                    lambda fig: pdf.savefig(fig, bbox_inches="tight"),
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
