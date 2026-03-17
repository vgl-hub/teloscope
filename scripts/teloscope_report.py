#!/usr/bin/env python3
"""
teloscope_report.py — Publication-ready figures from Teloscope output.

Usage:
    python teloscope_report.py <output_directory> [-o report.pdf]
    python teloscope_report.py <output_directory> --png

Reads Teloscope output files from the given directory and generates:
  Page 1: Assembly overview (classification donut + telomere length distribution)
  Page 2+: Per-chromosome telomere profiles (blocks, density, strand ratio)

Requires: Python 3.6+, matplotlib
"""

import sys
import os
import glob
import argparse
from collections import defaultdict, OrderedDict

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
        "savefig.dpi":        600,
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
    Parse a BEDgraph file -> OrderedDict[chrom -> list of (start, end, value)].
    """
    data = OrderedDict()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            if chrom not in data:
                data[chrom] = []
            data[chrom].append((int(parts[1]), int(parts[2]), float(parts[3])))
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
        for chrom, intervals in bedgraph_data.items():
            if intervals:
                sizes[chrom] = max(sizes.get(chrom, 0), max(e for _, e, _ in intervals))
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
    """Proper median calculation."""
    s = sorted(values)
    n = len(s)
    if n % 2 == 1:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2


def _bedgraph_to_step(intervals):
    """Convert BEDgraph intervals to step-plot coordinates."""
    xs, ys = [], []
    for start, end, val in intervals:
        xs.append(start)
        ys.append(max(val, 0))
    if intervals:
        xs.append(intervals[-1][1])
        ys.append(0)
    return xs, ys


def _blend_color(hex1, hex2, frac):
    """Blend hex1 and hex2 by *frac* (1=hex1, 0=hex2)."""
    r1, g1, b1 = int(hex1[1:3], 16), int(hex1[3:5], 16), int(hex1[5:7], 16)
    r2, g2, b2 = int(hex2[1:3], 16), int(hex2[3:5], 16), int(hex2[5:7], 16)
    r = int(r1 * frac + r2 * (1 - frac))
    g = int(g1 * frac + g2 * (1 - frac))
    b = int(b1 * frac + b2 * (1 - frac))
    return f"#{r:02x}{g:02x}{b:02x}"


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
            ax_hist.hist(
                data_stack, bins=bins, stacked=True,
                color=colors_stack, label=labels_stack,
                edgecolor="white", linewidth=0.4, alpha=0.85,
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
# Tier 2: Per-chromosome profile
# ---------------------------------------------------------------------------

def plot_chromosome_profile(chrom, chrom_size, blocks_list,
                            density_data=None, strand_data=None):
    """Multi-track chromosome plot: blocks, repeat density, strand ratio."""
    has_density = density_data is not None and len(density_data) > 0
    has_strand = strand_data is not None and len(strand_data) > 0

    n_tracks = 1
    height_ratios = [0.35]
    if has_density:
        n_tracks += 1
        height_ratios.append(1.0)
    if has_strand:
        n_tracks += 1
        height_ratios.append(0.5)

    fig_height = 0.8 + 0.9 * (n_tracks - 1)
    fig, axes = plt.subplots(
        n_tracks, 1, figsize=(FIG_WIDTH_DOUBLE, fig_height),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.08},
        sharex=True,
    )
    if n_tracks == 1:
        axes = [axes]

    fig.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.18)
    ax_idx = 0

    # ---- Track: Telomere blocks ----
    ax_blocks = axes[ax_idx]
    ax_idx += 1

    backbone_y = 0.5
    ax_blocks.plot([0, chrom_size], [backbone_y, backbone_y],
                   color="#d0d0d0", linewidth=2.5, solid_capstyle="round", zorder=1)

    for b in blocks_list:
        color = COLORS.get(b["label"], "#999999")
        rect = Rectangle(
            (b["start"], backbone_y - 0.3), b["length"], 0.6,
            facecolor=color, edgecolor="none", alpha=0.85, zorder=2,
        )
        ax_blocks.add_patch(rect)

    ax_blocks.set_ylim(-0.1, 1.1)
    ax_blocks.set_ylabel("Blocks", fontsize=7)
    ax_blocks.set_yticks([])
    ax_blocks.spines["bottom"].set_visible(False)
    ax_blocks.spines["left"].set_visible(False)
    ax_blocks.tick_params(bottom=False)

    labels_present = set(b["label"] for b in blocks_list)
    legend_elements = []
    for lab, name in [("p", "p-arm (fwd)"), ("q", "q-arm (rev)"), ("b", "balanced")]:
        if lab in labels_present:
            legend_elements.append(
                Line2D([0], [0], marker="s", color="w", markerfacecolor=COLORS[lab],
                       markersize=5, linestyle="none", label=name))
    if legend_elements:
        ax_blocks.legend(
            handles=legend_elements, loc="upper right", fontsize=6,
            frameon=False, ncol=len(legend_elements), handletextpad=0.2,
        )

    # ---- Track: Repeat density ----
    if has_density:
        ax_den = axes[ax_idx]
        ax_idx += 1
        xs, ys = _bedgraph_to_step(density_data)
        ax_den.fill_between(xs, ys, step="post", color=COLORS["density"],
                            alpha=0.35, linewidth=0)
        ax_den.plot(xs, ys, drawstyle="steps-post", color=COLORS["density"],
                    linewidth=0.6)
        ax_den.set_ylabel("Repeat\ndensity", fontsize=7)
        ax_den.set_ylim(0, 1.05)
        ax_den.tick_params(bottom=False)
        ax_den.spines["bottom"].set_visible(False)

    # ---- Track: Strand ratio ----
    if has_strand:
        ax_str = axes[ax_idx]
        ax_idx += 1
        for start, end, val in strand_data:
            if val < 0:
                ax_str.barh(0.5, end - start, left=start, height=1.0,
                            color=COLORS["no_data"], edgecolor="none", alpha=0.3)
            else:
                color = _blend_color(COLORS["fwd"], COLORS["rev"], val)
                ax_str.barh(0.5, end - start, left=start, height=1.0,
                            color=color, edgecolor="none", alpha=0.7)

        ax_str.set_ylim(0, 1)
        ax_str.set_yticks([])
        ax_str.set_ylabel("Strand\nratio", fontsize=7)

        ax_str.legend(
            handles=[
                Patch(facecolor=COLORS["fwd"], label="Forward"),
                Patch(facecolor="#808080", label="Mixed"),
                Patch(facecolor=COLORS["rev"], label="Reverse"),
                Patch(facecolor=COLORS["no_data"], alpha=0.6, label="No data"),
            ],
            loc="upper right", fontsize=5.5, frameon=False, ncol=4,
            handlelength=0.8, handleheight=0.7, handletextpad=0.2,
        )

    axes[-1].set_xlim(0, chrom_size)
    _format_bp_axis(axes[-1], chrom_size)

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
    parser.add_argument("--dpi", type=int, default=600,
                        help="DPI for raster output (default: 600)")
    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        sys.exit(f"Error: '{args.directory}' is not a directory.")

    files = find_files(args.directory)
    if "terminal" not in files:
        sys.exit(f"Error: No *_terminal_telomeres.bed found in '{args.directory}'.\n"
                 f"Run teloscope first to generate output files.")

    print(f"Found files: {', '.join(files.keys())}")

    blocks = parse_terminal_bed(files["terminal"])
    density_data = parse_bedgraph(files["density"]) if "density" in files else None
    strand_data  = parse_bedgraph(files["strand_ratio"]) if "strand_ratio" in files else None

    all_chroms = list(density_data.keys()) if density_data else None
    chrom_sizes = get_chrom_sizes(blocks, density_data)

    classifications = classify_chromosomes(blocks, all_chroms)
    total_chroms = sum(len(v) for v in classifications.values())
    total_telo = sum(len(blist) for blist in blocks.values())
    print(f"Chromosomes: {total_chroms}  |  Telomere blocks: {total_telo}  |  "
          f"Categories: {', '.join(f'{k}={len(v)}' for k, v in classifications.items())}")

    profile_chroms = sorted(blocks.keys(),
                            key=lambda c: chrom_sizes.get(c, 0), reverse=True)
    if density_data:
        no_telo = [c for c in density_data.keys() if c not in blocks]
        no_telo.sort(key=lambda c: chrom_sizes.get(c, 0), reverse=True)
        if len(profile_chroms) + len(no_telo) <= 40:
            profile_chroms.extend(no_telo)
        else:
            profile_chroms.extend(no_telo[:5])

    figures = []
    fig_overview = plot_assembly_overview(classifications, blocks, chrom_sizes)
    figures.append(("overview", fig_overview))

    for chrom in profile_chroms:
        csize = chrom_sizes.get(chrom, 0)
        if csize == 0:
            continue
        blist = blocks.get(chrom, [])
        den = density_data.get(chrom) if density_data else None
        strand = strand_data.get(chrom) if strand_data else None
        if not blist and den is None:
            continue
        fig_chr = plot_chromosome_profile(chrom, csize, blist, den, strand)
        figures.append((chrom, fig_chr))

    if not figures:
        sys.exit("No figures generated.")

    if args.png:
        out_dir = args.output or args.directory
        os.makedirs(out_dir, exist_ok=True)
        for name, fig in figures:
            safe_name = name.replace("/", "_").replace("\\", "_")
            path = os.path.join(out_dir, f"teloscope_{safe_name}.png")
            fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
            print(f"  Saved {path}")
            plt.close(fig)
    else:
        out_path = args.output or os.path.join(args.directory, "teloscope_report.pdf")
        with PdfPages(out_path) as pdf:
            for name, fig in figures:
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
        print(f"Report saved to {out_path}")


if __name__ == "__main__":
    main()
