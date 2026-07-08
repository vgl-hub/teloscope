#!/usr/bin/env python3
"""
plot_its.py — Single-locus figures for interstitial telomeric sequences (ITSs).

Interstitial telomeres (e.g. chromosome-fusion scars) have no end anchor, so they are
plotted as one square coordinate panel over a user-specified region, topped by a
whole-scaffold ideogram that marks where the zoom sits. Reuses the parsers,
track-drawing, styling and save/fallback machinery from teloscope_report.py without
modifying it.

Usage:
    python plot_its.py <output_directory> CHROM:START-END [CHROM:START-END ...]
    python plot_its.py <output_directory> CHROM                 # auto-center on largest ITS cluster
    python plot_its.py <output_directory> CHROM:S-E --canonical-matches --png

Coordinates accept commas/underscores and optional bp/kb/Mb suffixes.
"""

import os
import sys
import re
import glob
import argparse

import numpy as np

# Importing the report module first configures the Agg backend and Nature style.
from teloscope_report import (
    find_files, parse_terminal_bed, parse_bedgraph, parse_interval_bed,
    get_chrom_sizes, _draw_fraction_track, _draw_strand_track,
    _apply_terminal_x_axis, _save_figure_with_fallback, _fmt_bp,
    _sanitize_filename, _warn, _block_symbol_fits, _draw_block_symbol,
    _set_track_label, _hide_x_axis, COLORS, FIG_WIDTH_SINGLE, PANEL_LABEL_SIZE,
    LEGEND_TEXT_SIZE, FIGURE_SUMMARY_SIZE, AXIS_TICK_SIZE,
)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch, ConnectionPatch

# Orientation palette for interstitial blocks; balanced uses green per request.
ITS_ORIENT_COLORS = {"p": COLORS["p"], "q": COLORS["q"], "b": "#009E73"}
ITS_ORIENT_LABELS = (("p", "p (forward)"), ("q", "q (reverse)"), ("b", "balanced"))
ZOOM_COLOR = COLORS["Discordant"]  # red box/funnel marking the zoom region

# Every track is individually togglable; the scaffold ideogram is off by default.
ALL_TRACKS = ("ideogram", "blocks", "matches", "density", "canonical", "strand", "gc", "entropy")
DEFAULT_TRACKS = frozenset({"blocks", "density", "canonical", "strand", "gc", "entropy"})


# ---------------------------------------------------------------------------
# Region resolution
# ---------------------------------------------------------------------------

def _parse_coord(token):
    """Parse a coordinate token (commas/underscores + optional bp/kb/Mb suffix) to bp."""
    text = token.strip().replace(",", "").replace("_", "").replace(" ", "")
    m = re.fullmatch(r"([0-9]*\.?[0-9]+)(mbp|mb|kbp|kb|bp)?", text, flags=re.I)
    if not m:
        raise ValueError(f"invalid coordinate '{token}'")
    unit = (m.group(2) or "bp").lower()
    factor = {"bp": 1, "kb": 1_000, "kbp": 1_000, "mb": 1_000_000, "mbp": 1_000_000}[unit]
    return int(round(float(m.group(1)) * factor))


def _parse_region(spec, chrom_sizes):
    """Parse 'chrom:start-end' into a clamped (chrom, start, end); raise ValueError on bad input."""
    chrom, coords = spec.rsplit(":", 1)
    if "-" not in coords:
        raise ValueError(f"region '{spec}' is missing '-' between start and end")
    start_tok, end_tok = coords.split("-", 1)
    start, end = _parse_coord(start_tok), _parse_coord(end_tok)
    if chrom not in chrom_sizes:
        sample = ", ".join(list(chrom_sizes)[:5])
        raise ValueError(f"unknown sequence '{chrom}' (known: {sample} ...)")
    size = int(chrom_sizes[chrom])
    start = max(0, min(start, size))
    end = min(end, size)
    if end <= start:
        raise ValueError(f"region '{spec}' is empty after clamping to {size:,} bp")
    return chrom, start, end


def _auto_region_window(its_blocks_list, chrom_size, pad=None, merge_gap=50_000):
    """Window around the largest interstitial cluster on a scaffold; None if no ITS blocks."""
    if not its_blocks_list:
        return None
    blocks = sorted(its_blocks_list, key=lambda b: b["start"])
    clusters = []
    cur = {"start": blocks[0]["start"], "end": blocks[0]["end"], "bp": blocks[0]["length"]}
    for b in blocks[1:]:
        if b["start"] - cur["end"] <= merge_gap:
            cur["end"] = max(cur["end"], b["end"])
            cur["bp"] += b["length"]
        else:
            clusters.append(cur)
            cur = {"start": b["start"], "end": b["end"], "bp": b["length"]}
    clusters.append(cur)
    best = max(clusters, key=lambda c: c["bp"])
    span = best["end"] - best["start"]
    use_pad = pad if pad is not None else max(2 * span, 10_000)
    return max(0, best["start"] - use_pad), min(int(chrom_size), best["end"] + use_pad)


def _region_axis_spec(view_start, view_end):
    """Round genomic ticks for an absolute-position window; positions always in Mbp."""
    span_mb = max(view_end - view_start, 1) / 1e6
    dec = 1 if span_mb >= 2 else 2 if span_mb >= 0.2 else 3
    ticks = [t for t in ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(
        view_start, view_end) if view_start <= t <= view_end]
    if len(ticks) < 2:
        ticks = list(np.linspace(view_start, view_end, 5))
    return {
        "axis_start": view_start, "axis_end": view_end, "tick_pos": ticks,
        "tick_labels": [f"{t / 1e6:.{dec}f}" for t in ticks],
        "xlabel": "Position (Mbp)",
    }


# ---------------------------------------------------------------------------
# Track drawing
# ---------------------------------------------------------------------------

def _draw_ideogram(ax, chrom_size, view_start, view_end, blocks_list, its_blocks_list):
    """Whole-scaffold overview: terminal caps, all ITS ticks, and the zoom window boxed."""
    bar_y, bar_h = 0.60, 0.34
    bottom = bar_y - bar_h / 2.0
    ax.add_patch(Rectangle((0, bottom), chrom_size, bar_h, facecolor="#ececec",
                           edgecolor="#b9b9b9", linewidth=0.5, zorder=1))

    # Terminal telomeres as dark end caps (min width so they stay visible).
    cap_min = chrom_size * 0.004
    for b in blocks_list or []:
        w = max(b["end"] - b["start"], cap_min)
        x = min(b["start"], chrom_size - w) if b["start"] > chrom_size / 2.0 else b["start"]
        ax.add_patch(Rectangle((x, bottom), w, bar_h, facecolor=COLORS["terminal"],
                               edgecolor="none", zorder=2))

    # Every interstitial telomere on the scaffold as an orientation-colored tick.
    for b in its_blocks_list or []:
        xm = (b["start"] + b["end"]) / 2.0
        ax.plot([xm, xm], [bottom, bottom + bar_h],
                color=ITS_ORIENT_COLORS.get(b.get("label", ""), COLORS["its"]),
                linewidth=0.7, zorder=3)

    # Zoom window marker (min visible width).
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
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t / 1e6:.0f}" for t in ticks])
    ax.tick_params(axis="x", length=2, width=0.4, pad=1.0,
                   labelsize=AXIS_TICK_SIZE, colors="#666666")
    ax.set_xlabel("Scaffold (Mbp)", fontsize=AXIS_TICK_SIZE, color="#666666", labelpad=1.5)
    _set_track_label(ax, "Scaffold")
    return box_l, box_r, box_bottom


def _draw_its_blocks_track(ax, view_start, view_end, blocks_list, its_blocks_list,
                           gap_blocks_list, label="Blocks"):
    """Single-panel block track in genomic coordinates; ITS colored by orientation."""
    backbone_y = 0.5
    view_span = max(view_end - view_start, 1)
    ax.plot([view_start, view_end], [backbone_y, backbone_y],
            color="#d6d6d6", linewidth=0.9, solid_capstyle="round", zorder=1)

    def _draw(seq, color_fn, zorder):
        for b in seq or []:
            if b["end"] <= view_start or b["start"] >= view_end:
                continue
            ds, de = max(b["start"], view_start), min(b["end"], view_end)
            ax.add_patch(Rectangle((ds, backbone_y - 0.07), de - ds, 0.14,
                                   facecolor=color_fn(b), edgecolor="none",
                                   alpha=0.98, zorder=zorder))
            if _block_symbol_fits(de - ds, view_span, b.get("label", "")):
                _draw_block_symbol(ax, (ds + de) / 2, backbone_y, b.get("label", ""), zorder + 1)

    _draw(blocks_list, lambda b: COLORS["terminal"], 2)
    _draw(its_blocks_list, lambda b: ITS_ORIENT_COLORS.get(b.get("label", ""), COLORS["its"]), 3)
    _draw(gap_blocks_list, lambda b: COLORS["gap"], 5)

    ax.set_ylim(0.28, 0.72)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    _set_track_label(ax, label)
    _hide_x_axis(ax)
    return True


def _draw_matches_track(ax, view_start, view_end, matches_list, color, label="Canonical\nmatches"):
    """Thin tick track marking individual repeat-match positions in the window."""
    xs = [(m["start"] + m["end"]) / 2.0 for m in (matches_list or [])
          if m["end"] > view_start and m["start"] < view_end]
    ax.plot([view_start, view_end], [0.5, 0.5], color="#e6e6e6", linewidth=0.6, zorder=1)
    if xs:
        ax.vlines(xs, 0.2, 0.8, color=color, linewidth=0.5, alpha=0.85, zorder=2)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    _set_track_label(ax, label)
    _hide_x_axis(ax)
    return True


def plot_region_zoom(chrom, chrom_size, view_start, view_end, blocks_list,
                     its_blocks_list, density_data=None, canonical_data=None,
                     strand_data=None, gc_data=None, entropy_data=None,
                     gap_blocks_list=None, matches_list=None, enabled=DEFAULT_TRACKS):
    """Square single-locus figure of a zoomed region; each track is individually togglable."""
    def _has(bg):
        return bg is not None and len(bg[0]) > 0

    show_ideogram = "ideogram" in enabled
    zoom_specs = []
    if "blocks" in enabled:
        zoom_specs.append(("blocks", None))
    if "matches" in enabled and matches_list is not None:
        zoom_specs.append(("matches", matches_list))
    if "density" in enabled and _has(density_data):
        zoom_specs.append(("density", density_data))
    if "canonical" in enabled and _has(canonical_data):
        zoom_specs.append(("canonical", canonical_data))
    if "strand" in enabled and _has(strand_data):
        zoom_specs.append(("strand", strand_data))
    if "gc" in enabled and _has(gc_data):
        zoom_specs.append(("gc", gc_data))
    if "entropy" in enabled and _has(entropy_data):
        zoom_specs.append(("entropy", entropy_data))

    height_map = {"blocks": 0.20, "matches": 0.16, "density": 0.34, "canonical": 0.34,
                  "strand": 0.34, "gc": 0.34, "entropy": 0.34}
    height_ratios = ([0.30] if show_ideogram else []) + [height_map[n] for n, _ in zoom_specs]
    n_rows = len(height_ratios)
    if n_rows == 0:
        raise ValueError("no tracks enabled to plot")

    fig_h = 0.9 + 0.40 * n_rows
    fig, axes = plt.subplots(
        n_rows, 1, figsize=(FIG_WIDTH_SINGLE + 0.45, fig_h),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.42}, squeeze=False)
    # Reserve a fixed ~0.66 in header regardless of row count so the title block never crowds.
    fig.subplots_adjust(left=0.215, right=0.965, top=1 - 0.66 / fig_h, bottom=0.10)

    row_offset = 0
    if show_ideogram:
        ideo_ax = axes[0][0]
        box_l, box_r, box_bottom = _draw_ideogram(
            ideo_ax, chrom_size, view_start, view_end, blocks_list, its_blocks_list)
        row_offset = 1

    axis_spec = _region_axis_spec(view_start, view_end)
    visible_rows = []
    for i, (track_name, track_data) in enumerate(zoom_specs):
        row = row_offset + i
        ax = axes[row][0]
        if track_name == "blocks":
            visible = _draw_its_blocks_track(ax, view_start, view_end, blocks_list,
                                             its_blocks_list, gap_blocks_list)
        elif track_name == "matches":
            visible = _draw_matches_track(ax, view_start, view_end, track_data,
                                          COLORS["canonical"])
        elif track_name == "density":
            visible = _draw_fraction_track(ax, track_data, view_start, view_end, chrom_size,
                                           "full", COLORS["density"], label="Repeat\ndensity")
        elif track_name == "canonical":
            visible = _draw_fraction_track(ax, track_data, view_start, view_end, chrom_size,
                                           "full", COLORS["canonical"], label="Canonical\nratio")
        elif track_name == "gc":
            starts, ends, values = track_data
            visible = _draw_fraction_track(ax, (starts, ends, values / 100.0), view_start,
                                           view_end, chrom_size, "full", COLORS["gc"],
                                           label="GC\ncontent", y_max=1.0, y_min=0.0)
        elif track_name == "entropy":
            visible = _draw_fraction_track(ax, track_data, view_start, view_end, chrom_size,
                                           "full", COLORS["entropy"], label="Shannon\nentropy",
                                           y_max=2.0, y_min=1.5)
        else:
            visible = _draw_strand_track(ax, track_data, view_start, view_end, chrom_size,
                                         "full", label="Strand\nbias")
        if visible:
            ax.set_xlim(axis_spec["axis_start"], axis_spec["axis_end"])
            visible_rows.append(row)

    if visible_rows:
        for row in visible_rows[:-1]:
            _hide_x_axis(axes[row][0])
        _apply_terminal_x_axis(axes[visible_rows[-1]][0], axis_spec, "full")

    # Zoom funnel from the ideogram window box down to the first zoomed panel (when shown).
    if show_ideogram and visible_rows:
        target_ax = axes[visible_rows[0]][0]
        for x_ideo, x_ax in ((box_l, 0.0), (box_r, 1.0)):
            fig.add_artist(ConnectionPatch(
                xyA=(x_ideo, box_bottom), coordsA=ideo_ax.transData,
                xyB=(x_ax, 1.0), coordsB=target_ax.transAxes,
                color=ZOOM_COLOR, linewidth=0.5, linestyle=(0, (3, 2)), alpha=0.55, zorder=0))

    # Legend: only orientations/annotations actually present within the plotted window.
    def _in_view(b):
        return b["end"] > view_start and b["start"] < view_end
    present = {b.get("label", "") for b in (its_blocks_list or []) if _in_view(b)}
    handles = [Patch(facecolor=ITS_ORIENT_COLORS[k], label=name)
               for k, name in ITS_ORIENT_LABELS if k in present]
    if any(_in_view(b) for b in (blocks_list or [])):
        handles.append(Patch(facecolor=COLORS["terminal"], label="terminal"))
    if any(_in_view(b) for b in (gap_blocks_list or [])):
        handles.append(Patch(facecolor=COLORS["gap"], label="gap"))
    if handles:
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1 - 0.49 / fig_h),
                   ncol=min(len(handles), 4), fontsize=LEGEND_TEXT_SIZE, frameon=False,
                   handlelength=0.9, handletextpad=0.3, columnspacing=1.0)

    fig.suptitle(f"{chrom}:{view_start:,}-{view_end:,}", fontsize=PANEL_LABEL_SIZE,
                 fontweight="bold", y=1 - 0.14 / fig_h, x=0.5, ha="center", va="top")
    fig.text(0.5, 1 - 0.31 / fig_h, f"Scaffold size: {_fmt_bp(chrom_size)}",
             ha="center", va="top", fontsize=FIGURE_SUMMARY_SIZE, color="#444444")
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot interstitial telomeric sequences (ITSs) for specified regions.",
        epilog="Example: python plot_its.py output/ NC_052576.1:19250000-19280000",
    )
    parser.add_argument("directory", help="Teloscope output directory")
    parser.add_argument("regions", nargs="+",
                        help="CHROM:START-END for an exact window, or bare CHROM to "
                             "auto-center on its largest interstitial cluster")
    parser.add_argument("-o", "--output", default=None,
                        help="Output PDF path (or directory with --png)")
    parser.add_argument("--png", action="store_true",
                        help="Save individual PNG files instead of a single PDF")
    parser.add_argument("--canonical-matches", action="store_true",
                        help="Add a track of individual canonical repeat-match positions")
    parser.add_argument("--pad", type=int, default=None,
                        help="Symmetric padding in bp (default: auto for cluster windows)")
    parser.add_argument("--merge-gap", type=int, default=50_000,
                        help="Max gap (bp) to merge ITS blocks into one cluster for auto-centering")
    parser.add_argument("--dpi", type=int, default=450, help="DPI for raster output")
    parser.add_argument("--draft", action="store_true", help="Render at 150 DPI for fast iteration")
    parser.add_argument("--show", nargs="+", choices=ALL_TRACKS, default=[], metavar="TRACK",
                        help="Turn on tracks that are off by default (e.g. ideogram)")
    parser.add_argument("--hide", nargs="+", choices=ALL_TRACKS, default=[], metavar="TRACK",
                        help="Turn off tracks that are on by default")
    args = parser.parse_args()
    if args.draft:
        args.dpi = 150

    enabled = set(DEFAULT_TRACKS)
    if args.canonical_matches:
        enabled.add("matches")
    enabled |= set(args.show)
    enabled -= set(args.hide)

    if not os.path.isdir(args.directory):
        sys.exit(f"Error: '{args.directory}' is not a directory.")

    files = find_files(args.directory)
    if "interstitial" not in files:
        sys.exit(f"Error: no '*_interstitial_telomeres.bed' file in '{args.directory}'.")

    its_blocks = parse_terminal_bed(files["interstitial"])
    blocks = parse_terminal_bed(files["terminal"]) if "terminal" in files else {}
    density = parse_bedgraph(files["density"]) if "density" in files else None
    canonical = parse_bedgraph(files["canonical_ratio"]) if "canonical_ratio" in files else None
    strand = parse_bedgraph(files["strand_ratio"]) if "strand_ratio" in files else None
    gc = parse_bedgraph(files["gc"]) if "gc" in files else None
    entropy = parse_bedgraph(files["entropy"]) if "entropy" in files else None
    gaps = parse_interval_bed(files["gaps"]) if "gaps" in files else None

    matches = None
    if "matches" in enabled:
        hits = sorted(glob.glob(os.path.join(args.directory, "*_canonical_matches.bed")))
        if hits:
            matches = parse_interval_bed(hits[0])
        else:
            _warn("No '*_canonical_matches.bed' found; skipping canonical-matches track.")

    merged = {}
    for dataset in (blocks, its_blocks):
        for chrom, blist in (dataset or {}).items():
            merged.setdefault(chrom, []).extend(blist)
    chrom_sizes = get_chrom_sizes(merged, density, canonical, strand, gc, entropy)

    windows, seen = [], set()
    for spec in args.regions:
        try:
            if ":" in spec:
                chrom, start, end = _parse_region(spec, chrom_sizes)
                if args.pad:
                    start = max(0, start - args.pad)
                    end = min(int(chrom_sizes[chrom]), end + args.pad)
            else:
                chrom = spec
                if chrom not in chrom_sizes:
                    sample = ", ".join(list(chrom_sizes)[:5])
                    raise ValueError(f"unknown sequence '{chrom}' (known: {sample} ...)")
                win = _auto_region_window(its_blocks.get(chrom, []), chrom_sizes[chrom],
                                          pad=args.pad, merge_gap=args.merge_gap)
                if win is None:
                    raise ValueError(f"no interstitial telomeres on '{chrom}'; give CHROM:START-END")
                start, end = win
        except ValueError as exc:
            _warn(f"Skipping region '{spec}': {exc}")
            continue
        key = (chrom, start, end)
        if key in seen:
            _warn(f"Skipping duplicate region {chrom}:{start:,}-{end:,}.")
            continue
        seen.add(key)
        windows.append(key)

    if not windows:
        sys.exit("Error: no valid regions to plot.")

    def _region_fig(chrom, start, end):
        return plot_region_zoom(
            chrom, int(chrom_sizes[chrom]), start, end,
            blocks.get(chrom, []), its_blocks.get(chrom, []),
            density.get(chrom) if density else None,
            canonical.get(chrom) if canonical else None,
            strand.get(chrom) if strand else None,
            gc.get(chrom) if gc else None,
            entropy.get(chrom) if entropy else None,
            gaps.get(chrom) if gaps else None,
            matches.get(chrom, []) if matches is not None else None,
            enabled=enabled,
        )

    if args.png:
        out_dir = args.output or args.directory
        os.makedirs(out_dir, exist_ok=True)
        for chrom, start, end in windows:
            stem = _sanitize_filename(f"{chrom}_{start}-{end}")
            path = os.path.join(out_dir, f"teloscope_its_{stem}.png")
            ok, _ = _save_figure_with_fallback(
                lambda fig, path=path: fig.savefig(path, dpi=args.dpi),
                lambda chrom=chrom, start=start, end=end: _region_fig(chrom, start, end),
                f"{chrom}:{start:,}-{end:,}",
                "Failed to render this interstitial region. A placeholder image was written instead.",
            )
            print(f"{'[warning] ' if not ok else ''}{path}", file=sys.stderr)
    else:
        out_path = args.output or os.path.join(args.directory, "teloscope_its.pdf")
        with PdfPages(out_path) as pdf:
            for chrom, start, end in windows:
                ok, _ = _save_figure_with_fallback(
                    lambda fig: pdf.savefig(fig),
                    lambda chrom=chrom, start=start, end=end: _region_fig(chrom, start, end),
                    f"{chrom}:{start:,}-{end:,}",
                    "Failed to render this interstitial region. A placeholder page was written instead.",
                )
                print(f"{'[warning] ' if not ok else ''}{chrom}:{start:,}-{end:,}", file=sys.stderr)
        print(f"Saved {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
