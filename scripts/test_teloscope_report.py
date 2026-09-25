import atexit
import contextlib
import importlib.util
import io
import os
from collections import OrderedDict
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock
import sys

ROOT = Path(__file__).resolve().parents[1]
MPLCONFIGDIR = tempfile.mkdtemp(prefix="mplconfig_teloscope_test_")
os.environ["MPLCONFIGDIR"] = MPLCONFIGDIR
atexit.register(shutil.rmtree, MPLCONFIGDIR, ignore_errors=True)

from matplotlib.legend import Legend
import matplotlib.text as mtext
import numpy as np


MODULE_PATH = ROOT / "scripts" / "teloscope_report.py"
SPEC = importlib.util.spec_from_file_location("teloscope_report_under_test", MODULE_PATH)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def _its_row(chrom, start, end, label, telo_type, fwd_can=0, rev_can=0,
             fwd_noncan=0, rev_noncan=0, chrom_size=100_000, closest_end=None):
    closest_end = closest_end if closest_end is not None else ("q" if label == "q" else "p")
    telo_len = end - start
    return (f"{chrom}\t{start}\t{end}\t{telo_len}\t{label}\t{closest_end}\t"
            f"{fwd_can}\t{rev_can}\t{fwd_noncan}\t{rev_noncan}\t{chrom_size}\t{telo_type}\n")


def _its_frame(rows, tmpdir):
    bed_path = Path(tmpdir) / "its.bed"
    bed_path.write_text("".join(rows), encoding="utf-8")
    return REPORT.load_its_frame(str(bed_path))


def _gaps_frame(rows, tmpdir):
    path = Path(tmpdir) / "gaps.bed"
    path.write_text("".join(rows), encoding="utf-8")
    return REPORT.load_gaps_frame(str(path))


def _block(start, end, label, chrom_size, closest_end=None):
    start, end = int(start), int(end)
    return {
        "start": start,
        "end": end,
        "span": end - start,
        "length": end - start,
        "teloLen": end - start,
        "label": label,
        "closestEnd": closest_end if closest_end is not None else label,
        "pathSize": int(chrom_size),
    }


def _synthetic_terminal_dataset(chrom_size=20_000):
    chrom = "chrSynthetic"
    blocks = {
        chrom: [
            _block(0, 6_600, "p", chrom_size),
            _block(chrom_size - 6_600, chrom_size, "q", chrom_size),
        ]
    }
    chrom_sizes = {chrom: chrom_size}
    starts = np.arange(0, chrom_size, 1_000, dtype=np.int64)
    ends = np.minimum(starts + 1_000, chrom_size)
    values = np.linspace(0.2, 0.8, len(starts), dtype=np.float64)
    bedgraph = {chrom: (starts, ends, values)}
    return chrom, blocks, chrom_sizes, bedgraph


class TeloscopeReportTests(unittest.TestCase):
    def test_orient_kw_returns_vert_for_old_matplotlib(self):
        def fake_violinplot(**kwargs):
            pass
        fake_violinplot.__name__ = "violinplot_old"
        kw = REPORT._orient_kw(fake_violinplot, True)
        self.assertEqual(kw, {"vert": True})
        kw = REPORT._orient_kw(fake_violinplot, False)
        self.assertEqual(kw, {"vert": False})

    def test_orient_kw_returns_orientation_for_new_matplotlib(self):
        def fake_boxplot(**kwargs):
            pass
        fake_boxplot.__name__ = "boxplot_new"
        import inspect
        sig = inspect.signature(lambda orientation=None: None)
        REPORT._ORIENT_KW_CACHE["boxplot_new"] = "orientation" in sig.parameters
        kw = REPORT._orient_kw(fake_boxplot, True)
        self.assertIn("orientation", kw)
        self.assertEqual(kw["orientation"], "vertical")

    def test_parse_terminal_bed_reads_the_twelve_column_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "blocks.bed"
            bed_path.write_text(
                "chrTerm\t10\t70\t60\tp\tp\t3\t2\t5\t7\t100\tscaffold\n"
                "chrIts\t200\t260\t60\tb\tq\t4\t6\t1\t2\t500\tfusion\n",
                encoding="utf-8",
            )

            parsed = REPORT.parse_terminal_bed(str(bed_path))

        term = parsed["chrTerm"][0]
        self.assertEqual(term["length"], 60)
        self.assertEqual((term["label"], term["closestEnd"]), ("p", "p"))
        self.assertEqual((term["fwd"], term["rev"]), (8, 9))
        self.assertEqual((term["can"], term["noncan"]), (5, 12))
        self.assertEqual(
            (term["fwdCan"], term["revCan"], term["fwdNonCan"], term["revNonCan"]),
            (3, 2, 5, 7),
        )
        self.assertEqual((term["pathSize"], term["term"]), (100, "scaffold"))

        its = parsed["chrIts"][0]
        self.assertEqual((its["label"], its["closestEnd"]), ("b", "q"))
        self.assertEqual((its["pathSize"], its["term"]), (500, "fusion"))

    def test_parse_terminal_bed_skips_rows_with_the_wrong_column_count(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "blocks.bed"
            bed_path.write_text(
                "chrBroken\t20\t80\t60\tq\tq\t4\t6\t7\t3\t100\n",
                encoding="utf-8",
            )

            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                parsed = REPORT.parse_terminal_bed(str(bed_path))

        self.assertEqual(dict(parsed), {})
        self.assertIn("expected 12 BED columns", stderr.getvalue())

    def test_discordant_row_is_placed_at_its_closest_end_not_its_label(self):
        chrom_size = 20_000
        discordant_block = _block(0, 600, "q", chrom_size, closest_end="p")

        self.assertEqual(REPORT._block_end_distance(discordant_block, chrom_size), 0)

        # A p-window tight around the block, not spanning the whole chromosome (the old label="q" bug).
        p_window, q_window = REPORT.compute_view_windows([discordant_block], chrom_size)
        self.assertEqual(p_window, (0, 2_000))
        self.assertEqual(q_window, (18_000, 20_000))

    def test_fragmented_row_reports_length_as_teloLen_not_the_span(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "blocks.bed"
            bed_path.write_text(
                "chrFrag\t0\t1000\t600\tp\tp\t3\t2\t5\t7\t100000\tscaffold\n",
                encoding="utf-8",
            )

            parsed = REPORT.parse_terminal_bed(str(bed_path))

        self.assertEqual(parsed["chrFrag"][0]["length"], 600)

    def test_discordant_row_uses_the_reverse_glyph_and_sits_at_the_p_end(self):
        chrom_size = 20_000
        discordant_block = _block(0, 600, "q", chrom_size, closest_end="p")

        fig = REPORT.plot_terminal_zoom("chrDiscordant", chrom_size, [discordant_block])
        self.addCleanup(REPORT.plt.close, fig)

        self.assertEqual(fig.axes[0].get_title(), "p-arm")
        glyphs = {t.get_text() for t in fig.findobj(mtext.Text)} & {"<", ">"}
        self.assertEqual(glyphs, {">"})

    def test_balanced_row_is_not_dropped_from_the_p_arm_window(self):
        chrom_size = 20_000
        balanced_block = _block(0, 600, "b", chrom_size, closest_end="p")

        p_window, q_window = REPORT.compute_view_windows([balanced_block], chrom_size)
        self.assertEqual(p_window, (0, 2_000))

    def test_parse_report_maps_expected_categories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "synthetic_report.tsv"
            report_path.write_text(
                "+++\n"
                "pos\theader\ttelomeres\tlabels\tgaps\ttype\tanomaly\tgranular\n"
                "1\tchr_dis\t1\tq\t0\tincomplete\tdiscordant_q\tQ*\n"
                "2\tchr_frag\t1\tp\t0\tincomplete\tfragmented_p\tPp\n"
                "3\tchr_none\t0\tnone\t0\tnone\t.\t\n",
                encoding="utf-8",
            )

            parsed = REPORT.parse_report(str(report_path))

        self.assertEqual(parsed["Discordant"], ["chr_dis"])
        self.assertEqual(parsed["Fragmented"], ["chr_frag"])
        self.assertEqual(parsed["No telomeres"], ["chr_none"])
        self.assertNotIn("Misassembly", parsed)
        self.assertNotIn("Balanced", parsed)

    def test_compute_view_windows_rounds_to_displayed_kbp(self):
        chrom_size = 20_000
        blocks_list = [
            _block(0, 6_600, "p", chrom_size),
            _block(chrom_size - 6_600, chrom_size, "q", chrom_size),
        ]

        p_window, q_window = REPORT.compute_view_windows(blocks_list, chrom_size)

        self.assertEqual(p_window, (0, 14_000))
        self.assertEqual(q_window, (6_000, 20_000))

    def test_terminal_backbone_reaches_axis_limit(self):
        chrom, blocks, chrom_sizes, bedgraph = _synthetic_terminal_dataset()
        fig = REPORT.plot_terminal_zoom(
            chrom,
            chrom_sizes[chrom],
            blocks[chrom],
            density_data=bedgraph[chrom],
            canonical_data=bedgraph[chrom],
            strand_data=bedgraph[chrom],
        )
        self.addCleanup(REPORT.plt.close, fig)

        blocks_axis = fig.axes[0]
        xlim = tuple(float(v) for v in blocks_axis.get_xlim())
        backbone = blocks_axis.lines[0]
        line_min = float(np.min(backbone.get_xdata()))
        line_max = float(np.max(backbone.get_xdata()))

        self.assertEqual((line_min, line_max), (min(xlim), max(xlim)))

    def test_terminal_track_labels_right_align_multiline_rows(self):
        chrom, blocks, chrom_sizes, bedgraph = _synthetic_terminal_dataset()
        fig = REPORT.plot_terminal_zoom(
            chrom,
            chrom_sizes[chrom],
            blocks[chrom],
            density_data=bedgraph[chrom],
            canonical_data=bedgraph[chrom],
            strand_data=bedgraph[chrom],
        )
        self.addCleanup(REPORT.plt.close, fig)
        fig.canvas.draw()

        renderer = fig.canvas.get_renderer()
        for label_text in ("Repeat\ndensity", "Canonical\nratio", "Strand\nbias"):
            label = next(ax.yaxis.label for ax in fig.axes if ax.yaxis.label.get_text() == label_text)
            _, lines, _ = label._get_layout(renderer)
            right_edges = [float(x + size[0]) for _, size, x, _ in lines]

            self.assertEqual(label.get_ha(), "right")
            self.assertLess(max(right_edges) - min(right_edges), 0.1)

    def test_terminal_gap_intervals_are_rendered_as_light_gray_bars(self):
        chrom_size = 3_100
        fig = REPORT.plot_terminal_zoom(
            "chrGap",
            chrom_size,
            [
                _block(0, 600, "p", chrom_size),
                _block(2_500, 3_100, "q", chrom_size),
            ],
            gap_blocks_list=[{"start": 1_500, "end": 1_600}],
        )
        self.addCleanup(REPORT.plt.close, fig)

        gap_patches = [
            patch for patch in fig.axes[0].patches
            if float(patch.get_x()) == 1500.0 and float(patch.get_width()) == 100.0
        ]
        self.assertEqual(len(gap_patches), 1)
        self.assertEqual(gap_patches[0].get_facecolor(), (0.8392156862745098, 0.8392156862745098, 0.8392156862745098, 0.98))

    def test_overview_legends_form_centered_shared_group(self):
        blocks = {
            "chrA": [
                _block(0, 600, "p", 3_200),
                _block(2_600, 3_200, "q", 3_200),
            ]
        }
        chrom_sizes = {
            "chrA": 3_200,
            "chrB": 6_400,
            "chrC": 12_000,
            "chrD": 8_000,
            "chrE": 4_000,
        }
        classifications = OrderedDict(
            [
                ("T2T", ["chrA"]),
                ("Incomplete", ["chrB"]),
                ("Fragmented", ["chrC"]),
                ("Discordant", ["chrD"]),
                ("No telomeres", ["chrE"]),
            ]
        )

        fig = REPORT.plot_overview_page1(classifications, blocks, chrom_sizes)
        self.addCleanup(REPORT.plt.close, fig)
        fig.canvas.draw()

        legends = sorted(
            fig.findobj(Legend),
            key=lambda legend: legend.get_window_extent(fig.canvas.get_renderer()).x0,
        )
        self.assertEqual(len(legends), 2)

        legend_boxes = [
            legend.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
            for legend in legends
        ]
        legend_gap = legend_boxes[1].x0 - legend_boxes[0].x1
        self.assertGreaterEqual(legend_gap, 0.0)
        self.assertLess(legend_gap, 0.03)

        legend_group_center = (legend_boxes[0].x0 + legend_boxes[1].x1) / 2.0
        legend_axis = min(fig.axes, key=lambda ax: ax.get_position().y0)
        legend_axis_center = (legend_axis.get_position().x0 + legend_axis.get_position().x1) / 2.0
        self.assertAlmostEqual(legend_group_center, legend_axis_center, delta=0.02)

    def test_overview_summary_separates_scaffold_and_block_flag_counts(self):
        blocks = {
            "chrNear": [_block(0, 600, "p", 5_000)],
            "chrFar": [_block(1_400, 2_000, "p", 5_000)],
        }
        chrom_sizes = {
            "chrNear": 5_000,
            "chrFar": 5_000,
            "chrOther": 6_000,
        }
        classifications = OrderedDict(
            [
                ("Fragmented", ["chrNear", "chrMissing"]),
                ("Discordant", ["chrFar"]),
                ("No telomeres", ["chrOther"]),
            ]
        )

        fig = REPORT.plot_overview_page1(classifications, blocks, chrom_sizes)
        self.addCleanup(REPORT.plt.close, fig)

        summary_text = next(text.get_text() for text in fig.texts if text.get_text().startswith("Scaffolds ("))
        self.assertIn("Scaffolds (n=4, flagged=3)", summary_text)
        self.assertIn("telomere blocks (n=2, distance-flagged=1)", summary_text)
        self.assertNotIn("Scaffolds (n=4, distance-flagged=", summary_text)

    def test_text_floor_and_block_glyphs_are_nature_compliant(self):
        chrom, blocks, chrom_sizes, bedgraph = _synthetic_terminal_dataset()
        fig = REPORT.plot_terminal_zoom(
            chrom,
            chrom_sizes[chrom],
            blocks[chrom],
            density_data=bedgraph[chrom],
            canonical_data=bedgraph[chrom],
            strand_data=bedgraph[chrom],
        )
        self.addCleanup(REPORT.plt.close, fig)

        text_strings = []
        font_sizes = []
        for artist in fig.findobj(mtext.Text):
            text = artist.get_text().strip()
            if not text:
                continue
            text_strings.append(text)
            font_sizes.append(float(artist.get_fontsize()))

        self.assertIn("<", text_strings)
        self.assertIn(">", text_strings)
        self.assertGreaterEqual(min(font_sizes), REPORT.MIN_TEXT_SIZE)
        self.assertEqual(REPORT.COLORS["Fragmented"], "#E6AB02")
        self.assertEqual(REPORT.BLOCK_GLYPHS["b"], "<>")

    def test_flagged_scaffold_labels_are_right_aligned_without_category_suffix(self):
        blocks = {
            "Scaffold_1386.H1": [_block(1_200, 1_800, "p", 120_000)],
            "Longish_scaffold_alpha": [_block(1_400, 2_000, "p", 90_000)],
        }
        chrom_sizes = {
            "Scaffold_1386.H1": 120_000,
            "Longish_scaffold_alpha": 90_000,
        }
        classifications = OrderedDict(
            [
                ("Fragmented", ["Scaffold_1386.H1"]),
                ("Discordant", ["Longish_scaffold_alpha"]),
            ]
        )

        fig = REPORT.plot_overview_page1(classifications, blocks, chrom_sizes)
        self.addCleanup(REPORT.plt.close, fig)
        fig.canvas.draw()

        flagged_panel = next(ax for ax in fig.axes if ax.get_xlabel() == "Scaffold size (log10 bp)")
        labels = flagged_panel.texts
        self.assertTrue(any("\n" in text.get_text() for text in labels))

        renderer = fig.canvas.get_renderer()
        right_edges = [text.get_window_extent(renderer).x1 for text in labels]
        for text in labels:
            self.assertEqual(text.get_ha(), "right")
            self.assertNotIn("Fragmented", text.get_text())
            self.assertNotIn("Discordant", text.get_text())
        self.assertLess(max(right_edges) - min(right_edges), 1.0)

    def test_ranked_bar_panels_leave_left_spine_unclipped(self):
        page1_blocks = {"chrFlagged": [_block(1_200, 1_800, "p", 5_000)]}
        page1_sizes = {"chrFlagged": 5_000, "chrOther": 6_400}
        classifications = OrderedDict(
            [
                ("Fragmented", ["chrFlagged"]),
                ("Discordant", ["chrOther"]),
            ]
        )
        fig1 = REPORT.plot_overview_page1(classifications, page1_blocks, page1_sizes)
        self.addCleanup(REPORT.plt.close, fig1)
        page1_flagged = next(ax for ax in fig1.axes if ax.get_xlabel() == "Scaffold size (log10 bp)")
        self.assertIsNone(page1_flagged.spines["left"].get_bounds())
        self.assertTrue(page1_flagged.spines["bottom"].get_visible())

        page2_blocks = {"chrTel": [_block(1_400, 2_000, "p", 5_000)]}
        fig2 = REPORT.plot_overview_page2(page2_blocks, {"chrTel": 5_000})
        self.addCleanup(REPORT.plt.close, fig2)
        page2_flagged = next(ax for ax in fig2.axes if ax.get_xlabel() == "Flagged length (log10 bp)")
        self.assertIsNone(page2_flagged.spines["left"].get_bounds())
        self.assertTrue(page2_flagged.spines["bottom"].get_visible())

    def test_contig_rows_are_excluded_from_overview_statistics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "blocks.bed"
            bed_path.write_text(
                "chrMixed\t100\t700\t600\tp\tp\t3\t2\t5\t7\t20000\tscaffold\n"
                "chrMixed\t5000\t5400\t400\tq\tq\t2\t1\t3\t4\t20000\tcontig\n",
                encoding="utf-8",
            )

            parsed = REPORT.parse_terminal_bed(str(bed_path))

        self.assertEqual(len(parsed["chrMixed"]), 2)
        arm_blist = [b for b in parsed["chrMixed"] if b.get("term") == "scaffold"]
        contig_blist = [b for b in parsed["chrMixed"] if b.get("term") == "contig"]
        self.assertEqual(len(arm_blist), 1)
        self.assertEqual(len(contig_blist), 1)

        arm_blocks = {"chrMixed": arm_blist}
        contig_blocks = {"chrMixed": contig_blist}
        chrom_sizes = {"chrMixed": 20000}

        block_rows = REPORT._compute_block_rows(arm_blocks, chrom_sizes)
        total_arm_telomeres = len(block_rows)
        self.assertEqual(total_arm_telomeres, 1)

    def test_pair_fusions_pairs_a_q_then_p_within_distance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 1000, 1100, "q", "fusion", rev_can=5),
                _its_row("chrA", 1110, 1180, "p", "fragmentation", fwd_can=4),
            ], tmpdir)
            gaps = _gaps_frame([], tmpdir)

        pairs = REPORT.pair_fusions(df, gaps, 1000)

        self.assertEqual(len(pairs), 1)
        row = pairs.iloc[0]
        self.assertEqual((row["chr"], row["start"], row["end"]), ("chrA", 1000, 1180))
        self.assertEqual((row["q_bp"], row["p_bp"], row["min_arm"]), (100, 70, 70))
        self.assertEqual((row["combined_bp"], row["spacer_bp"]), (170, 10))
        # canonical_bp = matches x motif length (6); can_prop = canonical_bp / teloLen
        self.assertAlmostEqual(row["q_can_prop"], 30 / 100)
        self.assertAlmostEqual(row["p_can_prop"], 24 / 70)

    def test_pair_fusions_rejects_spacer_greater_than_d(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 1000, 1100, "q", "fusion", rev_can=5),
                _its_row("chrA", 3000, 3070, "p", "fragmentation", fwd_can=4),
            ], tmpdir)
            gaps = _gaps_frame([], tmpdir)

        pairs = REPORT.pair_fusions(df, gaps, 1000)

        self.assertEqual(len(pairs), 0)

    def test_pair_fusions_rejects_an_n_gap_between_the_pair(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 1000, 1100, "q", "fusion", rev_can=5),
                _its_row("chrA", 1110, 1180, "p", "fragmentation", fwd_can=4),
            ], tmpdir)
            gaps = _gaps_frame(["chrA\t1105\t1108\n"], tmpdir)

        pairs = REPORT.pair_fusions(df, gaps, 1000)

        self.assertEqual(len(pairs), 0)

    def test_pair_fusions_requires_the_c_engine_fusion_label(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 1000, 1100, "q", "single", rev_can=5),
                _its_row("chrA", 1110, 1180, "p", "fragmentation", fwd_can=4),
            ], tmpdir)
            gaps = _gaps_frame([], tmpdir)

        pairs = REPORT.pair_fusions(df, gaps, 1000)

        self.assertEqual(len(pairs), 0)

    def test_pair_fusions_ranks_by_min_arm_then_combined_bp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                # chrA: min_arm 50 (q 50 / p 200), combined 250
                _its_row("chrA", 0, 50, "q", "fusion"),
                _its_row("chrA", 60, 260, "p", "fragmentation"),
                # chrB: min_arm 90 (q 90 / p 90), combined 180
                _its_row("chrB", 0, 90, "q", "fusion"),
                _its_row("chrB", 100, 190, "p", "fragmentation"),
                # chrC: min_arm 90 (q 90 / p 95), combined 185 -- ties chrB on min_arm, wins on combined
                _its_row("chrC", 0, 90, "q", "fusion"),
                _its_row("chrC", 100, 195, "p", "fragmentation"),
            ], tmpdir)
            gaps = _gaps_frame([], tmpdir)

        pairs = REPORT.pair_fusions(df, gaps, 1000)

        self.assertEqual(list(pairs["chr"]), ["chrC", "chrB", "chrA"])
        self.assertEqual(list(pairs["min_arm"]), [90, 90, 50])

    def test_load_its_frame_computes_canonical_bp_and_can_prop(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 0, 100, "p", "single", fwd_can=4),   # engine floor: 4 matches
                _its_row("chrA", 200, 400, "q", "single", rev_can=10),
            ], tmpdir)

        floor_row, long_row = df.iloc[0], df.iloc[1]
        self.assertEqual(floor_row["canonical_bp"], 4 * 6)
        self.assertAlmostEqual(floor_row["can_prop"], 24 / 100)
        self.assertEqual(long_row["canonical_bp"], 10 * 6)
        self.assertAlmostEqual(long_row["can_prop"], 60 / 200)

    def test_load_its_frame_respects_a_non_default_motif_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "its.bed"
            bed_path.write_text(_its_row("chrA", 0, 100, "p", "single", fwd_can=4), encoding="utf-8")
            df = REPORT.load_its_frame(str(bed_path), motif_len=7)

        self.assertEqual(df.iloc[0]["canonical_bp"], 4 * 7)

    def test_rank_long_its_orders_by_canonical_bp_then_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 0, 500, "p", "single", fwd_can=10),  # canonical_bp=60, longest
                _its_row("chrB", 0, 100, "p", "single", fwd_can=10),  # canonical_bp=60, ties A, shorter
                _its_row("chrC", 0, 300, "p", "single", fwd_can=4),   # canonical_bp=24, lowest
            ], tmpdir)

        ranked = REPORT.rank_long_its(df, top_n=25)

        self.assertEqual(list(ranked["chr"]), ["chrA", "chrB", "chrC"])
        self.assertEqual(list(ranked["canonical_bp"]), [60, 60, 24])

    def test_compute_its_clusters_merges_within_the_gap_and_splits_beyond_it(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                _its_row("chrA", 0, 100, "p", "single"),
                _its_row("chrA", 40_000, 40_100, "p", "single"),      # 39_900 bp gap: same cluster
                _its_row("chrA", 200_000, 200_100, "p", "single"),    # 159_900 bp gap: new cluster
            ], tmpdir)

        clustered = REPORT.compute_its_clusters(df, merge_gap=50_000)

        self.assertEqual(list(clustered["cluster_id"]), [1, 1, 2])

    def test_summarize_its_clusters_filters_by_min_rows_and_ranks_by_its_bp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            df = _its_frame([
                # chrA: 3 rows within 50 kb, 300 bp total
                _its_row("chrA", 0, 100, "p", "single"),
                _its_row("chrA", 1_000, 1_100, "p", "single"),
                _its_row("chrA", 2_000, 2_100, "p", "single"),
                # chrB: only 2 rows -- below the min_rows=3 floor, dropped
                _its_row("chrB", 0, 5_000, "p", "single"),
                _its_row("chrB", 6_000, 11_000, "p", "single"),
                # chrC: 3 rows, 6_000 bp total -- ranked above chrA
                _its_row("chrC", 0, 2_000, "p", "single"),
                _its_row("chrC", 3_000, 5_000, "p", "single"),
                _its_row("chrC", 6_000, 8_000, "p", "single"),
            ], tmpdir)

        clusters = REPORT.summarize_its_clusters(df, merge_gap=50_000, min_rows=3)

        self.assertEqual(list(clusters["chr"]), ["chrC", "chrA"])
        self.assertEqual(list(clusters["rows"]), [3, 3])
        self.assertEqual(clusters.iloc[0]["its_bp"], 6_000)
        self.assertEqual(clusters.iloc[1]["its_bp"], 300)

    def test_pad_window_defaults_to_twice_the_span_floored_at_10kb(self):
        self.assertEqual(REPORT.pad_window(100_000, 100_200, 1_000_000), (90_000, 110_200))
        self.assertEqual(REPORT.pad_window(0, 100_000, 1_000_000), (0, 300_000))

    def test_read_params_parses_the_hash_params_line(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "report.tsv"
            report_path.write_text(
                "#teloscope version=0.1.6 commit=abc1234\n"
                "#params canonical=CCCTAA/TTAGGG patterns=38 window=1000 step=1000 "
                "terminal_limit=20000 max_match_dist=50 max_block_dist=500 min_block_len=300 "
                "min_block_density=0.5 min_block_counts=2 min_canonical_count=4 "
                "terminal_tolerance=3000 label_threshold=0.667 edit_distance=1 "
                "ultra_fast=false manual_curation=false\n"
                "#columns\tpos\theader\n",
                encoding="utf-8",
            )

            params = REPORT.read_params(str(report_path))

        self.assertEqual(params, {"max_block_dist": 500, "terminal_limit": 20000, "ultra_fast": False,
                                  "motif_len": 6, "min_canonical_count": 4,
                                  "known_params": {"max_block_dist", "terminal_limit", "motif_len", "min_canonical_count"}})

    def test_read_params_defaults_when_the_report_is_missing(self):
        defaults = {"max_block_dist": 1000, "terminal_limit": None, "ultra_fast": None,
                   "motif_len": 6, "min_canonical_count": 4, "known_params": set()}
        self.assertEqual(REPORT.read_params(None), defaults)
        self.assertEqual(REPORT.read_params("/no/such/report.tsv"), defaults)

    def test_read_params_keeps_other_keys_when_one_value_is_malformed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "report.tsv"
            report_path.write_text(
                "#params max_block_dist=500 terminal_limit=none ultra_fast=false\n",
                encoding="utf-8",
            )

            params = REPORT.read_params(str(report_path))

        self.assertEqual(params, {"max_block_dist": 500, "terminal_limit": None, "ultra_fast": False,
                                  "motif_len": 6, "min_canonical_count": 4, "known_params": {"max_block_dist"}})

    def test_read_params_derives_motif_length_from_the_canonical_pattern(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "report.tsv"
            report_path.write_text(
                "#params canonical=CCCTAA/TTAGGG min_canonical_count=4 ultra_fast=false\n",
                encoding="utf-8",
            )

            params = REPORT.read_params(str(report_path))

        self.assertEqual(params["motif_len"], 6)
        self.assertEqual(params["min_canonical_count"], 4)

    def test_load_gaps_frame_empty_file_keeps_int64_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gaps = _gaps_frame([], tmpdir)

        self.assertEqual(len(gaps), 0)
        self.assertEqual(str(gaps["start"].dtype), "int64")
        self.assertEqual(str(gaps["end"].dtype), "int64")


class ResilientITSReportTests(unittest.TestCase):
    def frame(self, rows):
        with tempfile.TemporaryDirectory() as tmpdir:
            return _its_frame(rows, tmpdir)

    def parts(self, df):
        pairs = REPORT.pair_fusions(df, REPORT.pd.DataFrame(columns=["chr", "start", "end"]), 1000)
        clusters = REPORT.summarize_its_clusters(df)
        sizes = {c: max(int(g["chrSize"].max()), int(g["end"].max()))
                 for c, g in df.groupby("chr")}
        df.attrs["display_labels"] = REPORT._its_labels(df, sizes)
        return pairs, clusters, sizes

    def draw(self, fig, name):
        self.addCleanup(REPORT.plt.close, fig)
        fig.canvas.draw()
        np.testing.assert_allclose(fig.get_size_inches(), [7.2, 3.7])
        renderer = fig.canvas.get_renderer()
        for text in fig.findobj(mtext.Text):
            if text.get_text().strip() and text.get_visible():
                self.assertGreaterEqual(text.get_fontsize(), 5, text.get_text())
        for text in fig.texts:
            box = text.get_window_extent(renderer)
            self.assertGreaterEqual(box.x0, -1, text.get_text())
            self.assertLessEqual(box.x1, fig.bbox.width + 1, text.get_text())
        for ax in fig.axes:
            for table in ax.tables:
                for cell in table.get_celld().values():
                    text_box = cell.get_text().get_window_extent(renderer)
                    self.assertLessEqual(text_box.width, cell.get_window_extent(renderer).width,
                                         cell.get_text().get_text())

    def test_bad_its_rows_are_rejected_individually(self):
        good = _its_row("chrA", 100, 200, "p", "single", fwd_can=4)
        rows = ["# comment\n", "track name=its\n", "browser position chrA\n", good,
                good.replace("\t100\t", "\tNaN\t", 1),
                _its_row("chrA", -1, 100, "q", "single"),
                _its_row("chrA", 200, 300, "p", "single", fwd_can=-1),
                "chrA\t1\t2\n",
                _its_row("chrA", 2**70, 2**70 + 100, "p", "single")]
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            df = self.frame(rows)
        self.assertEqual(len(df), 1)
        self.assertIn("Skipped 5 malformed", stderr.getvalue())

    def test_unknown_extent_is_not_a_q_end(self):
        df = self.frame([_its_row("chrA", 1000, 1100, "p", "single", chrom_size=0)])
        self.assertTrue(df["pos_frac"].isna().all())
        self.assertTrue(df["end_dist"].isna().all())
        pairs, clusters, sizes = self.parts(df)
        self.draw(REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes,
                                               REPORT.read_params(None)), "unknown_extent")

    def test_header_like_scaffold_names_and_conflicting_sizes(self):
        df = self.frame([_its_row(c, 100, 200, "p", "single") for c in
                         ("track_001", "browser_chr", "track", "browser")])
        self.assertEqual(len(df), 4)
        with contextlib.redirect_stderr(io.StringIO()):
            df = self.frame([_its_row("chrA", 100, 200, "p", "single", chrom_size=size)
                             for size in (1000, 2000)])
        self.assertTrue((df["chrSize"] == 0).all())
        self.assertTrue(df["pos_frac"].isna().all())

    def test_bedgraph_invalid_values_do_not_reach_plots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "track.bedgraph"
            path.write_text("track_001\t10\t20\t0.5\ntrack_001\t0\t10\t0.2\n"
                            "track_001\t20\t30\tnan\ntrack_001\t-1\t5\t0.9\n")
            with contextlib.redirect_stderr(io.StringIO()):
                data = REPORT.parse_bedgraph(str(path))
        np.testing.assert_array_equal(data["track_001"][0], [0, 10])
        np.testing.assert_allclose(data["track_001"][2], [0.2, 0.5])

    def test_populated_rank_tables_and_partial_scan_geometry(self):
        rows = []
        for i in range(8):
            chrom = f"very_long_sample_haplotype_identifier_scaffold_{i}"
            rows.extend([_its_row(chrom, 1000, 2000, "q", "fusion", rev_can=100),
                         _its_row(chrom, 2010, 4010, "p", "single", fwd_can=200),
                         _its_row(chrom, 8000, 8500, "b", "tail_to_tail", fwd_can=40)])
        df = self.frame(rows)
        pairs, clusters, sizes = self.parts(df)
        params = {"ultra_fast": True, "terminal_limit": 10_000, "motif_len": 6,
                  "min_canonical_count": 4, "max_block_dist": 1000,
                  "known_params": {"motif_len", "min_canonical_count", "max_block_dist"}}
        self.draw(REPORT.plot_its_composition_page(df, pairs, clusters, params), "full_candidates")
        self.draw(REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes, params), "partial_atlas")
        self.draw(REPORT.plot_its_statistics_page(df, pairs, clusters, params), "partial_statistics")

    def test_unknown_extent_is_disclosed_on_selected_locus(self):
        fig = REPORT.plot_its_loci_page([("L1", "chrA", 0, 1100)], {"chrA": 1100}, {},
                    {"chrA": [_block(1000, 1100, "p", 0)]}, {}, None, None, None, None, None,
                    unknown_extents={"chrA"})
        self.draw(fig, "unknown_locus")
        self.assertTrue(any("scaffold length unknown" in ax.get_xlabel() for ax in fig.axes))

    def test_sub_kilobase_region_ticks_are_distinguishable(self):
        for start, end in ((0, 1100), (100, 200), (0, 1), (100_000_000, 100_000_010)):
            ticks = REPORT._region_axis_spec(start, end)["tick_labels"]
            self.assertEqual(len(ticks), len(set(ticks)), (start, end, ticks))

    def test_failed_builder_closes_partial_figures_and_preserves_existing_ones(self):
        sentinel = REPORT.plt.figure()
        self.addCleanup(REPORT.plt.close, sentinel)
        saved = []
        def broken():
            REPORT.plt.figure()
            raise ValueError("synthetic rendering failure")
        with contextlib.redirect_stderr(io.StringIO()):
            ok, error = REPORT._save_figure_with_fallback(
                lambda fig: saved.append(fig.get_size_inches()), broken, "ITS", "Cannot render.")
        self.assertFalse(ok)
        self.assertIn("synthetic rendering failure", error)
        np.testing.assert_allclose(saved[0], [7.2, 3.7])
        self.assertEqual(REPORT.plt.get_fignums(), [sentinel.number])

    def test_gap_starting_before_spacer_rejects_fusion(self):
        df = self.frame([_its_row("chrA", 1000, 1100, "q", "fusion"),
                         _its_row("chrA", 1110, 1200, "p", "single")])
        gaps = REPORT.pd.DataFrame([["chrA", 1090, 1105]], columns=["chr", "start", "end"])
        self.assertTrue(REPORT.pair_fusions(df, gaps, 1000).empty)
        gaps["end"] = 1100  # half-open interval abuts spacer, does not overlap
        self.assertEqual(len(REPORT.pair_fusions(df, gaps, 1000)), 1)

    def test_overlapping_fusion_arrays_keep_engine_zero_distance(self):
        df = self.frame([_its_row("chrA", 1000, 1120, "q", "fusion"),
                         _its_row("chrA", 1110, 1200, "p", "single")])
        pairs = REPORT.pair_fusions(df, None, 1000)
        self.assertEqual(pairs.iloc[0]["spacer_bp"], 0)

    def test_ranks_do_not_depend_on_input_order(self):
        df = self.frame([_its_row(c, 100, 200, "p", "single", fwd_can=4)
                         for c in ("chrZ", "chrB", "chrA")])
        a = REPORT.rank_long_its(df)
        b = REPORT.rank_long_its(df.iloc[::-1])
        REPORT.pd.testing.assert_frame_equal(a, b)
        self.assertEqual(list(a["chr"]), ["chrA", "chrB", "chrZ"])

    def test_sparse_and_zero_canonical_pages_render_without_distortion(self):
        for label, rows in (
            ("empty", []),
            ("singleton", [_its_row("chrA", 100, 200, "b", "single", fwd_can=4)]),
            ("zero", [_its_row("chrA", 100, 200, "?", "new_class")]),
            ("identical", [_its_row("chrA", 100, 200, "p", "single", fwd_can=4)] * 800),
        ):
            with self.subTest(label=label):
                df = self.frame(rows)
                pairs, clusters, sizes = self.parts(df)
                params = REPORT.read_params(None)
                stats = REPORT.plot_its_statistics_page(df, pairs, clusters, params)
                self.draw(stats, label + "_statistics")
                self.draw(REPORT.plot_its_composition_page(df, pairs, clusters, params),
                          label + "_candidates")
                self.draw(REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes, params),
                          label + "_atlas")
                if label == "zero":
                    texts = [t.get_text() for t in stats.findobj(mtext.Text)]
                    self.assertIn("Positive n=0; zero canonical n=1", texts)
                    self.assertIn("unknown (n=1)", texts)
                    self.assertIn("other", texts)

    def test_low_canonical_observation_keeps_its_true_coordinate(self):
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=1)])
        fig, ax = REPORT.plt.subplots()
        self.addCleanup(REPORT.plt.close, fig)
        REPORT._draw_its_composition_panel(ax, df, 6, 4)
        self.assertEqual(float(ax.collections[0].get_offsets()[0, 1]), 6)

    def test_dense_population_and_all_scaffolds_are_accounted_for(self):
        n = 61_000
        i = np.arange(n)
        df = REPORT.pd.DataFrame({
            "chr": ["long_sample_haplotype_scaffold_" + str(j % 83) for j in i],
            "start": i * 1000, "end": i * 1000 + 100 + i % 300,
            "teloLen": 100 + i % 300, "teloLabel": np.array(["p", "q", "b"])[i % 3],
            "teloType": np.array(REPORT.CLASS_ORDER)[i % 4], "chrSize": 100_000_000,
            "canonical_bp": (i % 31) * 6, "can_prop": (i % 31) * 6 / (100 + i % 300),
            "pos_frac": i * 1000 / 100_000_000,
        })
        pairs = REPORT.pd.DataFrame(columns=REPORT._PAIR_COLUMNS)
        clusters = REPORT.summarize_its_clusters(df)
        sizes = {c: 100_000_000 for c in df["chr"].unique()}
        df.attrs["display_labels"] = REPORT._its_labels(df, sizes)
        counts = REPORT._its_position_counts(df, sizes)
        self.assertEqual(sum(v.sum() for _, v in counts.values()), n)
        summary = REPORT.its_scaffold_summary(df, {}, sizes)
        self.assertEqual(summary["rows"].sum(), n)
        self.assertEqual(summary["display_label"].nunique(), 83)
        chroms = REPORT._its_atlas_chroms(df, {}, sizes)
        seen = []
        for page_idx, start in enumerate(range(0, len(chroms), REPORT.ITS_ATLAS_ROWS)):
            chunk = chroms[start:start + REPORT.ITS_ATLAS_ROWS]
            seen.extend(chunk)
            fig = REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes,
                    {"ultra_fast": False}, chunk, page_idx + 1, 5, counts)
            self.draw(fig, f"dense_atlas_{page_idx + 1}")
        self.assertEqual(seen, chroms)
        stats = REPORT.plot_its_statistics_page(df, pairs, clusters, REPORT.read_params(None))
        self.draw(stats, "dense_statistics")
        hexagons = stats.axes[0].collections[0]
        self.assertEqual(int(hexagons.get_array().sum()), int((df["canonical_bp"] > 0).sum()))
        self.draw(REPORT.plot_its_composition_page(df, pairs, clusters, REPORT.read_params(None)),
                  "dense_candidates")

    def test_its_palette_and_pdf_geometry_match_terminal(self):
        for key in ("p", "q", "b"):
            self.assertEqual(REPORT.ITS_ORIENT_COLORS[key], REPORT.COLORS[key])
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=4)])
        pairs, clusters, _ = self.parts(df)
        fig = REPORT.plot_its_statistics_page(df, pairs, clusters, REPORT.read_params(None))
        self.addCleanup(REPORT.plt.close, fig)
        pdf = io.BytesIO()
        fig.savefig(pdf, format="pdf", dpi=450)
        self.assertIn(b"/MediaBox [ 0 0 518.4 266.4 ]", pdf.getvalue())
        self.assertIn(b"/FontFile2", pdf.getvalue())
        self.assertNotIn(b"/Subtype /Type3", pdf.getvalue())

    def test_single_locus_with_all_tracks_and_long_identifier(self):
        chrom = "sample#haplotype#extremely_long_scaffold_identifier"
        size = 1_000_000
        its = {chrom: [_block(400_000, 400_600, "b", size)]}
        track = {chrom: (np.array([390_000, 400_000]), np.array([400_000, 410_000]), np.array([0.2, 0.8]))}
        fig = REPORT.plot_its_loci_page([("L1", chrom, 390_000, 410_000)], {chrom: size},
                                        {}, its, {}, track, track, track, track, track)
        self.draw(fig, "locus_all_tracks")

    def test_its_only_cli_and_terminal_selection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "sample_interstitial_telomeres.bed").write_text(
                _its_row("chrA", 100, 200, "p", "single", fwd_can=4), encoding="utf-8")
            output = root / "its.pdf"
            stderr = io.StringIO()
            with mock.patch.object(sys, "argv", ["report", tmpdir, "--section", "its", "-o", str(output)]):
                with contextlib.redirect_stderr(stderr):
                    REPORT.main()
            self.assertTrue(output.is_file())
            self.assertNotIn("placeholder", stderr.getvalue())
            self.assertNotIn("Assembly overview", stderr.getvalue())
            exported = REPORT.pd.read_csv(root / "sample_its_rows.tsv", sep="\t")
            self.assertEqual(len(exported), 1)
            (root / "sample_terminal_telomeres.bed").write_text("", encoding="utf-8")
            with mock.patch.object(sys, "argv", ["report", tmpdir, "--section", "terminal",
                                                "-o", str(root / "terminal.pdf")]):
                with contextlib.redirect_stderr(io.StringIO()) as log:
                    REPORT.main()
            self.assertNotIn("ITS distributions", log.getvalue())
            self.assertEqual(REPORT.plt.get_fignums(), [])

    def test_default_split_combined_opt_in_and_empty_its_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "sample_terminal_telomeres.bed").write_text("", encoding="utf-8")
            (root / "sample_interstitial_telomeres.bed").write_text("", encoding="utf-8")
            base = root / "custom.pdf"
            with mock.patch.object(sys, "argv", ["report", tmpdir, "-o", str(base)]):
                with contextlib.redirect_stderr(io.StringIO()) as log:
                    REPORT.main()
            self.assertTrue((root / "custom_terminal.pdf").is_file())
            self.assertTrue((root / "custom_its.pdf").is_file())
            self.assertFalse(base.exists())
            self.assertNotIn("placeholder", log.getvalue())
            with mock.patch.object(sys, "argv", ["report", tmpdir, "--section", "all", "-o", str(base)]):
                with contextlib.redirect_stderr(io.StringIO()):
                    REPORT.main()
            self.assertTrue(base.is_file())
            (root / "sample_terminal_telomeres.bed").unlink()
            with mock.patch.object(sys, "argv", ["report", tmpdir, "--section", "all",
                                                "-o", str(root / "empty-its-only.pdf")]):
                with contextlib.redirect_stderr(io.StringIO()) as log:
                    REPORT.main()
            self.assertTrue((root / "empty-its-only.pdf").is_file())
            self.assertIn("ITS distributions", log.getvalue())
            self.assertNotIn("placeholder", log.getvalue())


if __name__ == "__main__":
    unittest.main()
