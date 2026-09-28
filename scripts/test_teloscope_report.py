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
        legend_axis = legends[0].axes
        legend_axis_center = (legend_axis.get_position().x0 + legend_axis.get_position().x1) / 2.0
        self.assertAlmostEqual(legend_group_center, legend_axis_center, delta=0.02)
        # Docked right under panel a's x-label rather than floating lower down.
        summary = next(ax for ax in fig.axes if ax.get_title() == "Scaffold classification")
        renderer = fig.canvas.get_renderer()
        label_bottom = summary.xaxis.label.get_window_extent(renderer).y0
        legend_top = max(legend.get_window_extent(renderer).y1 for legend in legends)
        self.assertLess(label_bottom - legend_top, 0.15 * fig.dpi)
        self.assertGreaterEqual(label_bottom - legend_top, 0)

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

        tile_ax = next(ax for ax in fig.axes if "Terminal telomeres" in [t.get_text() for t in ax.texts])
        texts = [t.get_text() for t in tile_ax.texts]
        tiles = dict(zip(texts[1::2], texts[0::2]))
        self.assertEqual(tiles["Scaffolds"], "4")
        self.assertEqual(tiles["Flagged scaffolds"], "3")
        self.assertEqual(tiles["Terminal telomeres"], "2")
        self.assertEqual(tiles["Distance-flagged telomeres"], "1")
        self.assertEqual(tiles["Median assembled array (bp)"], "600")

    def test_overview_tiles_lead_with_median_array_and_centre_on_the_page(self):
        blocks = {"chrA": [_block(0, 9_800, "p", 50_000)]}
        fig = REPORT.plot_overview_page1(OrderedDict([("Incomplete", ["chrA"])]), blocks, {"chrA": 50_000})
        self.addCleanup(REPORT.plt.close, fig)
        self.assertEqual(tuple(fig.get_size_inches()), REPORT.SLIDE_SIZE)
        tile_ax = next(ax for ax in fig.axes if "Terminal telomeres" in [t.get_text() for t in ax.texts])
        self.assertEqual([t.get_text() for t in tile_ax.texts][:10],
                         ["9,800", "Median assembled array (bp)", "1", "Scaffolds", "0", "Flagged scaffolds",
                          "1", "Terminal telomeres", "0", "Distance-flagged telomeres"])
        box = tile_ax.get_position()
        self.assertAlmostEqual(box.x0 + box.x1, 1.0, places=6)
        self.assertLessEqual(max(t.get_fontsize() for t in tile_ax.texts), 8.0)

    def test_ranked_bars_keep_one_thickness_and_pitch_for_any_count(self):
        heights, pitches = set(), set()
        for n in (1, 2, 10):
            fig = REPORT.plt.figure(figsize=(3.0, 1.6))
            self.addCleanup(REPORT.plt.close, fig)
            ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
            rows = [{"label": f"chr{i}", "value": 5.0, "color": "#888888"} for i in range(n)]
            REPORT._draw_ranked_bar_panel(ax, rows, "x", [0, 2, 4, 6])
            fig.canvas.draw()
            boxes = [p.get_window_extent() for p in ax.patches]
            heights.add(round(boxes[0].height, 3))
            if n > 1:
                pitches.add(round(boxes[0].y0 - boxes[1].y0, 3))
        self.assertEqual(len(heights), 1)
        self.assertEqual(len(pitches), 1)

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

        flagged_panel = next(ax for ax in fig.axes if ax.get_xlabel() == "Scaffold size (Mbp)")
        labels = [text for text in flagged_panel.texts if text.get_text() != "b"]
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
        page1_flagged = next(ax for ax in fig1.axes if ax.get_xlabel() == "Scaffold size (Mbp)")
        self.assertIsNone(page1_flagged.spines["left"].get_bounds())
        self.assertTrue(page1_flagged.spines["bottom"].get_visible())

        page2_blocks = {"chrTel": [_block(1_400, 2_000, "p", 5_000)]}
        fig2 = REPORT.plot_overview_page2(page2_blocks, {"chrTel": 5_000})
        self.addCleanup(REPORT.plt.close, fig2)
        page2_flagged = next(ax for ax in fig2.axes if ax.get_xlabel() == "Flagged length (kbp)")
        self.assertIsNone(page2_flagged.spines["left"].get_bounds())
        self.assertTrue(page2_flagged.spines["bottom"].get_visible())

    def test_positioning_panel_uses_plain_bp_decades_and_shares_length_scale(self):
        blocks = {"chrA": [_block(0, 6_000, "p", 90_000), _block(84_000, 90_000, "q", 90_000)],
                  "chrB": [_block(2_500, 14_000, "p", 90_000)]}
        fig = REPORT.plot_overview_page2(blocks, {"chrA": 90_000, "chrB": 90_000})
        self.addCleanup(REPORT.plt.close, fig)
        fig.canvas.draw()
        self.assertEqual(fig._suptitle.get_text(), "Terminal telomeres")
        self.assertEqual(tuple(fig.get_size_inches()), REPORT.SLIDE_SIZE)
        rain = next(ax for ax in fig.axes if ax.get_title() == "Length by arm")
        scatter = next(ax for ax in fig.axes if ax.get_title() == "Telomere positioning")
        self.assertEqual(scatter.get_xlabel(), "Distance to end (bp)")
        self.assertEqual(scatter.get_xscale(), "log")
        self.assertEqual([t.get_text() for t in scatter.get_xticklabels()], ["1", "10", "100", "1000", "10000"])
        self.assertEqual(rain.get_ylabel(), "Telomere length (kbp)")
        self.assertEqual(scatter.get_ylabel(), "Telomere length (kbp)")
        self.assertEqual(rain.get_ylim(), scatter.get_ylim())
        self.assertEqual(list(rain.get_yticks()), list(scatter.get_yticks()))

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
                                  "motif_len": 6, "min_canonical_count": 4, "label_threshold": 0.667,
                                  "known_params": {"max_block_dist", "terminal_limit", "motif_len",
                                                   "min_canonical_count", "label_threshold"}})

    def test_read_params_defaults_when_the_report_is_missing(self):
        defaults = {"max_block_dist": 1000, "terminal_limit": None, "ultra_fast": None,
                   "motif_len": 6, "min_canonical_count": 4, "label_threshold": 0.667, "known_params": set()}
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
                                  "motif_len": 6, "min_canonical_count": 4, "label_threshold": 0.667,
                                  "known_params": {"max_block_dist"}})

    def test_composition_class_matches_the_engine_thirds_rule(self):
        # fwdCan, revCan, fwdNonCan, revNonCan -> (strand, canonicity)
        cases = [((3, 0, 0, 0), (2, 2)), ((0, 3, 0, 0), (0, 2)), ((0, 0, 3, 0), (2, 0)),
                 ((0, 0, 0, 3), (0, 0)), ((1, 1, 1, 1), (1, 1)), ((2, 0, 0, 1), (1, 1)),
                 ((1, 0, 0, 2), (1, 1)), ((0, 0, 0, 0), (1, 1))]
        for counts, expected in cases:
            strand, canon = REPORT.composition_class(*[[v] for v in counts])
            self.assertEqual((int(strand[0]), int(canon[0])), expected, counts)

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

    def draw(self, fig, name, fixed_height=True):
        self.addCleanup(REPORT.plt.close, fig)
        fig.canvas.draw()
        if fixed_height:
            np.testing.assert_allclose(fig.get_size_inches(), [7.2, 3.7])
        else:
            self.assertEqual(fig.get_size_inches()[0], 7.2)
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
                                               REPORT.read_params(None)), "unknown_extent", False)

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
        self.draw(REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes, params), "partial_atlas", False)
        self.draw(REPORT.plot_its_summary_page(df, pairs, clusters, {}, params), "partial_statistics")

    def test_unknown_extent_is_disclosed_on_selected_locus(self):
        fig = REPORT.plot_its_loci_page([("L1", "chrA", 0, 1100)], {"chrA": 1100}, {},
                    {"chrA": [_block(1000, 1100, "p", 0)]}, {}, None, None, None, None, None,
                    unknown_extents={"chrA"})
        self.draw(fig, "unknown_locus")
        self.assertEqual(fig._suptitle.get_text(), "chrA  (≥1.1 kb)")
        self.assertTrue(any("Observed" in ax.get_ylabel() for ax in fig.axes))

    def test_sub_kilobase_region_ticks_are_distinguishable(self):
        for start, end in ((0, 1100), (100, 200), (0, 1), (100_000_000, 100_000_010)):
            for plain in (False, True):
                ticks = REPORT._region_axis_spec(start, end, plain)["tick_labels"]
                self.assertEqual(len(ticks), len(set(ticks)), (start, end, ticks))

    def test_locus_axis_and_title_share_one_unit(self):
        spec = REPORT._region_axis_spec(117_000, 171_000, plain=True)
        self.assertEqual(spec["xlabel"], "Position (kbp)")
        self.assertEqual(spec["tick_labels"], ["120", "130", "140", "150", "160", "170"])
        self.assertEqual(REPORT._fmt_region(117_000, 171_000), "117–171 kbp")
        self.assertEqual(REPORT._region_axis_spec(0, 20_000_000, plain=True)["xlabel"], "Position (Mbp)")
        self.assertEqual(REPORT._fmt_region(1_000_000, 11_000_000), "1–11 Mbp")
        self.assertEqual(REPORT._region_axis_spec(117_000, 171_000)["xlabel"], "Position (Mbp)")  # plot_its default

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
                stats = REPORT.plot_its_summary_page(df, pairs, clusters, {}, params)
                self.draw(stats, label + "_statistics")
                self.draw(REPORT.plot_its_summary_page(df, pairs, clusters, sizes, params), label + "_summary")
                self.draw(REPORT.plot_its_composition_page(df, pairs, clusters, params),
                          label + "_candidates")
                self.draw(REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes, params),
                          label + "_atlas", False)
                texts = [t.get_text() for t in stats.findobj(mtext.Text)]
                if label == "empty":
                    self.assertIn("No ITS", texts)
                    self.assertEqual(len(stats.axes), 0)

    def test_low_canonical_observation_keeps_its_true_coordinate(self):
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=1, rev_noncan=3)])
        pairs, clusters, _ = self.parts(df)
        fig = REPORT.plot_its_summary_page(df, pairs, clusters, {}, REPORT.read_params(None))
        self.draw(fig, "low_canonical")
        offsets = [c.get_offsets() for ax in fig.axes[3:4] for c in ax.collections if len(c.get_offsets())]
        np.testing.assert_allclose(offsets[0][0], [0.25, 0.25])

    def test_its_pages_never_say_rows(self):
        rows = [_its_row("chrA", 1000 + 1200 * i, 1500 + 1200 * i, "qp"[i % 2], "fusion",
                         fwd_can=10 * (i % 3), rev_noncan=5, fwd_noncan=i) for i in range(12)]
        df = self.frame(rows)
        pairs, clusters, _ = self.parts(df)
        params = REPORT.read_params(None)
        summary = lambda *a: REPORT.plot_its_summary_page(*a[:3], {"chrA": 100_000}, a[3])
        for build in (summary, REPORT.plot_its_composition_page):
            fig = build(df, pairs, clusters, params)
            self.draw(fig, "never_rows")
            for text in fig.findobj(mtext.Text):
                self.assertNotIn("row", text.get_text().lower(), text.get_text())
        self.assertFalse(pairs.empty)
        self.assertFalse(clusters.empty)

    def test_its_page_titles_are_single_bold_lines(self):
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=4)])
        pairs, clusters, _ = self.parts(df)
        params = REPORT.read_params(None)
        summary = lambda *a: REPORT.plot_its_summary_page(*a[:3], {}, a[3])
        for build, title in ((summary, "Assembly ITS summary"),
                             (REPORT.plot_its_composition_page, "ITS candidates")):
            fig = build(df, pairs, clusters, params)
            self.draw(fig, title)
            self.assertEqual(fig._suptitle.get_text(), title)
            self.assertEqual([t for t in fig.texts if t is not fig._suptitle], [])

    def test_statistics_page_scales_from_one_to_thousands(self):
        rng = np.random.default_rng(3)
        for n in (0, 1, 2000):
            with self.subTest(n=n):
                rows = [_its_row(f"chr{i % 7}", 1000 * i, 1000 * i + int(rng.integers(30, 9000)),
                                 "pqb"[i % 3], "single", fwd_can=int(rng.integers(0, 50)),
                                 rev_can=int(rng.integers(0, 50)), fwd_noncan=int(rng.integers(0, 50)),
                                 rev_noncan=int(rng.integers(0, 50)), chrom_size=10_000_000)
                        for i in range(n)]
                df = self.frame(rows)
                pairs, clusters, _ = self.parts(df)
                fig = REPORT.plot_its_summary_page(df, pairs, clusters, {}, REPORT.read_params(None))
                self.draw(fig, f"stats_{n}")
                if n:
                    points = sum(len(c.get_offsets()) for c in fig.axes[3].collections)
                    self.assertEqual(points, n)

    def test_candidates_page_renders_without_pairs_or_clusters(self):
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=4),
                         _its_row("chrB", 100, 900, "q", "single", rev_noncan=40)])
        pairs, clusters, _ = self.parts(df)
        self.assertTrue(pairs.empty and clusters.empty)
        fig = REPORT.plot_its_composition_page(df, pairs, clusters, REPORT.read_params(None))
        self.draw(fig, "candidates_sparse")
        self.assertEqual([t.get_text() for t in fig.findobj(mtext.Text)].count("None"), 2)

    def test_dense_population_and_all_scaffolds_are_accounted_for(self):
        n = 61_000
        i = np.arange(n)
        df = REPORT.pd.DataFrame({
            "chr": ["long_sample_haplotype_scaffold_" + str(j % 83) for j in i],
            "start": i * 1000, "end": i * 1000 + 100 + i % 300,
            "teloLen": 100 + i % 300, "teloLabel": np.array(["p", "q", "b"])[i % 3],
            "teloType": np.array(REPORT.CLASS_ORDER)[i % 4], "chrSize": 100_000_000,
            "fwdCan": i % 7, "revCan": i % 5, "fwdNonCan": i % 11, "revNonCan": i % 13,
            "canonical_bp": (i % 31) * 6, "can_prop": (i % 31) * 6 / (100 + i % 300),
            "pos_frac": i * 1000 / 100_000_000,
            "fwdCan": i % 31, "revCan": (i % 5) * (i % 2), "fwdNonCan": i % 7, "revNonCan": i % 11,
        })
        pairs = REPORT.pd.DataFrame(columns=REPORT._PAIR_COLUMNS)
        clusters = REPORT.summarize_its_clusters(df)
        sizes = {c: 100_000_000 for c in df["chr"].unique()}
        df.attrs["display_labels"] = REPORT._its_labels(df, sizes)
        cells = REPORT._its_atlas_cells(df)
        self.assertEqual(sum(len(starts) for starts, _, _ in cells.values()), n)
        summary = REPORT.its_scaffold_summary(df, {}, sizes)
        self.assertEqual(summary["rows"].sum(), n)
        self.assertEqual(summary["display_label"].nunique(), 83)
        chroms = REPORT._its_atlas_chroms(df, {}, sizes)
        seen = []
        pages = REPORT._paginate_atlas(chroms, sizes)
        for page_idx, chunk in enumerate(pages):
            self.assertGreaterEqual(REPORT._atlas_row_in(REPORT._atlas_panels(chunk, sizes)),
                                    REPORT.ATLAS_ROW_IN[0])
            seen.extend(chunk)
            fig = REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes,
                    {"ultra_fast": False}, chunk, page_idx + 1, len(pages), cells)
            self.draw(fig, f"dense_atlas_{page_idx + 1}", False)
            renderer = fig.canvas.get_renderer()
            for ax in fig.axes:
                boxes = [t.get_window_extent(renderer) for t in ax.get_yticklabels() if "\n" in t.get_text()]
                for i, box in enumerate(boxes):
                    self.assertFalse(any(box.overlaps(other) for other in boxes[:i]))

        self.assertEqual(sorted(seen), sorted(chroms))
        stats = REPORT.plot_its_summary_page(df, pairs, clusters, {}, REPORT.read_params(None))
        self.draw(stats, "dense_statistics")
        points = sum(len(c.get_offsets()) for c in stats.axes[3].collections)
        self.assertEqual(points, int((df[["fwdCan", "revCan", "fwdNonCan", "revNonCan"]].sum(axis=1) > 0).sum()))
        self.draw(REPORT.plot_its_composition_page(df, pairs, clusters, REPORT.read_params(None)),
                  "dense_candidates")

    def test_its_palette_and_pdf_geometry_match_terminal(self):
        for key in ("p", "q", "b"):
            self.assertEqual(REPORT.ITS_ORIENT_COLORS[key], REPORT.COLORS[key])
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=4)])
        pairs, clusters, _ = self.parts(df)
        fig = REPORT.plot_its_summary_page(df, pairs, clusters, {}, REPORT.read_params(None))
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

    def test_homolog_grouping_keeps_haplotypes_adjacent(self):
        names = ["chr1_mat", "chr33_mat", "chr2_pat", "chr33_pat", "chr1_pat", "hap1_chr5", "chr5_hap2",
                 "chr7_h1", "chr7_h2", "chr9.1", "chr9.2", "s#1#chr4", "s#2#chr4", "chrZ_PAT",
                 "ptg000001l.1", "ptg000001l.2"]
        sizes = {n: 1_000_000 - 1_000 * i for i, n in enumerate(names)}
        self.assertEqual(REPORT._homolog_key("chr33_mat"), REPORT._homolog_key("chr33_pat"))
        self.assertEqual(REPORT._homolog_key("chrZ_PAT"), "chrZ")
        groups = REPORT._homolog_groups(names, sizes)
        self.assertIn(["chr33_mat", "chr33_pat"], groups)
        for pair in (["hap1_chr5", "chr5_hap2"], ["chr7_h1", "chr7_h2"], ["s#1#chr4", "s#2#chr4"]):
            self.assertIn(sorted(pair), groups)
        # Version or piece suffixes are not haplotypes.
        for name in ("chr9.1", "chr9.2", "ptg000001l.1", "ptg000001l.2"):
            self.assertIn([name], groups)
        order = [c for g in groups for c in g]
        self.assertEqual(abs(order.index("chr33_mat") - order.index("chr33_pat")), 1)
        many = [f"chr{i}_{h}" for i in range(1, 16) for h in ("mat", "pat")]
        pages = REPORT._paginate_atlas(many, {c: 100 - int(c[3:].split("_")[0]) for c in many})
        self.assertGreater(len(pages), 1)
        for page in pages:
            self.assertEqual(len({REPORT._homolog_key(c) for c in page}) * 2, len(page))

    def test_candidates_colours_follow_composition_class(self):
        df = self.frame([_its_row("chrA", 1000, 1600, "q", "fusion", rev_can=90),
                         _its_row("chrA", 1700, 2300, "p", "fusion", fwd_can=90),
                         _its_row("chrB", 500, 900, "p", "single")])
        pairs, clusters, _ = self.parts(df)
        self.assertEqual(len(pairs), 1)
        fig = REPORT.plot_its_composition_page(df, pairs, clusters, REPORT.read_params(None))
        self.draw(fig, "candidates_colours")
        self.assertFalse(fig.legends)  # the segment legend lives on its own axes under the title
        texts = [t.get_text() for t in fig.findobj(mtext.Text)]
        # Zero-count ITS takes the (both, mixed) cell, as on the composition page.
        long_ax = fig.axes[0]
        faces = {REPORT.matplotlib.colors.to_hex(p.get_facecolor()) for p in long_ax.patches}
        self.assertIn(REPORT.DOUBLE_KEY_COLORS[1][1].lower(), faces)
        self.assertIn("NA", texts)
        fusion = {REPORT.matplotlib.colors.to_hex(p.get_facecolor()) for p in fig.axes[2].patches}
        self.assertEqual(fusion, {REPORT.DOUBLE_KEY_COLORS[2][0].lower(), REPORT.DOUBLE_KEY_COLORS[2][2].lower()})
        # A pair whose arms match no ITS falls back to the neutral grey, never a canonical colour.
        fig, ax = REPORT.plt.subplots()
        self.addCleanup(REPORT.plt.close, fig)
        REPORT._draw_its_fusion_glyphs(ax, df[df["chr"] == "chrB"], pairs, REPORT.LABEL_THRESHOLD)
        faces = {REPORT.matplotlib.colors.to_hex(p.get_facecolor()) for p in ax.patches}
        self.assertEqual(faces, {REPORT.COLORS["its"].lower()})

    def test_long_its_canonical_share_is_count_based(self):
        df = self.frame([_its_row("chrA", 100, 160, "p", "single", fwd_can=30, fwd_noncan=10)])
        df["canonical_bp"] = 180  # overlapping matches push the bp ratio past one
        fig, ax = REPORT.plt.subplots()
        self.addCleanup(REPORT.plt.close, fig)
        REPORT._draw_its_long_glyphs(ax, df)
        self.assertIn("75% canonical", [t.get_text() for t in ax.texts])

    def test_long_its_segments_run_fwd_before_rev(self):
        df = self.frame([_its_row("chrA", 100, 500, "b", "single", fwd_can=10, rev_can=20,
                                  fwd_noncan=30, rev_noncan=40)])
        fig, ax = REPORT.plt.subplots()
        self.addCleanup(REPORT.plt.close, fig)
        REPORT._draw_its_long_glyphs(ax, df)
        order = [REPORT.matplotlib.colors.to_hex(p.get_facecolor())
                 for p in sorted(ax.patches, key=lambda p: p.get_x())]
        key = REPORT.DOUBLE_KEY_COLORS
        self.assertEqual(order, [key[2][2], key[0][2], key[0][0], key[2][0]])
        np.testing.assert_allclose([p.get_width() for p in sorted(ax.patches, key=lambda p: p.get_x())],
                                   [40, 120, 160, 80])

    def test_candidates_legend_explains_segments_and_ticks_carry_no_units(self):
        df = self.frame([_its_row("chrA", 1000 + 3000 * i, 3500 + 3000 * i, "qp"[i % 2], "fusion",
                                  fwd_can=40 * (i % 2), rev_can=40 * (1 - i % 2), fwd_noncan=5)
                         for i in range(6)])
        pairs, clusters, _ = self.parts(df)
        self.assertFalse(pairs.empty or clusters.empty)
        fig = REPORT.plot_its_composition_page(df, pairs, clusters, REPORT.read_params(None))
        self.draw(fig, "candidates_legend")
        legends = [ax.get_legend() for ax in fig.axes if ax.get_legend()]
        self.assertEqual(len(legends), 1)
        labels = [t.get_text() for t in legends[0].get_texts()]
        for label in ("fwd canonical", "fwd non-canonical", "rev non-canonical", "rev canonical"):
            self.assertIn(label, labels)
        self.assertLess(labels.index("fwd canonical"), labels.index("rev canonical"))
        ticks = [t.get_text() for ax in fig.axes for t in ax.get_xticklabels() if t.get_text()]
        self.assertTrue(ticks)
        for text in ticks:
            self.assertRegex(text, r"^[0-9.,]+$")
        for ax in fig.axes[:3]:
            self.assertRegex(ax.get_xlabel(), r"\((bp|kbp|Mbp)\)")

    def test_panel_letters_share_the_terminal_left_margin(self):
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=4)])
        pairs, clusters, _ = self.parts(df)
        params = REPORT.read_params(None)
        figs = [REPORT.plot_its_summary_page(df, pairs, clusters, {}, params),
                REPORT.plot_its_summary_page(df, pairs, clusters, {}, params),
                REPORT.plot_its_composition_page(df, pairs, clusters, params),
                REPORT.plot_overview_page1(OrderedDict([("T2T", ["chrA"])]), {}, {"chrA": 1000}),
                REPORT.plot_overview_page2({}, {"chrA": 1000})]
        for fig in figs:
            self.addCleanup(REPORT.plt.close, fig)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            letters = [t for ax in fig.axes for t in ax.texts
                       if len(t.get_text()) == 1 and t.get_fontweight() == "bold"]
            self.assertTrue(letters)
            x0 = min(t.get_window_extent(renderer).x0 for t in letters) / fig.dpi
            self.assertAlmostEqual(x0, REPORT.PANEL_LETTER_X_IN, delta=0.03)

    def test_summary_page_tiles_and_empty_state(self):
        df = self.frame([_its_row("chr1_mat", 100, 400, "p", "single", fwd_can=40, chrom_size=5_000_000),
                         _its_row("chr1_pat", 100, 200, "q", "single", rev_can=15, chrom_size=4_000_000)])
        pairs, clusters, sizes = self.parts(df)
        fig = REPORT.plot_its_summary_page(df, pairs, clusters, sizes, REPORT.read_params(None))
        self.draw(fig, "summary_tiles")
        hist, per = fig.axes[1:3]
        self.assertEqual(hist.get_title(), "ITS length distribution")
        self.assertEqual(hist.get_ylabel(), "# ITSs")
        self.assertAlmostEqual(hist.bbox.width, hist.bbox.height)
        self.assertEqual(per.get_title(), "ITS scaffold distribution")
        self.assertEqual(per.get_ylabel(), "ITS density (bp/Mbp)")
        tiles = [t.get_text() for t in fig.axes[0].texts]
        self.assertEqual(tiles[1::2], ["Median ITS length (bp)", "ITS", "ITS content (kbp)",
                                       "Scaffolds with ITS", "Clusters", "Candidate fusions"])
        self.assertEqual(tiles[0::2], ["200", "2", "0.4", "2", "0", "0"])
        # One dot per scaffold, joined as homologs, on plain-number log axes.
        per = fig.axes[2]
        self.assertEqual(len(per.collections[0].get_offsets()), 2)
        self.assertEqual(len(per.lines), 1)
        for label in per.get_xticklabels() + per.get_yticklabels() + fig.axes[1].get_xticklabels():
            self.assertRegex(label.get_text(), r"^[0-9.]+$")
        empty = self.frame([])
        fig = REPORT.plot_its_summary_page(empty, *self.parts(empty)[:2], {}, REPORT.read_params(None))
        self.draw(fig, "summary_empty")
        self.assertEqual(len(fig.axes), 0)
        self.assertIn("No ITS", [t.get_text() for t in fig.texts])

    def test_summary_repel_and_composition_totals(self):
        df = self.frame([_its_row(f"chr{i}", 100, 400, "p", "single", fwd_can=40,
                                 chrom_size=5_000_000) for i in range(3)])
        pairs, clusters, sizes = self.parts(df)
        fig = REPORT.plot_its_summary_page(df, pairs, clusters, sizes, REPORT.read_params(None))
        self.draw(fig, "summary_coincident_labels")
        renderer = fig.canvas.get_renderer()
        per, joint = fig.axes[2:4]
        boxes = [mtext.Text.get_window_extent(t, renderer) for t in per.texts
                 if isinstance(t, mtext.Annotation)]
        self.assertEqual(len(boxes), 3)
        for i, box in enumerate(boxes):
            self.assertTrue(per.bbox.contains(box.x0, box.y0))
            self.assertTrue(per.bbox.contains(box.x1, box.y1))
            self.assertFalse(any(box.overlaps(other) for other in boxes[:i]))
        cells = [t.get_text() for t in joint.texts if "%" in t.get_text() and "Class" not in t.get_text()]
        self.assertEqual(len(cells), 9)
        self.assertIn("3\n100%", cells)
        self.assertEqual(sum(int(t.split("\n")[0].replace(",", "")) for t in cells), len(df))
        for text in joint.texts:
            self.assertIsNone(text.get_bbox_patch())
            self.assertGreaterEqual(text.get_window_extent(renderer).y0, 0)

    def test_forward_share_axis_is_mirrored(self):
        df = self.frame([_its_row("chrA", 100, 200, "p", "single", fwd_can=40),
                         _its_row("chrA", 900, 1000, "q", "single", rev_noncan=40)])
        pairs, clusters, _ = self.parts(df)
        fig = REPORT.plot_its_summary_page(df, pairs, clusters, {}, REPORT.read_params(None))
        self.draw(fig, "mirrored")
        joint, top = fig.axes[3], fig.axes[4]
        for ax in (joint, top):
            left, right = ax.get_xlim()
            self.assertGreater(left, right)
        self.assertEqual(joint.get_xlabel(), "← Forward share")
        # The fwd ITS sits left of the rev ITS in display space.
        xs = joint.transData.transform(joint.collections[0].get_offsets())[:, 0]
        colours = [REPORT.matplotlib.colors.to_hex(c) for c in joint.collections[0].get_facecolors()]
        self.assertLess(xs[colours.index(REPORT.DOUBLE_KEY_COLORS[2][2])],
                        xs[colours.index(REPORT.DOUBLE_KEY_COLORS[0][0])])

    def test_size_key_centres_its_references(self):
        for lengths in (np.array([50.0]), np.array([50.0, 900.0]), np.array([50.0, 9000.0])):
            fig, ax = REPORT.plt.subplots()
            self.addCleanup(REPORT.plt.close, fig)
            REPORT._draw_its_size_key(ax, lengths, False)
            xs = ax.collections[0].get_offsets()[:, 0]
            self.assertAlmostEqual(float(np.mean(xs)), 0.0)
            self.assertAlmostEqual(sum(ax.get_xlim()), 0.0)

    def test_top_five_loci_match_atlas_ranks(self):
        for n in (0, 1, 3, 5, 8):
            rows = []
            for i in range(n):
                rows.extend([_its_row(f"chr{i}", 1000, 2000, "q", "fusion", rev_can=100),
                             _its_row(f"chr{i}", 2010, 4010, "p", "single", fwd_can=200),
                             _its_row(f"chr{i}", 8000, 8500, "b", "tail_to_tail", fwd_can=40)])
            df = self.frame(rows)
            pairs, clusters, sizes = self.parts(df)
            loci = REPORT.resolve_its_loci(clusters, df, pairs, sizes)
            marks = REPORT._its_top_hits(df, pairs, clusters)
            expected = {f"{prefix}{rank}" for prefix, count in
                        (("L", len(df)), ("C", len(clusters)), ("F", len(pairs)))
                        for rank in range(1, min(count, 5) + 1)}
            self.assertEqual({tag for tag, _, _, _ in loci}, expected)
            self.assertEqual({tag for entries in marks.values() for tag, _ in entries}, expected)
            for tag, chrom, start, end in loci:
                midpoint = next(mid for label, mid in marks[chrom] if label == tag)
                self.assertLessEqual(start, midpoint)
                self.assertGreaterEqual(end, midpoint)

    def test_atlas_labels_keep_one_precision_and_mark_loci(self):
        df = self.frame([_its_row("chr1_mat", 1000, 1500, "p", "single", fwd_can=40, chrom_size=116_000_000),
                         _its_row("chr1_pat", 1000, 1500, "p", "single", fwd_can=40, chrom_size=115_910_000)])
        pairs, clusters, sizes = self.parts(df)
        fig = REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes, REPORT.read_params(None))
        self.draw(fig, "atlas_precision", False)
        ticks = [t.get_text() for t in fig.axes[0].get_yticklabels()]
        self.assertEqual([t.split("\n")[-1] for t in ticks], ["116.00 Mb", "115.91 Mb"])
        self.assertEqual([t.get_text() for t in fig.axes[0].get_xticklabels()][:3], ["0", "20", "40"])
        self.assertEqual(fig.axes[0].get_xlabel(), "Position (Mbp)")
        labels = [t.get_text() for leg in fig.legends for t in leg.get_texts()]
        self.assertIn("Scaffold", labels)
        self.assertIn("(L)  Longest canonical ITS", labels)
        texts = [t.get_text() for t in fig.findobj(mtext.Text)]
        self.assertFalse([t for t in texts if "top long ITS" in t or "—" in t], texts)
        self.assertIn("L1", texts)
        triangles = [l for l in fig.axes[0].get_lines() if l.get_marker() == "v"]
        self.assertEqual(len(triangles), 2)
        self.assertEqual(float(triangles[0].get_xdata()[0]), 1250.0)  # at the locus, not the scaffold end

    def test_atlas_page_is_one_slide_and_centred(self):
        def atlas(n):
            rows = [_its_row(f"chr{i}_{h}", 1000, 1500, "p", "single", fwd_can=40, chrom_size=5_000_000)
                    for i in range(n // 2) for h in ("mat", "pat")]
            df = self.frame(rows)
            pairs, clusters, sizes = self.parts(df)
            fig = REPORT.plot_its_overview_page(df, pairs, clusters, {}, sizes, REPORT.read_params(None))
            self.draw(fig, f"atlas_{n}", False)
            return fig
        two, twenty = atlas(2), atlas(20)
        pitch = []
        for fig in (two, twenty):
            self.assertEqual(tuple(fig.get_size_inches()), REPORT.SLIDE_SIZE)
            texts = [t.get_text() for t in fig.findobj(mtext.Text)]
            self.assertFalse([t for t in texts if "row" in t.lower()], texts)
            self.assertEqual(fig._suptitle.get_text(), "ITS atlas")
            self.assertEqual([t.get_text() for t in fig.texts], ["ITS atlas"])  # no footer
            box = fig.get_tightbbox(fig.canvas.get_renderer())
            self.assertAlmostEqual(box.x0, fig.get_size_inches()[0] - box.x1, delta=0.05)
            ax = fig.axes[0]
            pitch.append(ax.get_position().height * 3.7 / (ax.get_ylim()[0] - ax.get_ylim()[1]))
        self.assertAlmostEqual(pitch[0], REPORT.ATLAS_ROW_IN[1])  # few scaffolds: thick but capped
        self.assertGreaterEqual(pitch[1], REPORT.ATLAS_ROW_IN[0])
        self.assertLess(pitch[1], pitch[0])

    def test_locus_page_marks_the_zoom_with_a_box_only(self):
        chrom, size = "chr33_mat", 1_000_000
        its = {chrom: [dict(_block(400_000, 400_600, "p", size), fwdCan=0, revCan=0, fwdNonCan=90, revNonCan=0)]}
        fig = REPORT.plot_its_loci_page([("C1", chrom, 390_000, 410_000)], {chrom: size}, {}, its, {},
                                        None, None, None, None, None)
        self.draw(fig, "locus_box")
        self.assertFalse(fig.findobj(REPORT.ConnectionPatch))
        self.assertEqual(fig._suptitle.get_text(), "chr33_mat  (1 Mb)")
        texts = [t.get_text() for t in fig.findobj(mtext.Text)]
        self.assertIn("C1  390–410 kbp", texts)
        self.assertFalse([t for t in texts if "row" in t.lower()], texts)
        faces = {REPORT.matplotlib.colors.to_hex(c) for ax in fig.axes for coll in ax.collections
                 if isinstance(coll, REPORT.PatchCollection) for c in coll.get_facecolors()}
        self.assertIn(REPORT.DOUBLE_KEY_COLORS[0][2].lower(), faces)

    def test_locus_entropy_and_key_are_legible(self):
        chrom, size = "chr33_mat", 1_000_000
        its = {chrom: [_block(400_000, 400_600, "p", size)]}
        track = {chrom: (np.array([390_000, 400_000]), np.array([400_000, 410_000]), np.array([1.2, 1.9]))}
        fig = REPORT.plot_its_loci_page([("C1", chrom, 390_000, 410_000)], {chrom: size}, {}, its, {},
                                        None, None, None, None, track)
        self.draw(fig, "locus_entropy")
        ranges = [tuple(float(t) for t in ax.get_yticks()) for ax in fig.axes if ax.get_visible()]
        self.assertIn((1.0, 1.5, 2.0), ranges)  # the full 1-2 bit range, so dips are not clipped
        renderer = fig.canvas.get_renderer()
        key = next(ax for ax in fig.axes if ax.get_xlabel() == "← Forward share")
        self.assertAlmostEqual(key.get_position().width * 7.2, REPORT.LOCUS_KEY_SIDE, places=3)
        self.assertLessEqual(REPORT.LOCUS_KEY_SIDE, 0.45)
        self.assertLessEqual(key.get_tightbbox(renderer).x1, 0.97 * fig.bbox.width + 0.5)
        self.assertEqual({t.get_fontsize() for t in key.get_xticklabels()}, {REPORT.MIN_TEXT_SIZE})
        self.assertIn("Position (kbp)", [ax.get_xlabel() for ax in fig.axes])
        self.assertEqual(tuple(fig.get_size_inches()), REPORT.SLIDE_SIZE)
        ticks = [t.get_window_extent(renderer) for t in key.get_xticklabels()]
        for left, right in zip(ticks, ticks[1:]):
            self.assertLess(left.x1, right.x0)
        key_box = key.get_tightbbox(renderer)
        self.assertLessEqual(key_box.y1, fig.bbox.height)
        tracks_top = max(ax.get_position().y1 for ax in fig.axes if ax is not key) * fig.bbox.height
        self.assertGreater(key_box.y0, tracks_top)
        for text in (fig._suptitle, *[ax.title for ax in fig.axes if ax.get_title()]):
            self.assertFalse(key_box.overlaps(text.get_window_extent(renderer)), text.get_text())

    def test_shared_locus_tracks_keep_orientation_colours_by_default(self):
        size = 1_000_000
        its = [dict(_block(400_000, 400_600, k, size), fwdCan=0, revCan=0, fwdNonCan=90, revNonCan=0)
               for k in ("p", "q", "b")]
        fig, (ax_ideo, ax_blocks) = REPORT.plt.subplots(2, 1)
        self.addCleanup(REPORT.plt.close, fig)
        REPORT._draw_locus_ideogram(ax_ideo, size, 390_000, 410_000, [], its)
        REPORT._draw_its_blocks_track(ax_blocks, 390_000, 410_000, [], its, [])
        lines = {REPORT.matplotlib.colors.to_hex(c) for coll in ax_ideo.collections
                 for c in coll.get_colors()}
        faces = [REPORT.matplotlib.colors.to_hex(c) for coll in ax_blocks.collections
                 for c in coll.get_facecolors()]
        expected = {REPORT.ITS_ORIENT_COLORS[k] for k in ("p", "q", "b")}
        self.assertEqual(lines, {c.lower() for c in expected})
        self.assertEqual(set(faces), {c.lower() for c in expected})

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
            with mock.patch.object(sys, "argv", ["report", tmpdir, "--section", "its", "--png",
                                                "-o", str(root / "png")]):
                with contextlib.redirect_stderr(io.StringIO()) as log:
                    REPORT.main()
            self.assertTrue((root / "png" / "teloscope_its_summary.png").is_file())
            self.assertFalse((root / "png" / "teloscope_its_statistics.png").exists())
            self.assertNotIn("placeholder", log.getvalue())
            self.assertNotIn("Assembly overview", stderr.getvalue())
            exported = REPORT.pd.read_csv(root / "sample_its_rows.tsv", sep="\t")
            self.assertEqual(len(exported), 1)
            (root / "sample_terminal_telomeres.bed").write_text("", encoding="utf-8")
            with mock.patch.object(sys, "argv", ["report", tmpdir, "--section", "terminal",
                                                "-o", str(root / "terminal.pdf")]):
                with contextlib.redirect_stderr(io.StringIO()) as log:
                    REPORT.main()
            self.assertNotIn("ITS composition", log.getvalue())
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
            self.assertIn("ITS summary", log.getvalue())
            self.assertNotIn("placeholder", log.getvalue())


if __name__ == "__main__":
    unittest.main()
