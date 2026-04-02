import atexit
import importlib.util
import os
from collections import OrderedDict
from pathlib import Path
import shutil
import tempfile
import unittest

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


def _block(start, end, label, chrom_size):
    return {
        "start": int(start),
        "end": int(end),
        "length": int(end - start),
        "label": label,
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
    def test_parse_report_maps_expected_categories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "synthetic_report.tsv"
            report_path.write_text(
                "+++\n"
                "pos\theader\ttype\n"
                "1\tchr_dis\tdiscordant\n"
                "2\tchr_mis\tgapped_missassembly\n"
                "3\tchr_none\tnone\n",
                encoding="utf-8",
            )

            parsed = REPORT.parse_report(str(report_path))

        self.assertEqual(parsed["Discordant"], ["chr_dis"])
        self.assertEqual(parsed["Gapped Misassembly"], ["chr_mis"])
        self.assertEqual(parsed["No telomeres"], ["chr_none"])

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
                ("Misassembly", ["chrC"]),
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
                ("Misassembly", ["chrNear", "chrMissing"]),
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
        self.assertEqual(REPORT.COLORS["Misassembly"], "#8278F4")
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
                ("Misassembly", ["Scaffold_1386.H1"]),
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
            self.assertNotIn("Misassembly", text.get_text())
            self.assertNotIn("Discordant", text.get_text())
        self.assertLess(max(right_edges) - min(right_edges), 1.0)

    def test_ranked_bar_panels_leave_left_spine_unclipped(self):
        page1_blocks = {"chrFlagged": [_block(1_200, 1_800, "p", 5_000)]}
        page1_sizes = {"chrFlagged": 5_000, "chrOther": 6_400}
        classifications = OrderedDict(
            [
                ("Misassembly", ["chrFlagged"]),
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


if __name__ == "__main__":
    unittest.main()
