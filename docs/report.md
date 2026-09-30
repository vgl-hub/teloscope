[Back to README](index.md)

# Report generation

Teloscope can write terminal and interstitial (ITS) PDF reports during a FASTA run:

```sh
teloscope asm.fa -o results/ -r -e -g -i --plot-report
```

This writes `results/asm.fa_plot_report_terminal.pdf` and `results/asm.fa_plot_report_its.pdf`. In the default fast scan the ITS report covers only the scanned ends; add `-i` for the whole genome.

Use `-r` with `--plot-report` to get the repeat-density, canonical-ratio, and strand-ratio tracks. GC and entropy tracks are added when their files exist.

## What the report contains

Terminal report:

- page 1: assembly telomere summary, scaffold classes, and flagged scaffolds
- page 2: terminal telomere length and position by arm, and flagged telomeres
- next pages: one zoom per scaffold with its telomere blocks

Contig-terminal rows (`-n`) are drawn outline-only and left out of the counts.

ITS report, colored by a 3x3 key where hue is the strand (fwd blue, both purple, rev vermillion) and opacity is the canonical share:

- **Summary:** headline numbers, ITS length distribution, ITS count vs. density per scaffold (dot size is scaffold size), and ITS composition by strand and canonical share.
- **Atlas:** every scaffold with a terminal telomere or ITS, drawn to scale, homologs side by side. Up to one slide in two size groups, else a double slide in four, each group on its own axis. Rows spread to fill the page. ▼ marks the top (F) fusion candidates, (L) longest canonical ITS and (C) clusters of ITS. If rows do not fit even a double slide, the atlas becomes an unlabelled heatmap.
- **Candidates:** the top ITS clusters (3 or more ITS chained within 50 kbp, ranked by summed length), candidate fusions around their junction (ranked by the shorter array), and the top ITS by canonical bp, the bases in exact canonical repeats.
- **Loci:** one track page per top cluster, fusion, and ITS.

A candidate fusion is a `q` ITS followed by a `p` ITS within `-d` on the same scaffold, with no gap between them and at least one of the two classed `fusion` in the BED. It marks a possible end-to-end join, not a confirmed one.

## Standalone plotting

The plotting script also runs on an existing output directory, for runs made without `--plot-report`:

```sh
python3 scripts/teloscope_report.py results/
python3 scripts/teloscope_report.py results/ --png -o figures/
```

Every run in the directory that has no report yet gets `<input>_plot_report_terminal.pdf` and `<input>_plot_report_its.pdf` beside its files, as with `--plot-report`, several at once. Name one run by its file stem, `results/asm.fa`, to plot it again. Reads-mode runs are left out.

To plot a single ITS locus:

```sh
python3 scripts/plot_its.py results/ CHROM:START-END -o its.pdf
python3 scripts/plot_its.py results/ CHROM --png -o its_figures/
```

A bare `CHROM` centers on its largest ITS cluster. When the directory holds several runs, name one by its file stem.

Required input is `*_terminal_telomeres.bed` for the terminal report or `*_interstitial_telomeres.bed` for the ITS report. Gaps, window bedgraphs, and `*_report.tsv` are used when present.

## Standalone script options

| Flag | Meaning | Default |
| --- | --- | --- |
| `-o` | directory for the PDFs or PNGs | the output directory |
| `-j` | runs of a directory plotted at once, one process each | all cores |
| `--png` | write one PNG per page instead of the PDFs | `false` |
| `--dpi` | raster DPI | `450` |
| `--draft` | use 150 DPI for fast iteration | `false` |

## Requirements

Python 3 with `matplotlib` 3.5 or newer, `numpy`, and `pandas`. `--plot-report` looks for `scripts/teloscope_report.py` in the source checkout, or for `teloscope_report.py` next to the binary.
