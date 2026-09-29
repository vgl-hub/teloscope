[Back to README](index.md)

# Report generation

Teloscope can write terminal and interstitial (ITS) PDF reports during a FASTA run:

```sh
teloscope asm.fa -o results/ -r -e -g -i --plot-report
```

This writes `results/asm.fa_plot_report_terminal.pdf` and `results/asm.fa_plot_report_its.pdf`. The ITS report needs an ITS BED.

Use `-r` with `--plot-report` to get the repeat-density, canonical-ratio, and strand-bias tracks. GC and entropy tracks are added when their files exist.

Record filters apply before anything is written. For example, to keep only chromosomes listed by accession.version:

```sh
teloscope asm.fa -o results/ --include-bed chromosomes.ids -r --plot-report
```

`--include-prefix` also works, but check the FASTA IDs first, since prefixes differ between databases.

## What the report contains

Terminal report:

- page 1: assembly telomere summary, scaffold classes, and flagged scaffolds
- page 2: terminal telomere length and position by arm, and flagged telomeres
- next pages: one zoom per scaffold with its telomere blocks

Contig-terminal rows (`-n`) are drawn outline-only and left out of the counts.

ITS report. Colors follow a 3x3 key: hue is the strand (fwd blue, both purple, rev vermillion) and opacity is the canonical share.

- **Summary:** headline numbers, ITS length distribution, ITS count vs. density per scaffold (dot size is scaffold size), and ITS composition by strand and canonical share.
- **Atlas:** every scaffold with a terminal telomere or ITS, drawn to scale, homologs side by side. Up to one slide in two size groups, else a double slide in four, each group on its own axis. Rows spread to fill the page. ▼ marks the top (F) fusion candidates, (L) longest canonical ITS and (C) clusters of ITS. If rows do not fit even a double slide, the atlas becomes an unlabelled heatmap.
- **Candidates:** top clusters (>= 3 ITS within 50 kbp), candidate fusions around their junction, and the top ITS by canonical bp.
- **Loci:** one track page per top cluster, fusion, and ITS.

Candidate fusions are q→p pairs within the distance threshold with no gap between them, where at least one row is classed as fusion. They are not confirmed fusions. Fusions are sorted by shorter arm bp and clusters by summed ITS bp. Canonical bp is `(fwdCan+revCan) x` the canonical motif length, so it can exceed the ITS length when matches overlap.

## Standalone plotting

The plotting script also runs on an existing output directory:

```sh
python3 scripts/teloscope_report.py results/ -o report.pdf
python3 scripts/teloscope_report.py results/ --section its -o its.pdf
python3 scripts/teloscope_report.py results/ --png -o figures/
```

By default `-o report.pdf` writes `report_terminal.pdf` and `report_its.pdf`.

To plot a single ITS locus:

```sh
python3 scripts/plot_its.py results/ CHROM:START-END -o its.pdf
python3 scripts/plot_its.py results/ CHROM --png -o its_figures/
```

A bare `CHROM` centers on its largest ITS cluster.

Required input is `*_terminal_telomeres.bed` for the terminal report or `*_interstitial_telomeres.bed` for the ITS report. Gaps, window bedgraphs, and `*_report.tsv` are used when present.

## Standalone script options

| Flag | Meaning | Default |
| --- | --- | --- |
| `-o` | PDF stem in split mode; exact PDF path otherwise; PNG output directory | `<input_dir>/teloscope_report.pdf` stem |
| `--section` | `split`, `all` (combined), `terminal`, or `its` | `split` |
| `--png` | write one PNG per page instead of one PDF | `false` |
| `--dpi` | raster DPI | `450` |
| `--draft` | use 150 DPI for fast iteration | `false` |

## Requirements

Python 3 with `matplotlib` 3.5 or newer, `numpy`, and `pandas`.

## Regression check

```sh
python3 scripts/test_teloscope_report.py
```
