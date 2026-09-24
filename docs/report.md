[Back to README](index.md)

# Report generation

Teloscope can generate a PDF report during a FASTA run:

```sh
teloscope asm.fa -o results/ -r -e -g -i --plot-report
```

This writes `results/asm.fa_plot_report.pdf`.

`-r` is recommended with `--plot-report` so the report includes repeat-density, canonical-ratio, and strand-bias tracks. GC and entropy tracks are added automatically when their files are present.

Assembly record filters apply before these files are written. This example uses exact chromosome accession.version IDs derived from assembly metadata and omits every other record from the TSV, BED/BEDgraph files, and PDF:

```sh
teloscope asm.fa -o results/ --include-bed chromosomes.ids -r --plot-report
```

Prefix filters are also available, but database prefixes are not universal chromosome labels. Inspect the FASTA primary IDs before using `--include-prefix`.

## What the report contains

- page 1: assembly overview, scaffold classes, and flagged scaffolds
- page 2: telomere length summary and flagged telomere blocks
- next pages: one terminal zoom figure per scaffold with called telomere blocks
- last pages: three interstitial telomere (ITS) pages (when the interstitial BED has rows)

The three ITS pages, in order:
- **genome view**: an ideogram of every scaffold with a terminal telomere or ITS (split into
  long/short panels in full-scan mode, unfolded distance-to-end in fast mode), row/bp counts
  per junction class, and an ITS length distribution per strand
- **composition and top hits**: a length-vs-canonical-bp scatter (hexbin above ~20k rows) with
  the 100%-canonical diagonal and the engine's canonical-count floor, plus three ranked tables
  (longest ITS by canonical bp, ITS clusters, candidate fusions)
- **top loci**: terminal-zoom-style track columns for the top cluster, longest ITS row, and top
  candidate fusion

Each terminal zoom page can include:

- telomere block positions
- gap intervals
- repeat density
- canonical ratio
- strand bias
- GC content
- entropy

## Standalone plotting

The plotting script can be run on an existing Teloscope output directory:

```sh
python3 scripts/teloscope_report.py results/ -o report.pdf
python3 scripts/teloscope_report.py results/ --png -o figures/
```

The script auto-detects the Teloscope files in that directory.

For a single interstitial telomeric sequence locus, use `plot_its.py` on the same output directory:

```sh
python3 scripts/plot_its.py results/ CHROM:START-END -o its.pdf
python3 scripts/plot_its.py results/ CHROM --png -o its_figures/
```

A bare `CHROM` auto-centers on the largest interstitial telomere cluster on that scaffold.

Contig-terminal rows (written with `-n`) are excluded from telomere counts and length panels on the report, and are drawn outline-only on the terminal zoom pages.

Minimum required input:

- `*_terminal_telomeres.bed`

Common optional inputs:

- `*_gaps.bed`
- `*_window_repeat_density.bedgraph`
- `*_window_canonical_ratio.bedgraph`
- `*_window_strand_ratio.bedgraph`
- `*_window_gc.bedgraph`
- `*_window_entropy.bedgraph`
- `*_report.tsv`

## Standalone script options

| Flag | Meaning | Default |
| --- | --- | --- |
| `-o` | output file for PDF mode or output directory for PNG mode | `<input_dir>/teloscope_report.pdf` |
| `--png` | write one PNG per page instead of one PDF | `false` |
| `--dpi` | raster DPI | `450` |
| `--draft` | use 150 DPI for fast iteration | `false` |

## Requirements

- Python 3
- `matplotlib` 3.5 or newer
- `numpy`
- `pandas`

## Regression check

The report layout regression script lives in `scripts/`:

```sh
python3 scripts/test_teloscope_report.py
```

It checks the plotting code directly and does not need a full Teloscope run on disk.
