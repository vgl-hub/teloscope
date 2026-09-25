[Back to README](index.md)

# Report generation

Teloscope can generate separate terminal and interstitial (ITS) PDF reports during a FASTA run:

```sh
teloscope asm.fa -o results/ -r -e -g -i --plot-report
```

This writes `results/asm.fa_plot_report_terminal.pdf` and `results/asm.fa_plot_report_its.pdf`. Without an ITS BED no ITS report is written.

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

The ITS report has its own pages:

- **Distributions:** canonical bp versus row length (hexbin above 500 rows), length ECDFs by orientation, and class counts.
- **Atlas:** every scaffold with a terminal or ITS call, longest first, 20 per page; 100 bins per scaffold colored by ITS count on one shared log scale. Grey means no call; white means outside the initial end windows, not proven unscanned.
- **Selected candidates:** the top five rows by canonical bp, clusters (>= 3 rows within 50 kb), and candidate fusion pairs.
- **Selected loci:** one track page each for the top cluster, row, and fusion pair.

Canonical bp is the canonical match count times the motif length. Its ratio to row length (`can_prop`) can exceed one when matches overlap. Candidate fusions are q→p pairs within the distance threshold with an engine fusion label and no N-gap between them; they are not confirmed fusions. Missing metadata falls back to a 6 bp motif and a 1,000 bp threshold, and the report says so.

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
python3 scripts/teloscope_report.py results/ --section terminal -o terminal.pdf
python3 scripts/teloscope_report.py results/ --section its -o its.pdf
python3 scripts/teloscope_report.py results/ --section all -o combined.pdf
python3 scripts/teloscope_report.py results/ --png -o figures/
```

The script auto-detects the Teloscope files in that directory. By default `-o report.pdf` is a stem for `report_terminal.pdf` and `report_its.pdf`; `--section all` writes one combined PDF.

ITS output also includes `*_its_rows.tsv` (every accepted row), `*_its_scaffolds.tsv` (per-scaffold counts and display labels), and `*_its_top_hits.tsv`. Malformed BED rows are skipped with a warning.

For a single interstitial telomeric sequence locus, use `plot_its.py` on the same output directory:

```sh
python3 scripts/plot_its.py results/ CHROM:START-END -o its.pdf
python3 scripts/plot_its.py results/ CHROM --png -o its_figures/
```

A bare `CHROM` auto-centers on the largest interstitial telomere cluster on that scaffold.

Contig-terminal rows (written with `-n`) are excluded from telomere counts and length panels on the report, and are drawn outline-only on the terminal zoom pages.

Minimum required input:

- `*_terminal_telomeres.bed` for a terminal report, or
- `*_interstitial_telomeres.bed` for an ITS report (terminal BED optional).

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
| `-o` | PDF stem in split mode; exact PDF path otherwise; PNG output directory | `<input_dir>/teloscope_report.pdf` stem |
| `--section` | `split`, `all` (combined), `terminal`, or `its` | `split` |
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
