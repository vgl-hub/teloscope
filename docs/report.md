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

- page 1: assembly telomere summary tiles (median assembled array first), scaffold classes, and flagged scaffolds
- page 2: terminal telomere length by arm and positioning on one kbp scale, and flagged telomeres
- next pages: one terminal zoom figure per scaffold with called telomere blocks

The ITS report has its own pages. Colors follow one 3x3 key: hue is the strand (fwd blue, both purple, rev vermillion) and opacity is the canonical share. Forward (CCCTAA, the 5' end) is always drawn first, on the left.

- **Summary:** headline tiles, the ITS length distribution, ITS count against ITS density per scaffold (dot size follows log10 scaffold size, so dot chromosomes stand out), and composition per ITS.
  Panel c combines forward share (mirrored) versus canonical share with marginal histograms and ITS counts and % of ITS bp inside each of the nine classes. Both shares use the same thirds rule as `--label-threshold`, so the strand column matches the p/q/b label.
- **Atlas:** every scaffold with a terminal or ITS call, drawn to scale, homologs (mat/pat, hap1/hap2) side by side, in two columns at a fixed row height. The page grows to a double slide when one slide is not enough, and only then splits into more pages. Each ITS is one cell colored by its class. Scaffolds are grouped into size tiers, each with its own panel and scale, so every scaffold spans at least a quarter of its axis (macro, micro and dot chromosomes, or chromosomes and unplaced scaffolds, split at their widest size gaps). A tier keeps its scale across pages. A ▼ marks the locus of each locus page (C cluster, F fusion, L longest canonical ITS; ranks 1–5 match the candidate and locus pages). White means outside the initial end windows, not proven unscanned.
- **Candidates:** the clusters with most ITS bp (>= 3 ITS within 50 kbp), candidate fusion pairs drawn around their junction (q block left, p block right), and the top five ITS by canonical bp split into fwd/rev canonical and non-canonical.
- **Loci:** one track page for each of the top five clusters, fusion pairs, and ITS (or all available when fewer than five), with the same tracks as the terminal zoom pages.

Canonical bp is the canonical match count times the motif length. Its ratio to ITS length (`can_prop`) can exceed one when matches overlap. Candidate fusions are q→p pairs within the distance threshold with an engine fusion label and no N-gap between them; they are not confirmed fusions. Missing metadata falls back to a 6 bp motif and a 1,000 bp threshold, and the report says so.

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
