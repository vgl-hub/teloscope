[Back to README](index.md)

# Report generation

Teloscope can generate separate terminal and interstitial (ITS) PDF reports during a FASTA run:

```sh
teloscope asm.fa -o results/ -r -e -g -i --plot-report
```

This writes `results/asm.fa_plot_report_terminal.pdf` and
`results/asm.fa_plot_report_its.pdf`. An empty ITS BED produces an explicit
zero-observation report. A missing ITS BED means no ITS report is written.

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

The ITS report has its own page sequence:

- **Distributions:** canonical match bp versus row length, exact cumulative length
  distributions by orientation, and class counts and summed row lengths. Above 500
  positive observations, the composition plot uses hexagons with a logarithmic count
  scale. Every positive observation contributes; no random subsampling is used.
  Zero canonical counts are counted explicitly and excluded from logarithmic axes.
- **Atlas:** every scaffold represented by a terminal or ITS call, in decreasing
  length order, paginated at 20 scaffolds per page. This is not an inventory of
  scaffolds with no calls. Each row has 100 equal-width bins over its plotted extent;
  color shows the number of ITS row midpoints in each bin, with one shared logarithmic
  count scale across all atlas pages. The axis shows relative position and each row
  states its extent in bp/kb/Mb. Thus short scaffolds remain visible. Grey means no
  observed call; white marks positions outside the initial end windows. The engine
  can extend those windows while building terminal arrays, so white does not prove
  that a region was unscanned. Exact scanned extents are not recorded in the inputs.
- **Selected candidates:** up to five canonical-bp-ranked rows, five clusters, and
  five candidate fusion pairs. Rankings use deterministic scaffold/coordinate tie
  breaks. Clusters join rows separated by at most 50 kb and require at least three
  rows; these are spatial groups, not independent biological events.
- **Selected loci:** one terminal-style track page each for the top cluster, top
  canonical-bp-ranked row, and top candidate fusion pair, when present. These
  selected views can refer to the same region.

ITS pages use the terminal report's maximum-track geometry (7.20 × 3.70 inches),
palette, editable PDF text, and a minimum text size of 5 pt. Balanced orientation
uses the terminal yellow. Orientation ECDFs also use distinct line styles.

Canonical match bp is the canonical match count multiplied by the motif length.
The table's ratio (TSV field `can_prop`) divides this estimate by row length;
it can exceed one for overlapping matches and is not a genomic coverage fraction.
The diagonal indicates equal bp, not proven repeat purity. The configured count
floor is shown only when both motif and threshold metadata are available. Missing
motif metadata uses an explicitly disclosed six-base assumption. Missing scan
metadata is labeled unknown. Missing or conflicting scaffold sizes leave relative
coordinates undefined in the TSV, and plotted extents are labeled as inferred.

Candidate fusion pairs require q→p order, the configured distance threshold, at
least one engine fusion label, and no N-gap overlapping a positive spacer.
Overlapping arrays retain the engine's zero-distance convention. If the threshold
is missing, the report discloses the 1,000 bp default. Candidates are not confirmed
chromosome fusions.

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

The script auto-detects the Teloscope files in that directory. In the default split
mode, `-o report.pdf` supplies the stem for `report_terminal.pdf` and
`report_its.pdf`; `--section all` explicitly requests one combined PDF.
Combined mode omits empty ITS pages when a terminal report is available.
PNG mode writes one image per page, with distinct ITS atlas and locus filenames.

ITS output includes `*_its_rows.tsv` with every accepted row,
`*_its_scaffolds.tsv` with complete scaffold counts and the display-label mapping,
and `*_its_top_hits.tsv` with all candidate pairs, the top 25 rows, and all
qualifying clusters. Long scaffold identifiers receive unique display aliases;
full identifiers remain in the TSVs. Malformed BED rows are skipped individually
with diagnostics, using the same parsing rules for statistics and locus tracks.

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

It checks parsing, numerical accounting, section selection, geometry, PDF fonts,
and sparse/dense rendering without a full Teloscope run. To save visual proofs:

```sh
TELOSCOPE_REPORT_PROOFS=results/report_proofs python3 scripts/test_teloscope_report.py
```
