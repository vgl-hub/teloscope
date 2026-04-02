[Back to README](../README.md)

# Report generation

Teloscope can generate a PDF report during a FASTA run:

```sh
teloscope asm.fa -o results/ -r --plot-report
```

This writes `results/asm.fa_plot_report.pdf`.

`-r` is recommended with `--plot-report` so the report includes repeat-density, canonical-ratio, and strand-bias tracks. GC and entropy tracks are added automatically when their files are present.

## What the report contains

- page 1: assembly overview, scaffold classes, and flagged scaffolds
- page 2: telomere length summary and flagged telomere blocks
- later pages: one terminal zoom figure per scaffold with called telomere blocks

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
- `matplotlib`
- `numpy`
- `pandas`

## Regression check

The report layout regression script lives in `scripts/`:

```sh
python3 scripts/test_teloscope_report.py
```

It checks the plotting code directly and does not need a full Teloscope run on disk.
