[← Back to README](../README.md)

Report generation
============

The `--report` flag generates a publication-ready PDF after analysis:

```sh
teloscope asm.fa -o results/ -r --report
```

This produces `results/asm.fa_report.pdf` with an assembly overview (chromosome classification summary and telomere length distribution) followed by terminal zoom figures for each chromosome that has telomere blocks. Each figure shows two side-by-side panels (p-end and q-end) with block locations, repeat density, canonical ratio, and strand ratio tracks. Use `-r` alongside `--report` for the density and ratio tracks.

Requires Python 3.6+ with matplotlib, numpy, and pandas.

### Manual invocation

The report script can also be run separately:

```sh
python scripts/teloscope_report.py results/ -o report.pdf
python scripts/teloscope_report.py results/ --png -o figures/
```

The script auto-detects all Teloscope output files in the given directory. At minimum it needs the `*_terminal_telomeres.bed` file; repeat density, canonical ratio, and strand ratio BEDgraph files are used when available. All rendering is rasterized for speed, and figures are saved and closed one at a time to keep memory usage low. Figures follow Nature journal specifications (Arial, 7 pt, 183 mm width, 450 DPI, colorblind-safe palette).

| Flag | Description | Default |
|------|-------------|---------|
| `-o` | Output file (PDF) or directory (PNG) | `<directory>/teloscope_report.pdf` |
| `--png` | Save individual PNG files instead of a single PDF | `false` |
| `--dpi` | DPI for raster output | `450` |
| `--draft` | Render at 150 DPI for fast iteration | `false` |
