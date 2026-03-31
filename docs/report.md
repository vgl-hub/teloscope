[← Back to README](../README.md)

Report generation
============

The `--plot-report` flag generates a publication-ready PDF after analysis:

```sh
teloscope asm.fa -o results/ -r --plot-report
```

This produces `results/asm.fa_plot_report.pdf` with two overview pages followed by terminal zoom figures for each chromosome that has telomere blocks. The first overview page summarizes scaffold classification and flagged scaffolds, and the second summarizes telomere length and flagged telomeres. Terminal zoom figures show side-by-side p-end and q-end panels with block locations and any available repeat-density, canonical-ratio, strand-bias, GC-content, and entropy tracks. Use `-r` alongside `--plot-report` for the density and ratio tracks.

Requires Python 3.6+ with matplotlib, numpy, and pandas.

### Manual invocation

The report script can also be run separately:

```sh
python scripts/teloscope_report.py results/ -o report.pdf
python scripts/teloscope_report.py results/ --png -o figures/
```

The script auto-detects all Teloscope output files in the given directory. At minimum it needs the `*_terminal_telomeres.bed` file; repeat density, canonical ratio, strand ratio, GC, and entropy BEDgraph files are used when available. Figures are saved and closed one at a time to keep memory usage low. PDF export keeps text and core structural marks editable, while dense point layers may still be rasterized when appropriate. Figures follow Nature journal specifications (Arial/Helvetica-style sans serif, 5-7 pt body text, 183 mm width, 450 DPI raster output when needed, and a colorblind-safe palette).

| Flag | Description | Default |
|------|-------------|---------|
| `-o` | Output file (PDF) or directory (PNG) | `<directory>/teloscope_report.pdf` |
| `--png` | Save individual PNG files instead of a single PDF | `false` |
| `--dpi` | DPI for raster output | `450` |
| `--draft` | Render at 150 DPI for fast iteration | `false` |
