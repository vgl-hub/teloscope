[← Back to README](../README.md)

Report generation
============

A companion Python script generates publication-ready figures from Teloscope output. It requires Python 3.6+ and matplotlib (no other dependencies).

```sh
python scripts/teloscope_report.py output/dir/ -o report.pdf
```

This produces a multi-page PDF with an assembly overview (chromosome classification donut and telomere length distribution) followed by per-chromosome telomere profiles (block locations, repeat density, and strand ratio tracks). Individual PNGs can be generated instead:

```sh
python scripts/teloscope_report.py output/dir/ --png -o figures/
```

The script auto-detects all Teloscope output files in the given directory. At minimum it needs the `*_terminal_telomeres.bed` file; repeat density and strand ratio BEDgraph files are used when available. Figures follow Nature journal specifications (Arial, 7 pt, 183 mm width, 600 DPI, colorblind-safe palette).

| Flag | Description | Default |
|------|-------------|---------|
| `-o` | Output file (PDF) or directory (PNG) | `<directory>/teloscope_report.pdf` |
| `--png` | Save individual PNG files instead of a single PDF | `false` |
| `--dpi` | DPI for raster output | `600` |
