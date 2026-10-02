# Teloscope
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.cell.2026.07.018.svg)](https://doi.org/10.1016/j.cell.2026.07.018)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Teloscope is a telomere annotation tool. It rapidly matches, counts, and reports telomeric repeats in genome assemblies (FASTA, GFA) and reads (FASTQ, BAM).

## Install

```sh
conda install -c conda-forge -c bioconda teloscope
```

Or download a Linux, macOS, or Windows binary from [Releases](https://github.com/vgl-hub/teloscope/releases), or build from source with a C++17 compiler and zlib:

```sh
git clone --recursive https://github.com/vgl-hub/teloscope.git
cd teloscope
make -j
```

`--plot-report` needs Python 3 with `matplotlib`, `numpy`, and `pandas`, and finds `scripts/teloscope_report.py` in a source checkout or a Bioconda install; a release binary needs the script beside it.

## Quick start

| Task | Command |
| --- | --- |
| Scan a genome assembly (FASTA or FASTA.gz) | `teloscope asm.fa.gz` |
| Write every output and the PDF reports | `teloscope asm.fa -o results/ -r -g -e -m -a -i --plot-report` |
| Switch to a plant canonical repeat | `teloscope asm.fa -c CCCTAAA` |
| Search explicit motif variants | `teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG` |
| Add telomeres at contig ends to the terminal BED, e.g. for manual curation | `teloscope asm.fa -n` |
| Keep only records named like the longest one | `teloscope asm.fa --chr-only` |
| Keep only listed records | `teloscope asm.fa --include-bed ids.txt` |
| Annotate a graph for BandageNG | `teloscope asm.gfa -o results/` |
| Measure and subset telomeric reads (FASTQ or BAM) | `teloscope reads.fq.gz -j 32 -o results/` |
| Read decompressed stdin | <code>zcat asm.fa.gz &#124; teloscope -o results/</code> |

Each pattern is also searched as its reverse complement, and without `-p` the search set comes from `-c`. The default scan reads only sequence ends; `-r`, `-g`, `-e`, `-m`, or `-i` switch to a full scan. Compressed stdin must be BAM.

## Outputs

```text
results/
  asm.fa_terminal_telomeres.bed       telomeres at sequence ends
  asm.fa_interstitial_telomeres.bed   interstitial telomeres (ITS)
  asm.fa_gaps.bed                     assembly gaps
  asm.fa_report.tsv                   per-sequence classes and assembly summary, as on stdout
  asm.fa_terminal_telomeres.fa        -a
  asm.fa_window_*.bedgraph            -r, -g, -e
  asm.fa_*canonical_matches.bed       -m
  asm.fa_plot_report_*.pdf            --plot-report
```

A GFA run writes `asm.gfa.telo.annotated.gfa` and a BandageNG color file. A reads run writes the telomeric reads, a per-read telomere BED, and a report.

## Documentation

| Page | Covers |
| --- | --- |
| [Parameters](docs/parameters.md) | flags, defaults, record filters, and reads mode |
| [Outputs](docs/outputs.md) | every file and column |
| [Classification](docs/classification.md) | scaffold classes, anomalies, and block labels |
| [Algorithm](docs/algorithm.md) | how each mode works |
| [Report generation](docs/report.md) | PDF reports and standalone plotting |
| [Troubleshooting](docs/troubleshooting.md) | errors and surprising calls |
| [Testing](docs/testing.md) | test suites, validation, and repo layout |
| [Simulation](docs/simulation.md) | the synthetic benchmark |
| [Release checklist](docs/release.md) | GitHub, Bioconda, Galaxy, and Zenodo steps |

## Citation

Teloscope is part of the gfastar tool suite.

If you use Teloscope, cite:

The complete genome of a songbird  
Giulio Formenti, Nivesh Jain, Jack A. Medico, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Ines Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceicao, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O'Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnea Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis  
Cell, Volume 189, Issue 16, pages 4922-4945.e12, August 06, 2026  
https://doi.org/10.1016/j.cell.2026.07.018
