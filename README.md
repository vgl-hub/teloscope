Teloscope
============
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2025.10.14.682431.svg)](https://doi.org/10.1101/2025.10.14.682431)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Teloscope is a comprehensive telomere annotation tool. It rapidly matches, counts, and reports telomeric repeats from genome assemblies (.fa), (.fa.gz), or (.gfa). It groups repeat matches into telomere blocks, classifies each chromosome by its telomere completeness, and outputs annotations in standard BED/BEDgraph files along with a summary report to stdout.

Install
------------

From source:
```sh
git clone https://github.com/vgl-hub/teloscope.git --recursive
cd teloscope
make -j
```

From Bioconda:
```sh
conda install -c bioconda teloscope
```

Requirements: C++17 compiler, zlib, pthreads. The `gfalibs` submodule is fetched with `--recursive`. If missing, run `git submodule update --init`.

Quick start
------------

### Default (vertebrate TTAGGG)
    teloscope asm.fa

### Non-default organisms
    teloscope asm.fa -c CCCTAAA

### Fuzzy matching
    teloscope asm.fa -x 1
    teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG

### Full analysis
    teloscope -f asm.fa -o results/ -r -g -e -m -i -n

**Note:** Teloscope always searches for both the input patterns and their reverse complements. If no patterns are provided (`-p`), it defaults to the canonical repeat (`-c`) and its reverse complement. Patterns accept IUPAC ambiguity codes (e.g., `NNNGGG` expands to all 64 combinations of `[ACGT][ACGT][ACGT]GGG`).

Documentation
------------

- [Parameters](docs/parameters.md) — all flags with defaults, `-c` vs `-p` explainer
- [Outputs](docs/outputs.md) — output files, column specs, console summary
- [Classification](docs/classification.md) — chromosome classification table and decision logic
- [Algorithm](docs/algorithm.md) — how Teloscope works under the hood
- [Report generation](docs/report.md) — publication-ready figures from Teloscope output
- [Simulation](docs/simulation.md) — benchmarking framework for detection accuracy
- [Testing](docs/testing.md) — test suite, `make validate`, `make regenerate`

How to cite
------------

**Teloscope** is part of the gfastar tool suite.
If you use **Teloscope** in your research, please cite:

**The complete genome of a songbird**
Giulio Formenti, Nivesh Jain, **Jack A. Medico**, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Ines Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceicao, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O'Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnea Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis
bioRxiv 2025.10.14.682431; doi: https://doi.org/10.1101/2025.10.14.682431
