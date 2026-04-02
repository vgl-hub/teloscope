# Teloscope
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2025.10.14.682431.svg)](https://doi.org/10.1101/2025.10.14.682431)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Teloscope scans assembly ends for telomeric repeats. It reads `FASTA`, `FASTA.gz`, and `GFA` inputs, merges repeat matches into telomere blocks, classifies scaffolds in FASTA mode, and writes files that are easy to inspect in BED, TSV, BEDgraph, PDF, or GFA form.

In FASTA mode, Teloscope writes terminal telomere annotations, gap coordinates, and a summary table. In GFA mode, it writes an annotated graph for BandageNG. By default, synthetic telomere nodes are connected back to the assembly with GFA `J` jump records instead of hard-adjacency `L` links to reduce disruptive force-layout lines.

## What Teloscope writes

- FASTA mode: `*_terminal_telomeres.bed`, `*_gaps.bed`, `*_report.tsv`, plus optional window tracks, match BED files, ITS BED files, and a PDF report.
- GFA mode: `<input>.telo.annotated.gfa` with telomere placeholder segments and links attached to the original graph.

## Install

From source:

```sh
git clone https://github.com/vgl-hub/teloscope.git --recursive
cd teloscope
make -j
```

To build the helper binaries as well:

```sh
make all -j
```

From Bioconda:

```sh
conda install -c bioconda teloscope
```

Build requirements:

- C++17 compiler
- `zlib`
- `pthread`
- the `gfalibs` submodule

If the submodule is missing, run:

```sh
git submodule update --init --recursive
```

For `--plot-report`, install Python 3 with `matplotlib`, `numpy`, and `pandas`.

## Quick start

| Task | Command |
| --- | --- |
| Scan a vertebrate assembly with the default motif | `teloscope asm.fa` |
| Read compressed FASTA directly | `teloscope asm.fa.gz` |
| Write all optional FASTA outputs | `teloscope asm.fa -o results/ -r -g -e -m -i --plot-report` |
| Switch to a plant canonical repeat | `teloscope asm.fa -c CCCTAAA` |
| Search explicit motif variants | `teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG` |
| Annotate a graph for BandageNG | `teloscope asm.gfa -o results/` |
| Read decompressed stdin | `zcat asm.fa.gz | teloscope -o results/` |

Notes:

- Teloscope always searches both each input pattern and its reverse complement.
- If `-p` is omitted, Teloscope derives the search set from `-c`.
- Any genome-wide output flag (`-r`, `-g`, `-e`, `-m`, `-i`) disables ultra-fast mode automatically.
- GFA mode uses GFA `J` jump records for telomere connectors instead of direct `L` adjacency links.
- Gzipped stdin is not supported. Decompress before piping.

## Typical output layout

FASTA run:

```text
results/
  asm.fa_terminal_telomeres.bed
  asm.fa_gaps.bed
  asm.fa_report.tsv
  asm.fa_window_repeat_density.bedgraph
  asm.fa_window_canonical_ratio.bedgraph
  asm.fa_window_strand_ratio.bedgraph
  asm.fa_plot_report.pdf
```

GFA run:

```text
results/
  asm.gfa.telo.annotated.gfa
```

## Documentation

| Page | Covers |
| --- | --- |
| [Parameters](docs/parameters.md) | command-line flags, defaults, and flag interactions |
| [Outputs](docs/outputs.md) | every output file, naming rules, and column layouts |
| [Classification](docs/classification.md) | FASTA scaffold classes and the `granular` labels |
| [Algorithm](docs/algorithm.md) | how FASTA mode and GFA mode are processed |
| [Report generation](docs/report.md) | `--plot-report`, standalone plotting, and report inputs |
| [Simulation](docs/simulation.md) | the synthetic benchmark generator and evaluator |
| [Testing](docs/testing.md) | validator runs, report regression checks, and test regeneration |
| [Troubleshooting](docs/troubleshooting.md) | common build, input, and runtime failures |
| [Validation format](validateFiles/README.md) | the `.tst` harness, including directive-mode GFA cases |

## Repo layout

- `src/` and `include/`: main C++ implementation and headers
- `scripts/`: Python and shell helpers, including `teloscope_report.py` and its regression script
- `docs/`: user-facing documentation
- `testFiles/`: public FASTA and GFA fixtures plus expected outputs
- `validateFiles/`: `.tst` manifests used by `teloscope-validate`
- `gfalibs/`: graph I/O submodule used for GFA parsing and writing

## Validation

Build the validator and run the checked-in suite:

```sh
make validate
build/bin/teloscope-validate validateFiles
```

Run the report layout regression script:

```sh
python3 scripts/test_teloscope_report.py
```

Run the gap BED regression script:

```sh
bash scripts/test_gaps_bed.sh
```

## Citation

Teloscope is part of the gfastar tool suite.

If you use Teloscope, cite:

The complete genome of a songbird  
Giulio Formenti, Nivesh Jain, Jack A. Medico, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Ines Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceicao, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O'Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnea Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis  
bioRxiv 2025.10.14.682431  
https://doi.org/10.1101/2025.10.14.682431
