# Teloscope
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.cell.2026.07.018.svg)](https://doi.org/10.1016/j.cell.2026.07.018)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Teloscope scans assembly ends for telomeric repeats. It reads `FASTA`, `FASTA.gz`, and `GFA` inputs, merges repeat matches into telomere blocks, classifies scaffolds in FASTA mode, and writes files that are easy to inspect in BED, TSV, BEDgraph, PDF, or GFA form. FASTQ and BAM input is detected automatically (no flag): Teloscope writes the telomeric reads, a per-read telomere BED, and a report with an alignment-free, per-read telomere length estimate.

In FASTA mode, Teloscope writes terminal telomere annotations, gap coordinates, and a summary table. In GFA mode, it writes an annotated graph for BandageNG. Synthetic telomere nodes attach to the assembly with `L` links at `0M` overlap, the direct adjacency a cap represents, so BandageNG draws them as caps. `J` jump records stay reserved for real assembly gaps.

## What Teloscope writes

- FASTA mode: `*_terminal_telomeres.bed`, `*_interstitial_telomeres.bed`, `*_gaps.bed`, `*_report.tsv`, plus optional window tracks, match BED files, an optional telomere FASTA (`-a`), and a PDF report.
- GFA mode: `<input>.telo.annotated.gfa` with telomere placeholder segments linked to the original graph, plus `<input>.telo.annotated.colors.csv` that paints the caps green for BandageNG.
- Reads mode (FASTQ or BAM input, detected automatically): the telomeric reads (`*_telomeric.fastq` or `*_telomeric.bam`, unchanged records), `*_terminal_telomeres.bed` with one row per read telomere, and a `*_report.tsv` summary; an alignment-free estimate (see [Parameters](docs/parameters.md#reads-mode)).

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
| Write every optional output and the reports | `teloscope asm.fa -o results/ -r -g -e -m -a -i --plot-report` |
| Switch to a plant canonical repeat | `teloscope asm.fa -c CCCTAAA` |
| Search explicit motif variants | `teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG` |
| Annotate a graph for BandageNG | `teloscope asm.gfa -o results/` |
| Subset telomeric reads and estimate their telomere length, from FASTQ | `teloscope reads.fq.gz -j 32 -o results/` |
| Subset telomeric reads and estimate their telomere length, from BAM | `teloscope reads.bam -j 32 -o results/` |
| Also report telomeres at contig ends, e.g. before manual curation | `teloscope asm.fa -n` |
| Keep only records named like the longest one | `teloscope asm.fa --chr-only` |
| Read decompressed stdin | `zcat asm.fa.gz \| teloscope -o results/` |

Notes:

- Teloscope always searches both each input pattern and its reverse complement.
- If `-p` is omitted, Teloscope derives the search set from `-c`.
- Any of `-r`, `-g`, `-e`, `-m`, or `-i` forces the full scan instead of the fast end-only scan. In fast mode `-n` reads both end windows of every contig and adds contig-terminal rows to the terminal BED.
- GFA mode attaches telomere caps with `L` links at `0M` overlap; `J` records stay reserved for real assembly gaps.
- Gzipped FASTA/FASTQ on stdin or a pipe is not supported: pass the file, or decompress before piping. BAM on stdin or a pipe is supported (BAM mode reads BGZF directly).
- FASTQ or BAM input is detected by content, not by flag or file extension: `@` is FASTQ and the BAM magic is BAM, read from a file's inflated first bytes or from the raw first bytes of stdin or a pipe (a FIFO, `<(...)`, `/dev/stdin`), which reads like stdin. Anything else must start like FASTA or GFA or the run exits 1; only the first bytes are checked, so a GFA with missing columns can still stop the GFA reader. See [Parameters](docs/parameters.md#reads-mode).
- Reads mode always writes files under `-o` (default: the input's own directory; `.` for stdin or a pipe, so `./stdin_*` or `./<name>_*`); the report's rows are also printed to stdout. `-l` (default `300`, same as assembly) sets the measured BED rows and report only, adding reads to the subset only when set below the fixed 42 bp floor. See [Parameters](docs/parameters.md#reads-mode).
- Reads mode is an alignment-free estimate: reads pool chromosome ends by coverage, and a read broken inside a telomere looks complete and pulls the estimate down.
- BAM support has no external bioinformatics runtime dependency: it uses `zlib` directly and does not require HTSlib, `samtools`, or another converter.

## Filter assembly records

Filtering is off by default. The include/exclude flags can be repeated, and matching is case-sensitive. `--chr-only` is a single on/off switch.

| Flag | Effect | Default |
| --- | --- | --- |
| `--include-bed FILE` | keep IDs listed in column 1 | unset |
| `--exclude-bed FILE` | remove IDs listed in column 1 | unset |
| `--include-prefix LIST` | keep IDs with any comma-separated prefix | unset |
| `--exclude-prefix LIST` | remove IDs with any comma-separated prefix | unset |
| `--chr-only` | keep records named like the longest one | `false` |

`--chr-only` combines with the other filters. See [Parameters](docs/parameters.md#assembly-record-filters) for the naming rule.

Includes form a union and exclusions run last. Without an include flag, all records start selected. Prefixes are literal strings. Selector files accept one ID per line or BED3+ rows; column 1 selects a whole record, and BED coordinates never crop sequences.

FASTA matching uses the first token after `>`. GFA1 matching uses `P` names or, in a pathless graph, `S` names. Filters reject FASTQ and BAM input. Every selector must match, and an empty selection fails. See [Parameters](docs/parameters.md#assembly-record-filters) for validation and GFA limits.

For NCBI FASTA, prefer exact accession.version IDs from the [genome sequence report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/) or [assembly report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/). `GCA_` and `GCF_` name assemblies, not FASTA records, and prefixes such as `CM` or `NC_` are not universal chromosome tests.

## Typical output layout

FASTA run:

```text
results/
  asm.fa_terminal_telomeres.bed
  asm.fa_terminal_telomeres.fa
  asm.fa_interstitial_telomeres.bed
  asm.fa_gaps.bed
  asm.fa_report.tsv
  asm.fa_window_repeat_density.bedgraph
  asm.fa_window_canonical_ratio.bedgraph
  asm.fa_window_strand_ratio.bedgraph
  asm.fa_window_gc.bedgraph
  asm.fa_window_entropy.bedgraph
  asm.fa_canonical_matches.bed
  asm.fa_noncanonical_matches.bed
  asm.fa_plot_report_terminal.pdf
  asm.fa_plot_report_its.pdf
```

GFA run:

```text
results/
  asm.gfa.telo.annotated.gfa
  asm.gfa.telo.annotated.colors.csv
```

Reads run (FASTQ):

```text
results/
  reads.fq.gz_telomeric.fastq
  reads.fq.gz_terminal_telomeres.bed
  reads.fq.gz_report.tsv
```

Reads run (BAM; the subset file drops the input's extension):

```text
results/
  reads_telomeric.bam
  reads.bam_terminal_telomeres.bed
  reads.bam_report.tsv
```

## Documentation

| Page | Covers |
| --- | --- |
| [Parameters](docs/parameters.md) | command-line flags, defaults, and flag interactions |
| [Outputs](docs/outputs.md) | every output file, naming rules, and column layouts |
| [Classification](docs/classification.md) | FASTA scaffold classes and the `granular` labels |
| [Algorithm](docs/algorithm.md) | how FASTA mode and GFA mode are processed |
| [Report generation](docs/report.md) | `--plot-report`, ITS plotting, standalone plotting, and report inputs |
| [Simulation](docs/simulation.md) | the synthetic benchmark generator and evaluator |
| [Testing](docs/testing.md) | invariant and intent checks, validator runs, and test regeneration |
| [Troubleshooting](docs/troubleshooting.md) | common build, input, and runtime failures |
| [Release checklist](docs/release.md) | GitHub, Bioconda, and Zenodo release steps |
| [Validation format](https://github.com/vgl-hub/teloscope/blob/main/validateFiles/README.md) | the `.tst` harness, including directive-mode GFA cases |

## Repo layout

- `src/` and `include/`: main C++ implementation and headers
- `scripts/`: Python and shell helpers, including `teloscope_report.py`, `plot_its.py`, and plotting regression checks
- `docs/`: user-facing documentation
- `testFiles/`: public FASTA and GFA fixtures plus expected outputs
- `validateFiles/`: `.tst` manifests used by `teloscope-validate`
- `gfalibs/`: graph I/O submodule used for GFA parsing and writing

## Validation

Run the fixture, intent and invariant checks:

```sh
make test-synthetic
```

Build the validator and run the checked-in `.tst` suite:

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

Run the reads-mode (FASTQ/BAM) regression scripts:

```sh
make test-bam
make test-read-tl
```

## Citation

Teloscope is part of the gfastar tool suite.

If you use Teloscope, cite:

The complete genome of a songbird  
Giulio Formenti, Nivesh Jain, Jack A. Medico, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Ines Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceicao, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O'Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnea Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis  
Cell, Volume 189, Issue 16, pages 4922-4945.e12, August 06, 2026  
https://doi.org/10.1016/j.cell.2026.07.018
