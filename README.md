Teloscope
============
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2025.10.14.682431.svg)](https://doi.org/10.1101/2025.10.14.682431)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Introduction
------------
Teloscope is a comprehensive telomere annotation tool. It rapidly matches, counts, and reports telomeric repeats from genome assemblies (.fa), (.fa.gz), or (.gfa). Teloscope uses prefix trees and sliding windows to scan for user-defined or default repeat patterns, groups matches into telomere blocks, and classifies each chromosome by its telomere configuration. It reports all annotations in BED/BEDgraph files and produces a summary report to stdout.

To install **teloscope** from source:
```sh
git clone https://github.com/vgl-hub/teloscope.git --recursive
cd teloscope
make -j
```

Or from Bioconda:
```sh
conda install -c bioconda teloscope
```

Requirements: C++17 compiler, zlib, pthreads. The `gfalibs` submodule is fetched with `--recursive`. If missing, run `git submodule update --init`.

Usage
------------
    teloscope input.[fa|fa.gz|gfa] [options]
    teloscope -f input.[fa|fa.gz|gfa] -o output/dir/ [options]

The input file can be passed as the first positional argument or with `-f`.

**Note:** Teloscope automatically explores the input repeats and their reverse complements. If no patterns are provided, it scans for the canonical repeat and its reverse complement (default: CCCTAA/TTAGGG). Patterns accept IUPAC ambiguity codes (e.g., `NNNGGG` expands to all 64 combinations of `[ACGT][ACGT][ACGT]GGG`).

Examples
------------

* Example: Minimal run using Teloscope.

        teloscope asm.fa

* Example: Custom canonical repeat and non-canonical variants.

        teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG

* Example: Plant telomere (7-mer). No `-p` needed since patterns are derived from `-c`.

        teloscope asm.fa -c CCCTAAA

* Example: Edit distance mode to capture variants within Hamming distance 1.

        teloscope asm.fa -x 1

* Example: Set window and step sizes. Calculating window metrics.

        teloscope -f asm.fa -o results/ -j 16 -c TTAGGG -p NNNGGG -w 2000 -s 1000 -g -e -r --verbose

* Example: Full output with all optional flags.

        teloscope -f asm.fa -o results/ -j 16 -c TTAGGG -p TBAGGG,TTRGGG,YTAGGG -w 2000 -s 1000 -d 200 -l 1000 -x 1 -r -g -e -m -i -n -t 50000 --verbose --cmd

> **Why `-c` and `-p` are separate.**
> The canonical pattern (`-c`) defines *which* repeat is the reference motif. It determines how matches are classified as canonical vs. non-canonical, and how strand orientation is assigned (`p` vs. `q` labels). The pattern list (`-p`) defines *what* to search for. When `-p` is not provided, Teloscope derives the search patterns from `-c` and its reverse complement automatically. This separation lets you, for example, set `-c TTAGGG` as the reference while searching for dozens of variant patterns via `-p`, or switch to a plant canonical (`-c CCCTAAA`) without manually listing patterns. If you set `-c` to a non-default motif, you do **not** need to also pass `-p` unless you want to search for additional variants beyond the canonical pair.

Parameters
------------

To check out all options and flags, please use:
`teloscope -h`

### Input/Output

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-f` | `--input-sequence` | Input FASTA/GFA file (or pass as first positional argument) | *required* |
| `-o` | `--output` | Output directory | Input file directory |
| `-j` | `--threads` | Number of threads | All available |

### Pattern Configuration

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-c` | `--canonical` | Canonical telomere repeat | `TTAGGG` |
| `-p` | `--patterns` | Comma-separated list of repeat patterns to search | Canonical and its reverse complement |
| `-x` | `--edit-distance` | Hamming distance for pattern expansion (0, 1, or 2) | `0` |

Both the provided patterns and their reverse complements are searched automatically. When `-x` is set, all substitution variants of each seed pattern are generated and included in the search.

### Sliding Window

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-w` | `--window` | Window size (bp) | `1000` |
| `-s` | `--step` | Step size (bp) | `1000` (non-overlapping) |

**Note:** Window step defaults to window size (non-overlapping windows). This produces valid single-value BEDgraph files compatible with IGV and UCSC genome browsers. Use `-s` to set a smaller step for overlapping windows if needed.

### Block Detection

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-k` | `--max-match-distance` | Maximum distance (bp) for merging individual matches | `50` |
| `-d` | `--max-block-distance` | Maximum gap (bp) for extending adjacent blocks | `200` |
| `-l` | `--min-block-length` | Minimum block length (bp) to report | `500` |
| `-y` | `--min-block-density` | Minimum repeat density within a block (0.0 to 1.0) | `0.5` |
| `-t` | `--terminal-limit` | Distance from contig ends (bp) defining scaffold-terminal region | `50000` |

### Output Control

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-r` | `--out-win-repeats` | Output canonical/noncanonical repeats and density by window | `false` |
| `-g` | `--out-gc` | Output GC content for each window | `false` |
| `-e` | `--out-entropy` | Output Shannon entropy for each window | `false` |
| `-m` | `--out-matches` | Output all canonical and terminal non-canonical matches | `false` |
| `-i` | `--out-its` | Output assembly interstitial telomere (ITS) regions | `false` |
| `-u` | `--ultra-fast` | Only scan terminal telomeres at contig ends | `true` |
| `-n` | `--manual-curation` | Retain all terminal telomeres (contig + scaffold) in BED output | `false` (scaffold only) |
| `-v` | `--version` | Print current software version | |
| `-h` | `--help` | Print help | |
| | `--verbose` | Verbose output | |
| | `--cmd` | Print command line | |

**Note:** Enabling any genome-wide output flag (`-r`, `-g`, `-e`, `-m`, `-i`) automatically disables ultra-fast mode.

Outputs
------------
Teloscope always outputs telomere annotations in BED format:

* `terminal_telomeres.bed` Annotation of the full telomere blocks in the assembly. Columns: chr, start, end, length, label, fwdCount, revCount, canonCount, nonCanonCount, chrSize, type. The label column reports arm assignment: `p` (forward/p-arm, >66.6% forward matches), `q` (reverse/q-arm, <33.3% forward), or `b` (balanced). The type column reports `scaffold` or `contig` terminal classification. By default only scaffold-terminal telomeres are reported; use `--manual-curation` to include contig-terminal blocks.

Additional optional outputs (disabled in ultra-fast mode):

* `window_fwd.bedgraph` Forward repeat counts per window (`-r`).
* `window_rev.bedgraph` Reverse repeat counts per window (`-r`).
* `window_canonical.bedgraph` Canonical repeat counts per window (`-r`).
* `window_noncanonical.bedgraph` Non-canonical repeat counts per window (`-r`).
* `window_gc.bedgraph` GC content per window (`-g`).
* `window_entropy.bedgraph` Shannon entropy per window (`-e`).
* `canonical_matches.bed` Coordinates of canonical repeats throughout the assembly (`-m`).
* `noncanonical_matches.bed` Coordinates of non-canonical repeats in terminal regions (`-m`).
* `interstitial_telomeres.bed` Blocks of adjacent canonical repeat matches outside terminal regions, representing interstitial telomeres (ITSs) (`-i`).

All output filenames are prefixed with the input filename (e.g., `asm.fa_terminal_telomeres.bed`).

Console summary
------------
Teloscope prints a tab-delimited **Path Summary Report** followed by an **Assembly Summary Report** to stdout.

In ultra-fast mode, the per-path columns are: `pos`, `header`, `telomeres`, `labels`, `gaps`, `type`, `granular`. In full-scan mode, three additional columns are appended: `its`, `canonical`, `windows`.

The Assembly Summary includes total paths, gaps, and telomeres; telomere length statistics (mean, median, min, max); chromosome counts by telomere number (two, one, zero); and chromosome classification counts (see below).

Chromosome classification
------------
Teloscope classifies each chromosome based on its terminal telomere blocks and gap content:

| Type | Gapped variant | Description |
|------|---------------|-------------|
| `t2t` | `gapped_t2t` | Both p-arm and q-arm telomeres detected in correct orientation |
| `incomplete` | `gapped_incomplete` | Only one arm (p or q) has a scaffold-terminal telomere |
| `none` | `gapped_none` | No scaffold-terminal telomeres detected |
| `misassembly` | `gapped_misassembly` | Two telomeres in incorrect positional order, or duplicate arms |
| `discordant` | `gapped_discordant` | Terminal telomere with mixed forward/reverse orientation |

Gapped variants indicate the sequence contains internal gaps (Ns). The `granular` column in the path summary provides a per-block label string (e.g., `Pq` indicates a longest p-arm block plus a secondary q-arm block; `*` marks blocks with discordant orientation).

How it works
------------
Briefly, **Teloscope** reads an assembly and decomposes its parts. Input patterns (with optional IUPAC codes and edit distance) are expanded into all concrete variants along with their reverse complements, and a prefix tree (trie) is built for simultaneous multi-pattern matching. Each contig is scanned using the trie: in ultra-fast mode (default), only the terminal regions within `-t` bp of each end are scanned; otherwise, the full sequence is analyzed.

Individual matches within `-k` bp of each other are merged, and adjacent merged regions within `-d` bp are extended into telomere blocks. Blocks shorter than `-l` bp or with repeat density below `-y` are filtered out. Each block is labeled by strand orientation (`p`, `q`, or `b`), and blocks are used to classify the chromosome as T2T, incomplete, misassembled, discordant, or none. Terminal telomere BED annotations are always written. Optional window metrics, match coordinates, and ITS blocks are produced based on the enabled flags.

Simulation and validation
------------
Teloscope includes a simulation framework for benchmarking detection accuracy across mutation rates. Build with:
```sh
make simulate
```

This produces `teloscope-simulate`, a standalone tool with two modes:

**Generate mode** produces synthetic multi-contig FASTA files with known telomere positions:
```sh
teloscope-simulate -n 1000 -r 1e-4 -s 42 -o testFiles/simulate/rate_1e-4
```

Each synthetic sequence has the structure:
```
[p-canonical: R x CCCTAA] [p-TVR: T x HD1(CCCTAA)] [central: random ACGT] [q-TVR: T x HD1(TTAGGG)] [q-canonical: R x TTAGGG]
```

TVR (telomere variant repeat) regions contain repeats with exactly one substitution from the canonical motif, so they are detectable at `-x 1` but not at `-x 0`. Outputs: `sequences.fa` (synthetic FASTA) and `ground_truth.tsv` (per-contig telomere coordinates).

**Evaluate mode** compares teloscope BED output against the ground truth:
```sh
teloscope-simulate --evaluate -g ground_truth.tsv -b terminal_telomeres.bed
```

Output columns: `total_tips`, `detected`, `sensitivity`, `mean_bias_bp`, `mean_abs_err_bp`, `tvr_rate`, `fp_blocks`.

### Simulation parameters

| Flag | Description | Default |
|------|-------------|---------|
| `-n` | Number of sequences | `1000` |
| `-R` | Canonical repeats per tip | `2000` |
| `-T` | TVR repeats per tip | `100` |
| `-L` | Central segment length (bp) | `25000` |
| `-r` | Per-base mutation rate | `0.0` |
| `-s` | Random seed | `42` |
| `-o` | Output directory | `testFiles/simulate` |

### Pipeline script

A sweep across mutation rates (1e-6 through 1e-2) can be run with:
```sh
bash .github/workflows/val-simulate.sh
```

Configurable via environment: `SIM_N=1000000 SIM_SEED=42 bash .github/workflows/val-simulate.sh`

Testing
------------
Tests live under `testFiles/` and `validateFiles/`. Build and run with:
```sh
make validate
build/bin/teloscope-validate validateFiles
```

To regenerate expected test outputs from a known-good build:
```sh
make regenerate
build/bin/teloscope-generate-tests
```

How to cite
------------

**Teloscope** is part of the gfastar tool suite.
If you use **Teloscope** in your research, please cite:

**The complete genome of a songbird**
Giulio Formenti, Nivesh Jain, **Jack A. Medico**, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Ines Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceicao, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O'Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnea Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis
bioRxiv 2025.10.14.682431; doi: https://doi.org/10.1101/2025.10.14.682431
