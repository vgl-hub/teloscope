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
Teloscope is a comprehensive telomere annotation tool. It rapidly matches, counts, and reports telomeric repeats from genome assemblies (.fa), (.fa.gz), or (.gfa). It groups repeat matches into telomere blocks, classifies each chromosome by its telomere completeness, and outputs annotations in standard BED/BEDgraph files along with a summary report to stdout.

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

**Note:** Teloscope always searches for both the input patterns and their reverse complements. If no patterns are provided (`-p`), it defaults to the canonical repeat (`-c`) and its reverse complement. Patterns accept IUPAC ambiguity codes (e.g., `NNNGGG` expands to all 64 combinations of `[ACGT][ACGT][ACGT]GGG`).

Examples
------------

* Example: Minimal run.

        teloscope asm.fa

* Example: Custom canonical repeat and non-canonical variants.

        teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG

* Example: Plant telomere (7-mer). No `-p` needed; patterns are derived from `-c` automatically.

        teloscope asm.fa -c CCCTAAA

* Example: Allow up to 1 mismatch per repeat when matching patterns.

        teloscope asm.fa -x 1

* Example: Window metrics with custom window and step sizes.

        teloscope -f asm.fa -o results/ -j 16 -c TTAGGG -p NNNGGG -w 2000 -s 1000 -g -e -r --verbose

* Example: All outputs enabled.

        teloscope -f asm.fa -o results/ -j 16 -c TTAGGG -p TBAGGG,TTRGGG,YTAGGG -w 2000 -s 1000 -d 200 -l 1000 -x 1 -r -g -e -m -i -n -t 50000 --verbose --cmd

> **Why `-c` and `-p` are separate.**
> `-c` sets the reference motif (determines canonical vs. non-canonical classification and `p`/`q` labeling). `-p` sets what to search for. When `-p` is omitted, Teloscope derives search patterns from `-c` automatically. This means you can switch to a plant canonical (`-c CCCTAAA`) without also passing `-p`.

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
| `-c` | `--canonical` | Canonical telomere repeat used as the reference motif | `TTAGGG` |
| `-p` | `--patterns` | Comma-separated list of repeat patterns to search for | Canonical and its reverse complement |
| `-x` | `--edit-distance` | Allow up to this many mismatches per repeat (0, 1, or 2). Generates all substitution variants of each input pattern. | `0` |

Both the provided patterns and their reverse complements are searched automatically.

### Sliding Window

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-w` | `--window` | Window size in bp | `1000` |
| `-s` | `--step` | Step size in bp | `1000` (non-overlapping) |

**Note:** When step equals window size, windows are non-overlapping bins. This produces valid single-value BEDgraph files compatible with IGV and UCSC genome browsers. Use a smaller `-s` for overlapping windows.

### Block Detection

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-k` | `--max-match-distance` | Maximum gap (bp) between two repeat matches before they are treated as separate | `50` |
| `-d` | `--max-block-distance` | Maximum gap (bp) between two blocks before they stop being merged | `200` |
| `-l` | `--min-block-length` | Minimum block length (bp) to keep in the output | `500` |
| `-y` | `--min-block-density` | Minimum fraction of a block covered by repeat matches (0.0 to 1.0) | `0.5` |
| `-t` | `--terminal-limit` | How far from each contig end (bp) a telomere can be and still count as scaffold-terminal | `50000` |

### Output Control

| Flag | Long form | Description | Default |
|------|-----------|-------------|---------|
| `-r` | `--out-win-repeats` | Output per-window repeat density, canonical ratio, and strand ratio | `false` |
| `-g` | `--out-gc` | Output per-window GC content | `false` |
| `-e` | `--out-entropy` | Output per-window Shannon entropy | `false` |
| `-m` | `--out-matches` | Output individual match coordinates for canonical and terminal non-canonical repeats | `false` |
| `-i` | `--out-its` | Output interstitial telomere (ITS) blocks | `false` |
| `-u` | `--ultra-fast` | Only scan near contig ends; skip genome-wide analysis | `true` |
| `-n` | `--manual-curation` | Include contig-terminal telomeres in BED output (not just scaffold-terminal) | `false` |
| `-v` | `--version` | Print current software version | |
| `-h` | `--help` | Print help | |
| | `--verbose` | Verbose output | |
| | `--cmd` | Print command line | |

**Note:** Enabling any genome-wide output flag (`-r`, `-g`, `-e`, `-m`, `-i`) automatically disables ultra-fast mode.

Outputs
------------
Teloscope always produces the following file:

* `terminal_telomeres.bed` Telomere block annotations for the assembly.

    Columns: `chr`, `start`, `end`, `length`, `label`, `fwdCount`, `revCount`, `canonCount`, `nonCanonCount`, `chrSize`, `type`.

    The `label` column assigns each block to an arm: `p` (mostly forward-strand matches), `q` (mostly reverse-strand), or `b` (balanced). The `type` column is either `scaffold` (block is near a scaffold end) or `contig` (block is near a contig end but not a scaffold end). By default only scaffold-terminal telomeres are reported; use `-n`/`--manual-curation` to also include contig-terminal blocks.

Additional optional outputs (disabled in ultra-fast mode):

* `window_repeat_density.bedgraph` Fraction of each window covered by any repeat match, 0 to 1 (`-r`).
* `window_canonical_ratio.bedgraph` Canonical share of repeat density per window, 0 to 1; -1 when no matches (`-r`).
* `window_strand_ratio.bedgraph` Forward-strand share of repeat density per window, 0 to 1; -1 when no matches (`-r`).
* `window_gc.bedgraph` GC content per window (`-g`).
* `window_entropy.bedgraph` Shannon entropy per window (`-e`).
* `canonical_matches.bed` Coordinates of every canonical repeat match in the assembly (`-m`).
* `noncanonical_matches.bed` Coordinates of non-canonical repeat matches in terminal regions (`-m`).
* `interstitial_telomeres.bed` Telomere-like blocks found away from contig ends, i.e., interstitial telomeres or ITSs (`-i`).

All output filenames are prefixed with the input filename (e.g., `asm.fa_terminal_telomeres.bed`).

Console summary
------------
Teloscope prints two tables to stdout:

**Path Summary Report.** One row per sequence. Columns in ultra-fast mode: `pos`, `header`, `telomeres`, `labels`, `gaps`, `type`, `granular`. In full-scan mode, three columns are added: `its`, `canonical`, `windows`.

**Assembly Summary Report.** Aggregated counts including:
* Total paths, gaps, and telomeres
* Telomere length statistics (mean, median, min, max)
* Chromosome counts by telomere number (two, one, zero)
* Chromosome classification counts (see below)

Chromosome classification
------------
Teloscope classifies each chromosome based on its terminal telomere blocks and whether the sequence contains internal gaps (Ns):

| Type | Gapped variant | Description |
|------|---------------|-------------|
| `t2t` | `gapped_t2t` | Both a p-arm and a q-arm telomere detected, in the expected orientation |
| `incomplete` | `gapped_incomplete` | Only one arm (p or q) has a scaffold-terminal telomere |
| `none` | `gapped_none` | No scaffold-terminal telomeres detected |
| `misassembly` | `gapped_misassembly` | Two telomeres detected but in the wrong order, or the same arm appears twice |
| `discordant` | `gapped_discordant` | A p-arm telomere positioned on the q-side of the chromosome, or vice versa. Not biological; cannot be resolved by curation |

The `granular` column in the path summary shows the full block-by-block label for each sequence. Uppercase marks the longest block on each arm. For example, `Pq` means there is one large p-arm block (longest, uppercase) and one smaller q-arm block. An asterisk (`*`) after a label marks a block whose strand identity (p or q) does not match its position on the chromosome (positional discordance).

How it works
------------
Briefly, **Teloscope** reads an assembly and decomposes its parts. Input patterns (with optional IUPAC codes and edit distance) are expanded into all concrete variants along with their reverse complements, and a multi-pattern search tree is built for fast matching. Each contig is scanned: in ultra-fast mode (default), only the regions near each contig end are scanned; otherwise, the full sequence is analyzed.

Nearby repeat matches are merged into groups, and nearby groups are further merged into telomere blocks. Short blocks and blocks with low repeat density are discarded. Each surviving block is labeled by strand orientation (`p`, `q`, or `b`), and the set of blocks for a chromosome determines its classification (T2T, incomplete, misassembled, discordant, or none). Terminal telomere BED annotations are always written. Optional window metrics, match coordinates, and ITS blocks are produced when the corresponding flags are enabled.

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

Each synthetic sequence (~50 kbp) is built from five segments:

1. **p-arm canonical** 12 kbp or 2000 exact copies of CCCTAA.
2. **p-arm TVR** 600 bp or 100 copies of CCCTAA, each with one random substitution.
3. **Central segment** 25 kbp of random nucleotides.
4. **q-arm TVR** 600 bp or 100 copies of TTAGGG, each with one random substitution.
5. **q-arm canonical** 12 kbp or 2000 exact copies of TTAGGG.

After construction, per-base mutations are applied across the entire sequence at the specified rate (`-r`). The TVR regions have exactly one mismatch per repeat, so they are detectable at `-x 1` but not at `-x 0`. This tests both default and edit-distance modes.

Outputs: `sequences.fa` (synthetic FASTA) and `ground_truth.tsv` (per-contig telomere coordinates with exact boundaries for each segment).

**Evaluate mode** compares teloscope BED output against the ground truth:
```sh
teloscope-simulate --evaluate -g ground_truth.tsv -b terminal_telomeres.bed
```

Output columns: `total_tips`, `detected`, `sensitivity`, `mean_bias_bp`, `mean_abs_err_bp`, `tvr_rate`, `fp_blocks`.

### Simulation parameters

| Flag | Description | Default |
|------|-------------|---------|
| `-n` | Number of sequences to generate | `1000` |
| `-R` | Canonical repeats per tip (each repeat is 6 bp) | `2000` |
| `-T` | TVR (variant) repeats per tip (each repeat is 6 bp with 1 substitution) | `100` |
| `-L` | Central segment length in bp | `25000` |
| `-r` | Per-base mutation rate applied to the full sequence | `0.0` |
| `-s` | Random seed for reproducibility | `42` |
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
