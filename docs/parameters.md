[← Back to README](../README.md)

Parameters
============

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
| `-x` | `--edit-distance` | Allow up to this many mismatches per repeat (0, 1, or 2). Generates all substitution variants of each input pattern. | `1` |

Both the provided patterns and their reverse complements are searched automatically.

> **Why `-c` and `-p` are separate.**
> `-c` sets the reference motif (determines canonical vs. non-canonical classification and `p`/`q` labeling). `-p` sets what to search for. When `-p` is omitted, Teloscope derives search patterns from `-c` automatically. This means you can switch to a plant canonical (`-c CCCTAAA`) without also passing `-p`.

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
| | `--report` | Generate a PDF report after analysis (requires Python 3 + matplotlib) | `false` |
| `-v` | `--version` | Print current software version | |
| `-h` | `--help` | Print help | |
| | `--verbose` | Verbose output | |
| | `--cmd` | Print command line | |

**Note:** Enabling any genome-wide output flag (`-r`, `-g`, `-e`, `-m`, `-i`) automatically disables ultra-fast mode.

### Piping from stdin

Teloscope accepts input from stdin when no file is given. Gzipped input must be decompressed first (`zcat`).

```sh
cat asm.fa | teloscope -o results/
zcat asm.fa.gz | teloscope -o results/
```

Download and analyze directly from NCBI:
```sh
datasets download genome accession GCF_000001405.40 --include genome --filename ncbi.zip
unzip -p ncbi.zip '*.fna' | teloscope -o results/
```
