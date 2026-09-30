[Back to README](index.md)

# Troubleshooting

Rerun with `--verbose` first: it prints progress, which shows where a run stops.

## Errors

| Message | Fix |
| --- | --- |
| `file does not exist`, `No input file provided` | pass an existing input as the first argument or with `-f` |
| `Option -<flag> is missing a required argument` | give the flag a value |
| `input ... is not FASTA, GFA, FASTQ or BAM` | only the first bytes are checked, so a malformed GFA can still fail later |
| `Cannot create output directory`, `Output directory ... is not writable` | check the path, the permissions, and the free space |
| `Step size (...) cannot be larger than window size (...)` | keep `-s` at or below `-w` |
| `... must be > 0`, `... must be in the range ...` | `-w`, `-s`, `-t`, `-k`, `-d`, `-l`, `--terminal-tolerance`, and `--min-block-counts` take positive integers; `-y` takes `(0,1]`, `-x` `0` to `2`, and `--label-threshold` `(0.5,1]` |
| `... patterns is unusually high` | over 500 patterns: use fewer IUPAC codes in `-p` or a lower `-x` |
| `Could not locate teloscope_report.py` | see [Report](#report) |

## Build

If `gfalibs` is missing:

```sh
git submodule update --init --recursive
make -j
```

## Input

- Compressed input on a pipe fails, except BAM. Pass the file, or decompress into the pipe:

  ```sh
  teloscope asm.fa.gz                      # works
  zcat asm.fa.gz | teloscope -o results/   # works
  cat asm.fa.gz | teloscope -o results/    # fails
  ```

- Gzipped FASTQ on a pipe is read as BAM and fails with a hint; pass the file.
- BAM must be BGZF-compressed; SAM, CRAM, plain gzip, and uncompressed BAM are rejected. A missing BGZF end-of-file marker only warns.

## Outputs

- A FASTA run writes four files by default; the rest need flags ([Outputs](outputs.md)).
- A GFA run writes only the annotated graph and its color file.
- `-n` adds `contig` rows to the terminal BED (and to `-a`) and changes which arrays count as interstitial: a full scan loses them, and a fast scan can gain rows from the extra windows. `type`, `anomaly`, and the statistics stay the same.
- `-r`, `-g`, `-e`, `-m`, and `-i` force the full scan, and `-n` reads every contig end, so all of them cost time.
- `--chr-only` misses a chromosome whose name breaks the longest record's convention, such as `chromosome_X` beside `chr1`. Add it with `--include-bed`.

## Calls

Terminal calls depend most on `-c`, `--terminal-tolerance`, `-l`, `-y`, `-d`, and `-x`; interstitial rows on `-k`, `-d`, and `-y`.

- **No telomeres:** check that `-c` matches the organism and that the start zone (`-t`, `--terminal-tolerance`) reaches the telomere. Try a permissive run, then restore one threshold at a time:

  ```sh
  teloscope asm.fa -t 100000 --terminal-tolerance 100000 -l 200 -y 0.3 --verbose
  ```

- **A telomere in pieces:** raise `-d`, the longest stretch a telomere may bridge.
- **Wrong class:** `type` and `anomaly` come from the two scaffold arms only. Recheck `-c`, `-t`, and `--terminal-tolerance`, and don't compare them with `contig` rows.
- **Wrong `p`/`q` labels:** `-c` sets the canonical motif and with it the strand labels.

## GFA

- **No caps:** the segment sequence is `*`, the block fails `-l` or `-y`, the end is not path-terminal, or `-c` does not match.
- **Caps on some ends only:** with paths, only path-terminal segment ends are scanned.
- **Checking for caps:** `grep telomere_ results/asm.gfa.telo.annotated.gfa`

## Report

- `--plot-report` needs Python 3 with `matplotlib`, `numpy`, and `pandas`.
- On `Could not locate teloscope_report.py` or `Report generation failed`, run the script from a source checkout on the output directory:

  ```sh
  python3 scripts/teloscope_report.py results/
  ```

- The script needs `*_terminal_telomeres.bed` or `*_interstitial_telomeres.bed` in that directory, so point it at the outputs, not the repo root.
