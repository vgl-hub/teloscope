[Back to README](index.md)

# Parameters

```sh
teloscope <input> [options]
```

The input is FASTA, FASTA.gz, GFA, FASTQ, FASTQ.gz, or BAM, given as the first argument or with `-f`. Without one, Teloscope reads stdin. `teloscope -h` prints every flag.

## Input and output

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-f` | `--input-sequence` | input file | first argument |
| `-o` | `--output` | output directory | input file directory, `.` for stdin |
| `-j` | `--threads` | maximum worker threads | all available |

## Patterns

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-c` | `--canonical` | canonical telomere repeat | `TTAGGG` |
| `-p` | `--patterns` | comma-separated search patterns, IUPAC codes allowed | derived from `-c` |
| `-x` | `--edit-distance` | substitutions allowed per repeat, `0` to `2` | `1` |

Every pattern is also searched as its reverse complement. `-c` decides canonical versus variant counts and the `p`/`q` labels, even when `-p` is set.

## Block calling

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-t` | `--terminal-limit` | how far in from each end to look | `50000` |
|  | `--terminal-tolerance` | how far from an end a telomere may start | `3000` |
| `-k` | `--max-match-distance` | matches this close chain into one interstitial seed | `50` |
| `-d` | `--max-block-distance` | longest non-telomeric stretch inside a telomere | `1000` |
| `-l` | `--min-block-length` | shortest piece kept | `300` |
| `-y` | `--min-block-density` | repeat-covered fraction a piece needs, in `(0,1]` | `0.5` |
|  | `--min-block-counts` | canonical matches a block needs | `2` |
|  | `--label-threshold` | forward-strand share for a `p` or `q` label, in `(0.5,1]` | `0.667` |

The start zone is the smaller of `--terminal-tolerance` and `-t`, counted in called bases. A value exactly at `-l`, `-y`, or `--min-block-counts` passes.

## Windows

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-w` | `--window` | window size in bp | `1000` |
| `-s` | `--step` | step in bp, at most `-w` | `1000` |

## Outputs

| Flag | Long form | Writes | Default |
| --- | --- | --- | --- |
| `-r` | `--out-win-repeats` | repeat density, canonical ratio, and strand ratio per window | off |
| `-g` | `--out-gc` | GC content per window | off |
| `-e` | `--out-entropy` | Shannon entropy per window | off |
| `-m` | `--out-matches` | canonical and terminal non-canonical matches | off |
| `-a` | `--out-fasta` | terminal telomere sequences | off |
| `-i` | `--out-its` | interstitial telomeres from a whole-sequence scan | off |
| `-n` | `--manual-curation` | telomeres at contig ends too | off |
| `-u` | `--ultra-fast` | a scan of sequence ends only | on |
|  | `--plot-report` | terminal and ITS PDF reports | off |

`-r`, `-g`, `-e`, `-m`, or `-i` switch to a full scan. `-n` keeps the fast scan but reads both ends of every contig. [Outputs](outputs.md) lists the files.

## Record filters

| Flag | Effect | Default |
| --- | --- | --- |
| `--include-bed FILE` | keeps IDs in column 1 | unset |
| `--exclude-bed FILE` | removes IDs in column 1 | unset |
| `--include-prefix LIST` | keeps IDs starting with any comma-separated prefix | unset |
| `--exclude-prefix LIST` | removes IDs starting with any comma-separated prefix | unset |
| `--chr-only` | keeps records named like the longest one | off |

- Filtering is off unless a flag is set. Include and exclude flags repeat. Includes, `--chr-only` among them, form a union; exclusions run last.
- Matching is case-sensitive, and prefixes are literal, not globs. FASTA uses the first word after `>`, accession version included. GFA1 uses `P` names, or `S` names in a graph without paths.
- A selector file holds one ID per line or BED3+ rows. Column 1 selects the whole record; coordinates must be valid but never crop it. Blank, `#`, `track`, and `browser` lines are skipped.
- Every ID and prefix must match a record, and the selection cannot be empty. Duplicate FASTA IDs are rejected. The input and selected counts go to stderr and the report.
- Filters work on FASTA and on GFA1 made of `H`, `S`, `L`, `J`, and `P` records; anything else is rejected, as are FASTQ and BAM. A filtered GFA keeps every record and only limits which ends are scanned. Filtered stdin is read as FASTA, so a filtered GFA needs a `.gfa` file name.
- Filters apply after loading, so they cut scan time and output size, not memory.

`--chr-only` takes the longest record as a chromosome. A record is kept when its name starts with the same letters and has the same number of separators (characters that are neither letters nor digits). stderr prints the rule it used.

| Longest record | Kept | Dropped |
| --- | --- | --- |
| `CM093074.1` | other `CM` IDs | `JAXXXX010000001.1` |
| `NC_000001.11` | other `NC_` IDs | `NW_` IDs |
| `OZ124247.1` | other `OZ` IDs | `CAUPLK010000001.1` |
| `SUPER_1` | `SUPER_Z`, `SUPER_1A` | `SUPER_1_unloc_3`, `SCAFFOLD_12` |
| `chr1` | `chrX`, `chrM` | `chrUn_KI270302v1`, `chr1_KI270706v1_random` |

A mitochondrion is kept when it follows the convention (`CM010492.2`). With bare numbers (`1` to `22`, `X`), only the separator count decides. To add a missed chromosome, list it with `--include-bed`.

For NCBI FASTA, use accession.version IDs from the [genome sequence report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/) or [assembly report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/). `GCA_` and `GCF_` name assemblies, not records, and `CM` or `NC_` are not universal chromosome prefixes.

## Reads mode

FASTQ or BAM input is detected from its first bytes. Two defaults change: `-t` becomes `2000`, the read tile, and `--terminal-tolerance` becomes `300`. `-l` sets which telomeres are measured. A read is kept when it has a measured telomere or a telomeric block of at least 42 bp. BAM is read directly, without HTSlib or samtools. Assembly-only output flags are ignored with a warning. [Outputs](outputs.md#reads-mode) describes the files.

## Stdin and pipes

Input that is not a regular file (a pipe, `<(...)`, `/dev/stdin`) is read as a stream, and outputs default to `.`:

```sh
zcat asm.fa.gz | teloscope -o results/
cat reads.bam | teloscope -o results/
```

Compressed FASTA or FASTQ on a stream is not supported; BAM is.

## Other flags

| Flag | Long form | Meaning |
| --- | --- | --- |
| `-v` | `--version` | print the version |
| `-h` | `--help` | print the help |
|  | `--verbose` | print progress |
|  | `--cmd` | print the resolved command line |
