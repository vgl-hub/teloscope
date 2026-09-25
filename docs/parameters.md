[Back to README](index.md)

# Parameters

For the full built-in help:

```sh
teloscope -h
```

Basic usage:

```sh
teloscope input.fa [options]
teloscope input.fa.gz [options]
teloscope input.gfa [options]
teloscope -f input.fa [options]
teloscope input.fq.gz [options]
teloscope input.bam [options]
```

## Input and output

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-f` | `--input-sequence` | input FASTA, FASTA.gz, GFA, FASTQ, or BAM file | required unless passed positionally |
| `-o` | `--output` | output directory | input file directory (`.` for stdin) |
| `-j` | `--threads` | maximum worker threads | all available |

## Assembly record filters

| Long form | Value | Meaning | Default |
| --- | --- | --- | --- |
| `--include-bed` | `FILE` | keep IDs listed in column 1 | unset |
| `--exclude-bed` | `FILE` | remove IDs listed in column 1 | unset |
| `--include-prefix` | `LIST` | keep IDs with any comma-separated prefix | unset |
| `--exclude-prefix` | `LIST` | remove IDs with any comma-separated prefix | unset |
| `--chr-only` | — | keep records named like the longest one | `false` |

Filtering is off when all flags are unset. The include/exclude flags can be repeated. `--chr-only` is a single on/off switch. Includes form a union and exclusions run last; without includes, all records start selected. Matching is case-sensitive, and prefixes are literal strings rather than globs or regular expressions.

`--chr-only` takes the longest record as a chromosome and reads its name as the assembly's naming convention. A record is kept when its name starts with the longest record's leading run of letters and has the same number of separator characters (characters that are neither letters nor digits). GenBank `CM093074.1` keeps `CM…` and drops `JAXXXX010000001.1`. RefSeq `NC_…` drops `NW_…`. ENA `OZ124247.1` drops `CAUPLK010000001.1`. `SUPER_1` keeps `SUPER_Z` and `SUPER_1A` and drops `SUPER_1_unloc_3` and `SCAFFOLD_12`. `chr1` keeps `chrX` and `chrM` and drops `chrUn_KI270302v1` and `chr1_KI270706v1_random`. `Chr01` drops `scaffold123`. The mitochondrion is kept when it follows the convention (for example `CM010492.2`). In a bare-number scheme (`1`…`22`, `X`) the letter run is empty, so only the separator count decides. When every record follows one convention, everything is kept. stderr prints the longest record, its prefix and separator count, and how many records were selected. `--include-bed` and `--exclude-bed` combine with it as usual: includes form a union and exclusions run last.

FASTA matching uses the first token after `>`, including accession versions. GFA1 matching uses `P` path names or `S` segment names when no paths exist.

Selector files accept one-ID rows or BED3+ rows and skip blank, `#`, `track`, and `browser` lines. BED coordinates must be unsigned integers with start no greater than end, but they never crop or remap a sequence.

Every exact ID and prefix must match at least one input name. Teloscope rejects unmatched selectors, duplicate FASTA IDs, invalid or empty selector files, invalid BED coordinates, and empty final selections. It reports selected and input counts to stderr and in the FASTA summary.

Filters require FASTA or supported GFA1 input. GFA filtering preserves supported graph records and only limits scanned terminal ends. Filtered GFA accepts `H`, `S`, `L`, `J`, and `P` records; it rejects GFA2, `C`, `W`, and unknown records. Filtered stdin is parsed as FASTA, so GFA input needs a `.gfa` or `.gfa.gz` filename. FASTQ and BAM input reject filters.

Filtering occurs after input loading. It reduces scanning and output size, including `--plot-report`, but not parsing or peak memory. FASTA BED, BEDgraph, TSV, and report outputs contain selected records only.

For database FASTA, use exact accession.version IDs from the NCBI [genome sequence report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/) or [assembly report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/). `GCA_` and `GCF_` are assembly accessions, and prefixes such as `CM` or `NC_` are not universal chromosome rules. Use prefixes only after checking the input headers.

## Pattern control

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-c` | `--canonical` | reference telomere repeat | `TTAGGG` |
| `-p` | `--patterns` | comma-separated search patterns | derived from `-c` |
| `-x` | `--edit-distance` | allowed mismatches per repeat unit (`0`–`2`) | `1` |

Notes:

- Reverse complements are always searched automatically.
- IUPAC ambiguity codes are allowed in `-p`.
- `-c` controls canonical versus non-canonical counting even when `-p` is set explicitly.

## Windowing

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-w` | `--window` | window size in bp | `1000` |
| `-s` | `--step` | step size in bp | `1000` |

When `-s` equals `-w`, window outputs are non-overlapping BEDgraph bins.

## Block calling

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-k` | `--max-match-distance` | matches this close chain into one interstitial seed | `50` |
| `-d` | `--max-block-distance` | maximum non-telomeric stretch inside a telomere | `1000` |
| `-l` | `--min-block-length` | minimum piece length to keep | `300` |
| `-y` | `--min-block-density` | minimum repeat-covered fraction for a piece, in `(0,1]` | `0.5` |
| `-t` | `--terminal-limit` | how far in from each end to look | `50000` |
|  | `--terminal-tolerance` | how far from an end a telomere may start | `3000` |
|  | `--label-threshold` | forward-strand fraction for the p/q label; `b` between, in `(0.5,1]` | `0.667` |
|  | `--min-block-counts` | minimum canonical matches per block | `2` |

The start zone is the smaller of `--terminal-tolerance` and `-t`, counted in called bases.

A piece or row exactly at a threshold — `-l`, `-y`, or `--min-block-counts` — is kept, not rejected.

## Output flags

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-r` | `--out-win-repeats` | write repeat density, canonical ratio, and strand ratio tracks | `false` |
| `-g` | `--out-gc` | write GC BEDgraph | `false` |
| `-e` | `--out-entropy` | write entropy BEDgraph | `false` |
| `-m` | `--out-matches` | write canonical and terminal non-canonical match BED files | `false` |
| `-i` | `--out-its` | scan whole sequences for interstitial telomeres (the interstitial BED is always written) | `false` |
| `-u` | `--ultra-fast` | scan sequence ends only | `true` |
| `-n` | `--manual-curation` | also report telomeres at contig ends | `false` |
| `-a` | `--out-fasta` | write the terminal telomere sequences as FASTA | `false` |
|  | `--plot-report` | write separate terminal and ITS PDF reports after the run | `false` |

Any of `-r`, `-g`, `-e`, `-m`, or `-i` forces the full scan. `-n` keeps the fast scan but reads both end windows of every contig and adds contig-terminal rows to the terminal BED. These are assembly-only outputs: with FASTQ or BAM input they are ignored, with a warning.

## Informational flags

| Flag | Long form | Meaning |
| --- | --- | --- |
| `-v` | `--version` | print the program version |
| `-h` | `--help` | print help text |
|  | `--verbose` | print extra progress messages |
|  | `--cmd` | print the resolved command line |

## Common command patterns

Default vertebrate run:

```sh
teloscope asm.fa
```

Filter assembly records:

```sh
teloscope asm.fa --include-bed primary.ids --exclude-bed do_not_plot.ids -o results/
teloscope asm.fa --include-prefix hap1_chr,hap2_chr -o results/
```

Plant canonical repeat:

```sh
teloscope asm.fa -c CCCTAAA
```

Explicit search patterns:

```sh
teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG
```

All optional FASTA outputs plus the report:

```sh
teloscope asm.fa -o results/ -r -g -e -m -i --plot-report
```

Graph annotation:

```sh
teloscope asm.gfa -o results/
```

Reads mode:

```sh
teloscope reads.fq.gz -j 32 -o results/
teloscope reads.bam -j 32 -o results/
```

## Reads mode

FASTQ or BAM input is detected from its first bytes; there is no flag, and `--fastq-subset`/`--bam-subset` were removed. `-l` sets the measured rows; reads are kept with a fixed 42 bp floor, so `-l` below 42 also keeps more reads. See [Outputs](outputs.md#reads-mode-outputs).

## Stdin and pipes

Teloscope reads from stdin when no input file is given. A path that is not a regular file (a named pipe, `<(...)`, `/dev/stdin`) is read the same way, and its outputs default to `.`.

```sh
cat asm.fa | teloscope -o results/
zcat asm.fa.gz | teloscope -o results/
zcat reads.fq.gz | teloscope -o results/
cat reads.bam | teloscope -o results/
```

Compressed FASTA or FASTQ on stdin or a pipe is not supported; BAM is.
