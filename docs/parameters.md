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

There is no flag for FASTQ or BAM input: Teloscope sniffs the content (the BAM magic bytes, or a FASTQ `@` header) and always writes the telomeric reads, a per-read telomere BED, and a report. See [Reads mode](#reads-mode).

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

In reads mode, `-l` (default `300`, same as assembly) sets the measured BED rows and report only, adding reads to the subset only when set below the fixed 42 bp floor. See [Reads mode](#reads-mode) for how that relates to the kept FASTQ/BAM subset.

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

Reads mode, from a file or piped in:

```sh
teloscope reads.fq.gz -j 32 -o results/
teloscope reads.bam -j 32 -o results/
```

## Reads mode

FASTQ or BAM input is detected by content, not by a flag. A regular file is sniffed for the BAM magic bytes or a FASTQ `@` header; gzipped FASTQ files are sniffed the same way, since decompression happens first. Stdin is only peeked one byte, so it recognizes plain FASTQ (`@`) or BAM (BGZF's `0x1f`); gzipped FASTQ on stdin looks like BAM and fails with a hint to pass the file path instead. Anything else is read as an assembly. The two removed flags (`--fastq-subset`, `--bam-subset`) exit 1 with a message pointing at this.

One pass measures and subsets every read together: each read is scanned like an assembly contig end (the same block rule, and the same `-d`, `-y`, `-x`, `-c`, and `-p`), which produces the BED row(s) and report; separately, a read is kept in the output FASTQ/BAM when it has a BED row or when a much more permissive scan — the whole read, a fixed 42 bp floor (seven repeats of the default 6 bp motif) instead of `-l` — finds a block on its own. `-l` (default `300`, same as assembly) sets the measurement; the subset already holds every read that passes the fixed 42 bp scan, so `-l` can only add reads to it through new BED rows.

Two defaults differ from assembly mode unless set explicitly. `--terminal-tolerance` is `300` bp: how far from a read tip a telomere may start. `-t` is `2000` bp: the first tile scanned at each end, raised to `--terminal-tolerance` if smaller. The tile grows while the telomere continues, so it never shortens one.

A telomere is "complete" in the report when its strand matches its end (C-rich at the read start, G-rich at the read end) and the read continues at least `-d` past it; "reaching read end" when the strand matches but the read ends within `-d`; "discordant" when the strand does not match the end. Every row is written to the BED regardless of which of the three it is.

This is an alignment-free estimate. Reads pool chromosome ends by coverage, and a read broken inside a telomere looks complete and pulls the estimate down. Report it with its read counts, not as a replacement for alignment-based tools such as Telogator2.

Assembly output flags are ignored with a warning; assembly record filters (`--include-bed`/`--exclude-bed`/`--include-prefix`/`--exclude-prefix`/`--chr-only`) are rejected outright.

BAM-specific rules: the keep scan runs on every record, including secondary and supplementary ones, so they can still end up in the kept BAM. The measured scan skips secondary/supplementary records (flag `0x900`) and any record hard-clipped at either end of its CIGAR; both counts are printed on stderr. A `0x10` (reverse-strand) record's stored `SEQ` is reverse-complemented before measuring, so an aligned BAM's row matches what the same read's FASTQ row would show.

Outputs are written under `-o` (default: the input's own directory; `./stdin_*` for stdin): `<name>_telomeric.fastq` or `<stem>_telomeric.bam` (the BAM name drops the input's extension), `<name>_terminal_telomeres.bed`, and `<name>_report.tsv`; the report's rows are also printed to stdout. See [Outputs](outputs.md#reads-mode-outputs). On any failure, all three are removed.

## Stdin

Teloscope reads from stdin when no input file is given:

```sh
cat asm.fa | teloscope -o results/
zcat asm.fa.gz | teloscope -o results/
zcat reads.fq.gz | teloscope -o results/
cat reads.bam | teloscope -o results/
```

Compressed FASTA or FASTQ stdin is not supported (gzipped FASTQ on stdin fails with a hint to pass the file path instead). BAM stdin is supported because BAM mode handles BGZF directly.
