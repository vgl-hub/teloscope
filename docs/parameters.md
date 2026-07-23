[Back to README](../README.md)

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
teloscope --fastq-subset input.fq.gz [options] > telomeric.fq
teloscope --bam-subset input.bam [options] > telomeric.bam
```

## Input and output

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-f` | `--input-sequence` | input FASTA, FASTA.gz, GFA, FASTQ, or BAM file | required unless passed positionally |
| `-o` | `--output` | output directory | input file directory |
| `-j` | `--threads` | maximum worker threads | all available |
|  | `--fastq-subset` | stream FASTQ reads with Teloscope-valid telomeric blocks to stdout, or to a file with `-o` | `false` |
|  | `--bam-subset` | stream BAM records with Teloscope-valid telomeric blocks to stdout, or to a file with `-o` | `false` |

## Assembly record filters

| Long form | Value | Meaning | Default |
| --- | --- | --- | --- |
| `--include-bed` | `FILE` | keep IDs listed in column 1 | unset |
| `--exclude-bed` | `FILE` | remove IDs listed in column 1 | unset |
| `--include-prefix` | `LIST` | keep IDs with any comma-separated prefix | unset |
| `--exclude-prefix` | `LIST` | remove IDs with any comma-separated prefix | unset |

Filtering is off when all four flags are unset. Flags can be repeated. Includes form a union and exclusions run last; without includes, all records start selected. Matching is case-sensitive, and prefixes are literal strings rather than globs or regular expressions.

FASTA matching uses the first token after `>`, including accession versions. GFA1 matching uses `P` path names or `S` segment names when no paths exist.

Selector files accept one-ID rows or BED3+ rows and skip blank, `#`, `track`, and `browser` lines. BED coordinates must be unsigned integers with start no greater than end, but they never crop or remap a sequence.

Every exact ID and prefix must match at least one input name. Teloscope rejects unmatched selectors, duplicate FASTA IDs, invalid or empty selector files, invalid BED coordinates, and empty final selections. It reports selected and input counts to stderr and in the FASTA summary.

Filters require FASTA or supported GFA1 input. GFA filtering preserves supported graph records and only limits scanned terminal ends. Filtered GFA accepts `H`, `S`, `L`, `J`, and `P` records; it rejects GFA2, `C`, `W`, and unknown records. Filtered stdin is parsed as FASTA, so GFA input needs a `.gfa` or `.gfa.gz` filename. FASTQ, BAM, and both read-subset modes reject filters.

Filtering occurs after input loading. It reduces scanning and output size, including `--plot-report`, but not parsing or peak memory. FASTA BED, BEDgraph, TSV, and report outputs contain selected records only.

For database FASTA, use exact accession.version IDs from the NCBI [genome sequence report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/) or [assembly report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/). `GCA_` and `GCF_` are assembly accessions, and prefixes such as `CM` or `NC_` are not universal chromosome rules. Use prefixes only after checking the input headers.

## Pattern control

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-c` | `--canonical` | reference telomere repeat | `TTAGGG` |
| `-p` | `--patterns` | comma-separated search patterns | derived from `-c` |
| `-x` | `--edit-distance` | allowed mismatches per repeat unit | `1` |

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
| `-k` | `--max-match-distance` | max gap between matches before splitting them | `50` |
| `-d` | `--max-block-distance` | max gap between nearby repeat groups before extension stops | `500` |
| `-l` | `--min-block-length` | minimum block length to keep | `300` for assembly, `42` for read subsets |
| `-y` | `--min-block-density` | minimum repeat-covered fraction for a block | `0.5` |
| `-t` | `--terminal-limit` | distance from a sequence end that still counts as terminal | `50000` |

## Output flags

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-r` | `--out-win-repeats` | write repeat density, canonical ratio, and strand ratio tracks | `false` |
| `-g` | `--out-gc` | write GC BEDgraph | `false` |
| `-e` | `--out-entropy` | write entropy BEDgraph | `false` |
| `-m` | `--out-matches` | write canonical and terminal non-canonical match BED files | `false` |
| `-i` | `--out-its` | write interstitial telomere BED | `false` |
| `-u` | `--ultra-fast` | scan sequence ends only | `true` |
| `-n` | `--manual-curation` | include contig-terminal blocks in BED output | `false` |
|  | `--plot-report` | write a PDF report after the run | `false` |

Any of `-r`, `-g`, `-e`, `-m`, or `-i` disables ultra-fast mode automatically.

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

FASTQ read subset before mapping:

```sh
teloscope --fastq-subset reads.fq.gz -j 32 | minimap2 -ax map-hifi ref.fa -
```

BAM record subset:

```sh
teloscope --bam-subset reads.bam -j 32 > telomeric.bam
```

In read subset modes, the default `-l` is `42` bp. This is intended to retain reads with at least about seven telomeric repeat units after Teloscope's block and density filters. Assembly annotation keeps the stricter `300` bp default.

## Stdin

Teloscope reads from stdin when no input file is given:

```sh
cat asm.fa | teloscope -o results/
zcat asm.fa.gz | teloscope -o results/
zcat reads.fq.gz | teloscope --fastq-subset > telomeric.fq
cat reads.bam | teloscope --bam-subset > telomeric.bam
```

Compressed FASTA or FASTQ stdin is not supported. BAM stdin is supported because BAM mode handles BGZF directly.
