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
| `--include-bed FILE` | BED or one-ID-per-line file | analyze whole records whose exact primary IDs occur in column 1 | unset |
| `--exclude-bed FILE` | BED or one-ID-per-line file | exclude whole records whose exact primary IDs occur in column 1 | unset |
| `--include-prefix LIST` | comma-separated literal prefixes | analyze IDs that start with any listed prefix | unset |
| `--exclude-prefix LIST` | comma-separated literal prefixes | exclude IDs that start with any listed prefix | unset |

Exact files keep metadata-derived selections explicit and auditable. Prefix lists cover known accession families without requiring a generated file or regular-expression syntax. These options are long-only because `-i` and `-e` already control ITS and entropy output.

Filtering is off when all four options are unset. With no include option, every record starts selected. Repeated include files and prefix lists form one union. All exclusion matches are removed afterward, so exclusion wins. Exact IDs and prefixes are case-sensitive. Prefixes are literal strings, not globs or regular expressions; use `sample_`, not `sample_*`.

For FASTA, matching uses the first whitespace-delimited token after `>`, including an accession version such as `.11`. For GFA1 `P` records, matching uses path names. A pathless GFA1 graph uses segment names.

Selector files may mix one-column ID rows and BED3+ rows. Blank lines, `#` comments, `track` lines, and `browser` lines are skipped. BED start and end must be unsigned integers with start no greater than end. Coordinates are validated but never crop, join, or remap a sequence; selection always applies to the complete record.

Every exact ID and every prefix must match at least one input-domain ID. Teloscope exits on an unmatched selector, duplicate FASTA primary ID, invalid or empty selector file, invalid BED coordinates, or an empty final selection. It reports the selected and input counts to stderr and in the FASTA summary.

Filtering a GFA does not delete supported original graph records. It scans terminal segment ends reached from selected paths, or selected segments in a pathless graph. Caps and links are graph-level annotations, so an annotated segment end shared by selected and excluded paths is visible to both. Filtered GFA2, GFA1 `C` containment and `W` walk records, and unknown GFA record types are rejected because the current writer cannot preserve them safely.

Filtering happens after the assembly is loaded. It reduces scanning and output size, including `--plot-report`, but not input parsing or peak loader memory. Filters require FASTA or supported GFA1 input; uncompressed filtered stdin is treated as FASTA, and filtered GFA input must use a `.gfa` or `.gfa.gz` filename. Filtered FASTQ/BAM input and both read-subset modes are rejected.

Database prefixes are format-specific. `GCA_` and `GCF_` are NCBI assembly accessions, not sequence IDs. `CM` is a GenBank/INSDC scaffold/CON family rather than a universal chromosome label, and `NC_` works only when that text starts the actual FASTA primary ID. For chromosome-only selection, prefer exact accession.version IDs derived from the NCBI sequence or assembly report. See [Database FASTA names](../README.md#database-fasta-names).

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

One inspected submitter naming family:

```sh
teloscope asm.fa --include-prefix hap1_chr -o results/
```

Exact chromosome IDs derived from assembly metadata:

```sh
teloscope asm.fa --include-bed chromosomes.ids -o results/
```

Exact include list with a final exclusion list:

```sh
teloscope asm.fa --include-bed primary.ids --exclude-bed do_not_plot.ids -o results/
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
