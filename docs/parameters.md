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
| `-d` | `--max-block-distance` | max gap between nearby repeat groups before extension stops | `200` |
| `-l` | `--min-block-length` | minimum block length to keep | `500` for assembly, `60` for read subsets |
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

In read subset modes, the default `-l` is `60` bp. This is intended to retain reads with at least about ten telomeric repeat units after Teloscope's block and density filters. Assembly annotation keeps the stricter `500` bp default.

## Stdin

Teloscope reads from stdin when no input file is given:

```sh
cat asm.fa | teloscope -o results/
zcat asm.fa.gz | teloscope -o results/
zcat reads.fq.gz | teloscope --fastq-subset > telomeric.fq
cat reads.bam | teloscope --bam-subset > telomeric.bam
```

Compressed FASTA or FASTQ stdin is not supported. BAM stdin is supported because BAM mode handles BGZF directly.
