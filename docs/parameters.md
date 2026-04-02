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
```

## Input and output

| Flag | Long form | Meaning | Default |
| --- | --- | --- | --- |
| `-f` | `--input-sequence` | input FASTA, FASTA.gz, or GFA file | required unless passed positionally |
| `-o` | `--output` | output directory | input file directory |
| `-j` | `--threads` | maximum worker threads | all available |

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
| `-l` | `--min-block-length` | minimum block length to keep | `500` |
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

## Stdin

Teloscope reads from stdin when no input file is given:

```sh
cat asm.fa | teloscope -o results/
zcat asm.fa.gz | teloscope -o results/
```

Compressed stdin is not supported. Decompress first.
