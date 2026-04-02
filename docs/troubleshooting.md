[Back to README](../README.md)

# Troubleshooting

Use this page when a run fails, writes less than expected, or produces calls that do not match what you expected.

Start with:

```sh
teloscope asm.fa --cmd --verbose
```

`--cmd` shows the resolved command line. `--verbose` makes it easier to see whether the run is failing during input parsing, scanning, or report generation.

## Setup problems

### Build fails because `gfalibs` is missing

Teloscope needs the `gfalibs` submodule for GFA support.

```sh
git submodule update --init --recursive
make -j
```

### `--plot-report` fails immediately

The report step needs Python 3 plus:

- `matplotlib`
- `numpy`
- `pandas`

If those packages are missing, run Teloscope without `--plot-report` or install them in the Python environment that will run `scripts/teloscope_report.py`.

## Input and argument errors

### No input file was provided

If you see `Error: Input sequence file is required` or `Error: No input file provided`, pass the input as either:

```sh
teloscope asm.fa
teloscope -f asm.fa
```

### A flag is missing its value

If you see `Error: Option -<flag> is missing a required argument`, one of the value-taking flags was given without a value.

Common cases:

- `-f`
- `-o`
- `-c`
- `-p`
- `-j`
- `-w`
- `-s`
- `-t`
- `-k`
- `-d`
- `-l`
- `-y`
- `-x`

### Compressed stdin is not supported

This does not work:

```sh
cat asm.fa.gz | teloscope -o results/
```

Use one of these:

```sh
teloscope asm.fa.gz
zcat asm.fa.gz | teloscope -o results/
```

### Output directory is not writable

If you see `Cannot create output directory`, `Output directory ... is not writable`, or `Could not open ... for writing`, check:

- the directory exists or can be created
- you have write permission
- the filesystem has free space
- no other process is writing the same outputs at the same time

Quick check:

```sh
mkdir -p results
test -w results && echo ok
```

## Invalid flag values

These flags must be positive integers:

- `-w`
- `-s`
- `-t`
- `-k`
- `-d`
- `-l`

`-y` must be a number from `0` to `1`.

`-x` must be `0`, `1`, or `2`.

If Teloscope says a value is invalid, start from a known-good baseline:

```sh
teloscope asm.fa -w 1000 -s 1000 -t 50000 -k 50 -d 200 -l 500 -y 0.5 -x 1
```

Also keep `-s <= -w`. A larger step than window size is rejected.

## Missing or unexpected outputs

### FASTA run wrote fewer files than expected

By default, FASTA mode writes:

- `*_terminal_telomeres.bed`
- `*_gaps.bed`
- `*_report.tsv`

Everything else depends on flags:

- `-r`: repeat-density, canonical-ratio, strand-ratio BEDgraphs
- `-g`: GC BEDgraph
- `-e`: entropy BEDgraph
- `-m`: canonical and non-canonical match BED files
- `-i`: interstitial telomere BED
- `--plot-report`: PDF report

### GFA run did not write BED or TSV files

That is expected. GFA mode writes only:

```text
<input>.telo.annotated.gfa
```

### `-n` did not change classification

`-n/--manual-curation` only keeps contig-terminal telomeres in the BED output. It does not change the scaffold classification logic in `*_report.tsv`.

### Runtime increased after enabling output flags

`-r`, `-g`, `-e`, `-m`, and `-i` disable ultra-fast end-only scanning and force full-sequence analysis. That slowdown is expected.

## Calls look wrong

### No telomeres were called

Check these first:

- `-c` matches the organism
- `-t` is large enough to include the terminal block
- `-l` is not too strict
- `-y` is not too strict
- `-x 0` is not too strict for the assembly

A useful debugging run is:

```sh
teloscope asm.fa -t 100000 -l 200 -y 0.3 -x 1 --verbose
```

Then tighten one parameter at a time.

### Telomeres are split into too many small blocks

`-k` and `-d` are probably too small.

- `-k` controls how far apart repeat matches can be before they stop merging
- `-d` controls how far apart nearby repeat groups can be before block extension stops

Increase them gradually if blocks are fragmenting.

### Too much sequence merged into one block

`-k` or `-d` is probably too large. Lower them if unrelated repeat runs are being fused.

### Classification looks wrong

Check:

- the canonical motif in `-c`
- whether the telomere block is still within `-t`
- whether you are comparing scaffold-terminal logic to contig-terminal BED output

Remember:

- `*_report.tsv` classifies from scaffold-terminal blocks
- `-n` affects BED retention, not classification

### Pattern matching became very slow

This usually means the search set exploded because of:

- many entries in `-p`
- IUPAC ambiguity codes
- a larger `-x`

If Teloscope warns that pattern count is unusually high, reduce ambiguity in `-p` or lower `-x`.

### `p` and `q` labels look reversed for the organism

`-c` defines the canonical repeat used for canonical counting and strand labeling. If `-c` is wrong, the labels can look wrong even if repeat matching itself succeeded.

## GFA-specific behavior

### The annotated GFA has no telomere nodes

The usual causes are:

- segment sequence is `*`, so there is nothing to scan
- the detected block fails `-l` or `-y`
- the segment end is not path-terminal in a graph with paths
- `-c` does not match the assembly motif

### Only some segment ends were annotated

When paths are present, Teloscope annotates path-terminal segment ends. It does not annotate every segment end in the graph.

### You want to confirm that telomere nodes were added

Check the annotated graph directly:

```sh
rg "telomere_" results/asm.gfa.telo.annotated.gfa
```

## Report generation

### Teloscope could not find `teloscope_report.py`

If you see `Warning: Could not locate teloscope_report.py` or `Warning: Report generation failed`, run the script directly:

```sh
python3 scripts/teloscope_report.py results/
```

### The standalone report script says files are missing

Point it at a Teloscope output directory, not the repo root:

```sh
python3 scripts/teloscope_report.py results/
```

At minimum, that directory must contain `*_terminal_telomeres.bed`.

## Quick checklist

If a run looks wrong, check these in order:

1. Confirm whether you are running FASTA mode or GFA mode.
2. Confirm the canonical motif in `-c`.
3. Confirm which optional outputs you actually requested.
4. Check whether `-r`, `-g`, `-e`, `-m`, or `-i` forced full-sequence scanning.
5. Revisit `-t`, `-l`, `-y`, `-k`, and `-d`.
6. Rerun once with `--cmd --verbose`.
