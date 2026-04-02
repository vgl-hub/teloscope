[Back to README](../README.md)

# Troubleshooting

This page is organized by the symptom you see at the command line or in the output directory. The examples use `teloscope asm.fa`, but the same checks apply to `asm.fa.gz` and `asm.gfa` unless noted.

## Start with the built-in debug flags

If a run is behaving differently than expected, these two flags help first:

```sh
teloscope asm.fa --cmd --verbose
```

`--cmd` prints the resolved command line. `--verbose` prints more progress messages while the run is active.

## Build and setup problems

### Build fails because `gfalibs` is missing

Teloscope needs the `gfalibs` submodule for GFA support.

Fix:

```sh
git submodule update --init --recursive
make -j
```

## Input and argument errors

### No input file was provided

You may see one of these:

```text
Error: Input sequence file is required. Use -f or --input-sequence.
Error: No input file provided. Use -f or pass as positional argument.
```

Use one of these forms:

```sh
teloscope asm.fa
teloscope -f asm.fa
teloscope asm.gfa -o results/
```

### A flag is missing its value

You may see:

```text
Error: Option -o is missing a required argument
```

This happens when a value-taking flag is followed by nothing or by another flag. Common cases:

- `-f`
- `-o`
- `-c`
- `-p`
- `-j`
- `-t`
- `-k`
- `-d`
- `-l`
- `-y`
- `-x`
- `-w`
- `-s`

Example of a bad command:

```sh
teloscope asm.fa -o
```

Correct form:

```sh
teloscope asm.fa -o results/
```

### Extra positional arguments are being ignored

You may see:

```text
Warning: Ignoring extra positional argument '...'
```

Teloscope accepts one input file. If you need to process multiple assemblies, run Teloscope once per file.

### Compressed stdin is not supported

This fails:

```sh
cat asm.fa.gz | teloscope -o results/
```

Use one of these instead:

```sh
teloscope asm.fa.gz
zcat asm.fa.gz | teloscope -o results/
```

## Output directory and file surprises

### The output directory is not writable

You may see one of these:

```text
Error: Cannot create output directory '...'
Error: Output directory '...' is not writable.
Error: Could not open '...' for writing.
```

Check:

- the directory exists or can be created
- you have write permission
- the filesystem has free space
- another process is not writing the same files at the same time

Quick check:

```sh
mkdir -p results
test -w results && echo ok
```

### You expected optional files, but they were not written

Some outputs are only written when their flags are enabled.

Examples:

- `-r` writes repeat-density, canonical-ratio, and strand-ratio BEDgraphs
- `-g` writes GC BEDgraph
- `-e` writes entropy BEDgraph
- `-m` writes canonical and non-canonical match BED files
- `-i` writes interstitial telomere BED
- `--plot-report` writes the PDF report

If you only run:

```sh
teloscope asm.fa -o results/
```

you should expect the default core outputs only:

- `*_terminal_telomeres.bed`
- `*_gaps.bed`
- `*_report.tsv`

### You expected FASTA-style BED files from a GFA run

GFA mode does not write BED, BEDgraph, or TSV summary files. It writes one annotated graph:

```text
<input>.telo.annotated.gfa
```

Example:

```sh
teloscope asm.gfa -o results/
ls results/*.telo.annotated.gfa
```

## Pattern and motif problems

### `-c` or `-p` contains digits

You may see:

```text
Error: Canonical pattern '...' contains numerical characters.
Error: Pattern '...' contains numerical characters.
```

Teloscope expects nucleotide strings, not copy counts or labels.

Bad:

```sh
teloscope asm.fa -c TTAGGG2
teloscope asm.fa -p TTAGGG,TVR1
```

Good:

```sh
teloscope asm.fa -c TTAGGG
teloscope asm.fa -p TTAGGG,TCAGGG,TGAGGG
```

### Empty `-c` or `-p` fell back to defaults

You may see:

```text
Warning: Empty canonical pattern provided, using default vertebrate TTAGGG.
Warning: Empty pattern list provided, using default: TTAGGG, CCCTAA
Warning: No valid patterns supplied via -p. Using canonical: ...
```

This usually means a shell quoting or variable-expansion mistake produced an empty argument.

Common cause:

```sh
teloscope asm.fa -p "$PATTERNS"
```

where `PATTERNS` is empty. Check the variable before running the command.

### The `p` and `q` labels look wrong for your organism

`-c` sets the reference motif used for canonical counting and strand labeling. If `-c` is wrong, labels and classifications can look wrong even when repeat matching worked.

Examples:

```sh
teloscope asm.fa -c TTAGGG
teloscope asm.fa -c CCCTAAA
```

Use the canonical motif for the organism you are analyzing. Do not assume vertebrate `TTAGGG` for plant-like assemblies.

### Pattern count exploded and the run became slow

You may see:

```text
Warning: 1180 patterns is unusually high and may be slow on large genomes.
  Consider fewer IUPAC wildcards or a lower -x value.
```

Pattern count grows fast when you combine:

- many entries in `-p`
- IUPAC ambiguity codes
- larger `-x`

Examples:

| Input | `-x` | Concrete patterns after expansion |
| --- | --- | --- |
| `TTAGGG` | `1` | `38` |
| `TTAGGG` | `2` | `308` |
| `NNNGGG` | `0` | `127` |
| `NNNGGG` | `1` | `1180` |
| `NNNGGG` | `2` | `3367` |
| `NNNNNN` | `2` | `4096` |

If runtime is high:

- reduce ambiguity in `-p`
- lower `-x`
- keep `-c` fixed and search a smaller explicit motif set

### Edit distance is rejected

You may see:

```text
Error: Edit distance (-x/--edit-distance) must be in the range [0,2].
Error: Invalid edit distance '...'. Must be a number [0,2].
```

Allowed values are `0`, `1`, and `2`.

Rule of thumb:

- `-x 0`: exact repeats only
- `-x 1`: sensible default for most runs
- `-x 2`: broader search, slower, more permissive

## Window and binning problems

### Window size or step size is invalid

You may see:

```text
Error: Window size (-w or --window) must be > 0.
Error: Step size (-s or --step) must be > 0.
Error: Invalid window size '...'. Must be a number.
Error: Invalid step size '...'. Must be a number.
```

Use positive integers:

```sh
teloscope asm.fa -w 1000 -s 1000
```

### A threshold flag is not a number or is out of range

You may see one of these:

```text
Error: Terminal limit (-t or --terminal-limit) must be > 0.
Error: Invalid terminal limit '...'. Must be a number.
Error: Max match distance (-k/--max-match-distance) must be > 0.
Error: Invalid max match distance '...'. Must be a number.
Error: Min block length (-l or --min-block-length) must be > 0.
Error: Invalid min block length '...'. Must be a number.
Error: Max block distance (-d or --max-block-distance) must be > 0.
Error: Invalid max block distance '...'. Must be a number.
Error: Min block density (-y/--min-block-density) must be in the range [0,1].
Error: Invalid min block density '...'. Must be a number [0,1].
```

Use:

- positive integers for `-t`, `-k`, `-d`, and `-l`
- a numeric fraction from `0` to `1` for `-y`

Examples:

```sh
teloscope asm.fa -t 50000 -k 50 -d 200 -l 500 -y 0.5
```

### Step size is larger than window size

You may see:

```text
Error: Step size (2000) cannot be larger than window size (1000).
```

Keep `-s <= -w`.

Use:

- `-w 1000 -s 1000` for non-overlapping windows
- `-w 1000 -s 250` for overlapping windows

### Your BEDgraph tracks look denser than expected

If `-s` is smaller than `-w`, windows overlap. That is valid, but it changes how the tracks look in IGV or UCSC.

If you want one value per bin, keep:

```sh
teloscope asm.fa -r -w 1000 -s 1000
```

## Threshold and classification surprises

### No telomeres were called

If the summary says `none` or the BED file is empty, the usual causes are:

- `-t` is too small, so terminal blocks fall outside the allowed end window
- `-l` is too large, so shorter blocks are discarded
- `-y` is too strict, so lower-density blocks are discarded
- `-x 0` is too strict for a noisy assembly
- the organism-specific canonical motif in `-c` is wrong

Start with a permissive debugging run:

```sh
teloscope asm.fa -t 100000 -l 200 -y 0.3 -x 1 --verbose
```

Then tighten one parameter at a time.

### Telomeres are split into too many small blocks

The usual cause is that merging is too strict.

Relevant flags:

- `-k`: maximum distance between repeat matches
- `-d`: maximum distance between repeat groups during block extension

If `-k` or `-d` is too small, nearby matches stop merging. Try increasing them gradually:

```sh
teloscope asm.fa -k 100 -d 400
```

### Too much sequence merged into one block

The opposite problem can happen when `-k` or `-d` is too large. If unrelated repeat runs are being fused, lower one or both values:

```sh
teloscope asm.fa -k 25 -d 100
```

### Blocks disappear after changing `-l` or `-y`

These two flags are the most common reason a previously visible block vanishes.

- `-l` filters by absolute block length
- `-y` filters by repeat-covered fraction

If a run suddenly loses borderline blocks, compare:

```sh
teloscope asm.fa -l 500 -y 0.5
teloscope asm.fa -l 200 -y 0.3
```

### Classification says `incomplete`, `misassembly`, or `discordant` when you expected `t2t`

Check these first:

- Is `-c` correct for the organism?
- Is `-t` large enough to catch the terminal block?
- Are you looking at scaffold-terminal blocks only?
- Did you expect contig-terminal blocks to appear in the BED file?

`*_report.tsv` classifies sequences from scaffold-terminal blocks. Use `-n/--manual-curation` if you also want contig-terminal blocks retained in the BED output for review.

### You expected more BED rows with `-n`

`-n/--manual-curation` changes BED retention, not the core scaffold classification rules. It keeps contig-terminal telomeres in the BED output, but it does not change how `t2t`, `incomplete`, `misassembly`, `discordant`, or `none` are assigned.

## FASTA output flag interactions

### `-u` and genome-wide output flags seem to disagree

`-u/--ultra-fast` scans sequence ends only. That is the default.

The genome-wide output flags:

- `-r`
- `-g`
- `-e`
- `-m`
- `-i`

automatically disable ultra-fast mode because they require full-sequence analysis.

If runtime suddenly increases after enabling one of those flags, that is expected.

### ITS file is missing

`*_interstitial_telomeres.bed` is only written with `-i`.

Use:

```sh
teloscope asm.fa -i -o results/
```

### Match BED files are missing

`*_canonical_matches.bed` and `*_noncanonical_matches.bed` are only written with `-m`.

Use:

```sh
teloscope asm.fa -m -o results/
```

### Repeat-density and ratio tracks are missing

These files are only written with `-r`:

- `*_window_repeat_density.bedgraph`
- `*_window_canonical_ratio.bedgraph`
- `*_window_strand_ratio.bedgraph`

Use:

```sh
teloscope asm.fa -r -o results/
```

## GFA-specific problems

### The GFA run produced no telomere nodes

Common causes:

- the segment sequence is `*`, so there is nothing to scan
- the detected block fails `-l` or `-y`
- paths are present and the segment end is not path-terminal
- the motif in `-c` does not match the assembly

Check the output graph:

```sh
rg "telomere_" results/asm.gfa.telo.annotated.gfa
```

### Only some segment ends were annotated

When GFA paths are present, Teloscope annotates path-terminal segment ends. It does not annotate every segment end just because that segment contains terminal-like repeat content.

If you expected every segment to be scanned, inspect whether the graph has `P` lines and whether the segment end is terminal in any path.

### You expected BED or report files from a GFA run

GFA mode writes the annotated graph only. This is normal:

```sh
teloscope asm.gfa -o results/
```

Expected output:

```text
results/asm.gfa.telo.annotated.gfa
```

## Report generation problems

### `--plot-report` failed because Python packages are missing

The report generator needs:

- `matplotlib`
- `numpy`
- `pandas`

Install them in the Python environment that Teloscope will use.

### Teloscope could not locate `teloscope_report.py`

You may see:

```text
Warning: Could not locate teloscope_report.py.
Warning: Report generation failed.
```

This usually means the script is not available relative to the binary location or working tree.

Manual workaround:

```sh
python3 scripts/teloscope_report.py results/
```

### The standalone report script says files are missing

You may see:

```text
Error: Missing teloscope output files in '...'
```

Point the script at the Teloscope output directory, not the repo root:

```sh
python3 scripts/teloscope_report.py results/
```

At minimum the directory must contain `*_terminal_telomeres.bed`.

## Threading and performance

### The run is slower than expected on a shared machine

By default `-j` uses all available threads. That is fine on a dedicated node, but it can oversubscribe laptops or shared servers.

Cap the thread count explicitly:

```sh
teloscope asm.fa -j 4
```

### The run is fast without optional flags, then slow with them

That usually means you crossed from terminal-only scanning into full-sequence scanning. This happens when you enable any of:

- `-r`
- `-g`
- `-e`
- `-m`
- `-i`

If you only need terminal calls and the summary table, keep the run minimal:

```sh
teloscope asm.fa -o results/
```

## Quick diagnosis checklist

If a run looks wrong, check these in order:

1. Confirm the input mode: FASTA or GFA.
2. Confirm the organism-specific canonical motif with `-c`.
3. Confirm which optional outputs you actually requested.
4. Check whether `-r`, `-g`, `-e`, `-m`, or `-i` disabled ultra-fast mode.
5. Check whether `-t`, `-l`, `-y`, `-k`, or `-d` became too strict.
6. Rerun once with `--cmd --verbose`.
