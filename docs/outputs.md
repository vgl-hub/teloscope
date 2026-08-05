[Back to README](index.md)

# Outputs

Teloscope writes different outputs in assembly and read subset modes.

## File naming

FASTA outputs keep the input file name as a prefix:

- `asm.fa_terminal_telomeres.bed`
- `asm.fa_gaps.bed`
- `asm.fa_report.tsv`

GFA mode writes one graph file:

- `asm.gfa.telo.annotated.gfa`

BAM subset mode with `-o results/` writes:

- `results/reads_telomeric.bam`

## FASTA mode outputs

Always written:

| File | Purpose |
| --- | --- |
| `*_terminal_telomeres.bed` | terminal telomere blocks |
| `*_gaps.bed` | gap intervals from runs of `N`/`n`/`X`/`x` |
| `*_report.tsv` | per-sequence and assembly summaries |

Optional:

| File | Flag | Purpose |
| --- | --- | --- |
| `*_window_repeat_density.bedgraph` | `-r` | repeat density per window |
| `*_window_canonical_ratio.bedgraph` | `-r` | canonical-share track per window |
| `*_window_strand_ratio.bedgraph` | `-r` | forward-strand share per window |
| `*_window_gc.bedgraph` | `-g` | GC content per window |
| `*_window_entropy.bedgraph` | `-e` | Shannon entropy per window |
| `*_canonical_matches.bed` | `-m` | canonical repeat matches |
| `*_noncanonical_matches.bed` | `-m` | terminal non-canonical repeat matches |
| `*_interstitial_telomeres.bed` | `-i` | interstitial telomere-like blocks |
| `*_plot_report.pdf` | `--plot-report` | PDF summary report |

## Provenance headers

Every BED, BEDGraph, and `*_report.tsv` file written in FASTA mode starts with three comment lines:

```text
#teloscope version=0.1.6 commit=<short-commit-or-unknown>
#params canonical=<forward>/<reverse> patterns=<count> window=<bp> step=<bp> terminal_limit=<bp> max_match_dist=<bp> max_block_dist=<bp> min_block_len=<bp> min_block_density=<fraction> edit_distance=<n> ultra_fast=<true-or-false> manual_curation=<true-or-false>
#columns	<tab-separated column names>
```

`patterns` is the number of search patterns after expansion and deduplication. `commit` is the short commit checked out when the binary was built, or `unknown` when Git metadata was unavailable. It does not indicate whether that commit was built from a clean worktree. The headers are written to files only; the normal report on standard output is unchanged. BEDGraph `track` lines follow these headers.

## Telomere block BED files

`*_terminal_telomeres.bed` and `*_interstitial_telomeres.bed` have the same 15-column layout. Coordinates are zero-based, half-open, so `length = end - start`.

| Column | Name | Terminal BED | Interstitial BED |
| ---: | --- | --- | --- |
| 1 | `chr` | FASTA record ID | FASTA record ID |
| 2 | `start` | block start | block start |
| 3 | `end` | block end | block end |
| 4 | `length` | block span in bp | block span in bp |
| 5 | `label` | `p` or `q`, set by the terminal scan direction | `p`, `q`, or `b`, set from all forward/reverse matches |
| 6 | `fwdCount` | forward-oriented matches | forward-oriented matches |
| 7 | `revCount` | reverse-oriented matches | reverse-oriented matches |
| 8 | `canonCount` | exact canonical matches | exact canonical matches |
| 9 | `nonCanonCount` | variant matches | variant matches |
| 10 | `chrSize` | full FASTA record length | full FASTA record length |
| 11 | `blockType` | `scaffold` or `contig` | `interstitial` |
| 12 | `canFwd` | exact forward canonical matches | exact forward canonical matches |
| 13 | `canRev` | exact reverse canonical matches | exact reverse canonical matches |
| 14 | `canCov` | canonical matched bp divided by `length` | canonical matched bp divided by `length` |
| 15 | `repCov` | all matched bp divided by `length` | all matched bp divided by `length` |

`canCov` and `repCov` are written to four decimal places. Matched bp is the sum of match lengths, not the size of their interval union. Matches may overlap, so `repCov` can exceed `1.0`; `canCov` can also exceed `1.0` for a self-overlapping custom canonical motif.

The count columns obey these identities:

```text
canFwd + canRev = canonCount
fwdCount + revCount = canonCount + nonCanonCount
```

For an interstitial block, `p` means more than 66.6% of all matches are forward-oriented, `q` means less than 33.3%, and `b` is the interval between those thresholds. Terminal blocks retain `p` or `q` from the end-specific scan; their labels are not recalculated from the counts.

“Forward” is a sequence-family convention, not a reference `+` strand annotation. Teloscope orders the canonical motif and its reverse complement lexicographically and calls the smaller string forward. With the default motif pair, forward is `CCCTAA`, normally seen at a chromosome start, and reverse is `TTAGGG`, normally seen at a chromosome end. Each concrete seed after IUPAC expansion is assigned to the closer canonical orientation (ties go to forward), and its edit-distance variants inherit that orientation.

`blockType=scaffold` marks a block within `terminal_limit` of a FASTA record end. `blockType=contig` marks a block found by a terminal scan of an internal ungapped-segment end; these rows are emitted only with `-n/--manual-curation`. The shared layout means the two files can be combined without reshaping, while column 11 retains the block context.

## `*_gaps.bed`

Columns:

1. `chr`
2. `start`
3. `end`

Each row marks one contiguous run of `N`, `n`, `X`, or `x`. Teloscope splits a FASTA record into ungapped segments at every such run and builds blocks within one segment at a time. A terminal or interstitial block therefore cannot span a row in `*_gaps.bed`.

To attach the nearest gap to each terminal and interstitial block with BEDTools:

```sh
awk '!/^#/' asm.fa_terminal_telomeres.bed asm.fa_interstitial_telomeres.bed > blocks.data.bed
awk '!/^#/' asm.fa_gaps.bed > gaps.data.bed
bedtools closest -a blocks.data.bed -b gaps.data.bed -d > blocks_with_nearest_gap.tsv
```

The output contains the 15 block columns, the three gap columns, and the block-to-gap distance. This join is optional: column 11 already distinguishes scaffold-terminal, contig-terminal, and interstitial blocks.

## `*_report.tsv`

The report contains two tables.

Path Summary columns in ultra-fast mode:

- `pos`
- `header`
- `telomeres`
- `labels`
- `gaps`
- `type`
- `granular`

Full-scan mode adds:

- `its`
- `canonical`
- `windows`

Assembly Summary reports totals and counts for:

- paths
- gaps
- telomeres
- telomere length statistics
- chromosomes with two, one, or zero telomeres
- each scaffold class

When an assembly record filter is active, the summary also reports the number of input and selected paths. `Total paths` and every other statistic describe the selected paths only.

## GFA mode output

GFA mode writes the original graph plus synthetic telomere segments connected back to the carrier segment with `L telomere_...` links at `0M` overlap. It also writes `<input>.telo.annotated.colors.csv`, which paints every cap green (`#008000`) for BandageNG.

Each synthetic telomere segment includes:

- `LN:i:6`
- `RC:i:6000`
- `TL:i:<detected_block_length_bp>`

Each telomere link points from the synthetic node (`+`) to the assembly segment, with the segment-side orientation set so the cap sits on the correct physical end (`+` at a start, `-` at an end), and carries `RC:i:0`. `J` jump records stay reserved for real assembly gaps. With GFA1 `P` paths, only terminal segment ends reached from selected paths are scanned and orientation follows the path context. In a pathless GFA1 graph, selected segments are scanned independently. Caps are graph-level, so a shared annotated segment end is visible from selected and excluded paths. Filtered GFA2, GFA1 `C` containment and `W` walk records, and unknown GFA record types are rejected.

No BED, BEDgraph, or TSV files are written in GFA mode. Record filters do not delete supported original graph records; they limit which terminal segment ends Teloscope scans for new annotations.

## BAM subset output

BAM subset output contains the original BAM header and unchanged passing records in input order. Records without `SEQ` are omitted. The output has a standard BGZF EOF marker but no `.bai` or `.csi` index.
