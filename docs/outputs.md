[Back to README](index.md)

# Outputs

Every file name starts with the input file name, as in `asm.fa_report.tsv`. Only the kept-reads BAM drops the extension: `reads_telomeric.bam`.

## FASTA mode

| File | Flag | Content |
| --- | --- | --- |
| `*_terminal_telomeres.bed` |  | telomeres at sequence ends |
| `*_interstitial_telomeres.bed` |  | interstitial telomeres (ITS) |
| `*_gaps.bed` |  | runs of `N`, `n`, `X`, or `x` |
| `*_report.tsv` |  | per-sequence and assembly summaries, as printed on stdout |
| `*_window_repeat_density.bedgraph` | `-r` | repeat density per window |
| `*_window_canonical_ratio.bedgraph` | `-r` | canonical share of the repeats per window |
| `*_window_strand_ratio.bedgraph` | `-r` | forward-strand share of the repeats per window |
| `*_window_gc.bedgraph` | `-g` | GC content per window |
| `*_window_entropy.bedgraph` | `-e` | Shannon entropy per window |
| `*_canonical_matches.bed` | `-m` | canonical matches |
| `*_noncanonical_matches.bed` | `-m` | non-canonical matches near sequence ends |
| `*_terminal_telomeres.fa` | `-a` | terminal telomere sequences |
| `*_plot_report_terminal.pdf`, `*_plot_report_its.pdf` | `--plot-report` | [PDF reports](report.md) |

The fast scan fills the interstitial file only from the ends it reads; `-i` scans everything ([scanning modes](algorithm.md#fasta-scanning-modes)). BED files have no header and BEDgraph files start with a `track` line; run metadata lives in the report.

## Telomere block BED

Both block files share one 12-column layout, with zero-based, half-open coordinates. Columns 1 to 3 are standard BED and the rest are Teloscope fields, so a bigBed conversion needs an AutoSql file for BED3+9.

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `chr` | FASTA record ID |
| 2 | `start` | block start |
| 3 | `end` | block end |
| 4 | `teloLen` | sum of the piece lengths |
| 5 | `teloLabel` | strand: `p` forward (`CCCTAA`), `q` reverse (`TTAGGG`), `b` mixed (interstitial rows only) |
| 6 | `closestEnd` | end the telomere belongs to, `p` or `q`; for an interstitial row, the nearer end, `p` at the midpoint |
| 7 | `fwdCan` | exact forward canonical matches |
| 8 | `revCan` | exact reverse canonical matches |
| 9 | `fwdNonCan` | forward variant matches |
| 10 | `revNonCan` | reverse variant matches |
| 11 | `chrSize` | record length |
| 12 | `teloType` | terminal file: `scaffold` for an arm, `contig` for a contig end (`-n`); interstitial file: junction class |

`teloLen` equals `end - start` except on a fragmented telomere, one built from more than one piece. No block crosses a gap.

The junction class compares an interstitial row with the nearest row within `-d` on the same contig: `fusion` is a reverse array then a forward one, `tail_to_tail` forward then reverse, `fragmentation` two arrays of the same strand, and `single` means nothing is that close.

Coverage counts each base under a match once, however many matches overlap it, so counts alone do not give it. A terminal piece needs `-y` canonical coverage over its length; an interstitial row needs `-y` coverage from all its matches.

Forward and reverse are a sequence convention, not a `+` strand: of the canonical motif and its reverse complement, the lexicographically smaller one is forward. By default that is `CCCTAA`, found at chromosome starts, and reverse is `TTAGGG`, found at ends. Each expanded seed takes the closer orientation (ties go to forward), and its variants inherit it.

## `*_report.tsv`

Three comment lines come first:

```text
#teloscope version=0.1.6 commit=<short-commit-or-unknown>
#params canonical=<forward>/<reverse> patterns=<count> window=<bp> step=<bp> terminal_limit=<bp> max_match_dist=<bp> max_block_dist=<bp> min_block_len=<bp> min_block_density=<fraction> min_block_counts=<n> min_canonical_count=<n> terminal_tolerance=<bp> label_threshold=<fraction> edit_distance=<n> ultra_fast=<true-or-false> manual_curation=<true-or-false>
#columns	<column names>
```

`commit` is the commit the binary was built from, and `patterns` counts the search patterns after expansion. Then comes the per-sequence table:

| Column | Meaning |
| --- | --- |
| `pos` | position in the input, from 1; skips records removed by a filter |
| `header` | sequence ID |
| `telomeres` | arms found: `0`, `1`, or `2` |
| `labels` | ends with an arm, in order: `p`, `q`, `pq`, or `none` |
| `gaps` | number of gaps |
| `type` | `t2t`, `incomplete`, or `none` ([Classification](classification.md)) |
| `anomaly` | `.`, or a comma-separated list of `discordant_p`, `discordant_q`, `fragmented_p`, `fragmented_q` |
| `its` | interstitial rows, full scan only |

The assembly summary follows:

| Section | Lines |
| --- | --- |
| Assembly Summary | paths, gaps, scaffold and contig N50, telomeres; ITS blocks in a full scan; input and selected paths under a filter |
| Telomere Statistics | mean, median, min, and max arm length |
| Chromosome Telomere Counts | sequences with two, one, or zero arms |
| Chromosome Telomere/Gap Completeness | `t2t`, `incomplete`, and `none`, each split by gaps |
| Scaffold Anomalies | flagged and clean sequences, and the discordant and fragmented arms behind them |

Contig rows (`-n`) never count in the report. Under a filter, every number describes the selected paths.

## `*_gaps.bed`

One BED3 row per run of `N`, `n`, `X`, or `x`. To find the gap nearest each block with BEDTools, which needs sorted input:

```sh
cut -f1-3 asm.fa_terminal_telomeres.bed asm.fa_interstitial_telomeres.bed | sort -k1,1 -k2,2n > blocks.bed
sort -k1,1 -k2,2n asm.fa_gaps.bed > gaps.bed
bedtools closest -a blocks.bed -b gaps.bed -d > nearest_gap.tsv
```

Each output row holds the block, the gap, and their distance.

## `*_terminal_telomeres.fa`

One record per row of the terminal BED, named `>{chr}:{start}-{end}` with the row's own coordinates. The sequence is the block span in uppercase on one line; a fragmented telomere is written whole, spacers included.

## GFA mode

`<input>.telo.annotated.gfa` is the input graph plus one telomere segment per detected end, joined to its segment by an `L` link at `0M` overlap. `J` records stay reserved for real gaps. `<input>.telo.annotated.colors.csv` paints every cap green (`#008000`) for BandageNG.

Each cap carries `LN:i:6`, `RC:i:6000`, and `TL:i:<telomere length>`. Its link runs from the cap (`+`) to the segment end (`+` at a start, `-` at an end) and carries `RC:i:0`.

With `P` paths, only path-terminal segment ends are scanned, oriented by the path. Without paths, every segment is scanned on its own. Caps belong to the graph, so a shared segment end shows its cap in every path. GFA mode writes no BED or report files.

## Reads mode

- `*_telomeric.fastq` or `*_telomeric.bam`: the kept reads, unchanged and in input order.
- `*_terminal_telomeres.bed`: one row per measured read telomere, in the [block layout](#telomere-block-bed). `chr` is the read name, `chrSize` the read length, and `teloType` is `read`. A read can have a `p` row, a `q` row, both, or none.
- `*_report.tsv`: the comment lines, then `label<TAB>value` lines, also printed on stdout.

| Line | Meaning |
| --- | --- |
| `Reads measured` | FASTQ: every read; BAM: primary records with a sequence and no hard clip at either end |
| `Reads kept` | records written to the subset; can exceed `Reads measured` on an aligned BAM |
| `Read telomeres` | BED rows written |
| `Complete` | strand matches the end (C-rich at a read start, G-rich at a read end), and the read continues at least `-d` past the telomere |
| `Reaching read end` | strand matches, but the read ends within `-d` of the telomere |
| `Discordant` | strand does not match the end |
| `Mean length` to `Max length` | over `Complete` rows: mean, median, 25th, 75th, and 90th percentiles (linear interpolation, two decimals), min, and max; left out when there are none |

The estimate needs no alignment, so a read broken inside a telomere looks complete and pulls the lengths down.
