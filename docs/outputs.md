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
| `*_terminal_telomeres.bed` | terminal telomere rows |
| `*_interstitial_telomeres.bed` | interstitial telomere-like rows |
| `*_gaps.bed` | gap intervals from runs of `N`/`n`/`X`/`x` |
| `*_report.tsv` | per-sequence and assembly summaries |

`*_interstitial_telomeres.bed` is always written. In fast mode it holds only what the head and tail scan windows found; `-i` or `-n` forces the full scan so it holds every array.

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
| `*_plot_report.pdf` | `--plot-report` | PDF summary report |

## Run provenance

The always-written `*_report.tsv` starts with three comment lines:

```text
#teloscope version=0.1.6 commit=<short-commit-or-unknown>
#params canonical=<forward>/<reverse> patterns=<count> window=<bp> step=<bp> terminal_limit=<bp> max_match_dist=<bp> max_block_dist=<bp> min_block_len=<bp> min_block_density=<fraction> min_block_counts=<n> min_canonical_count=<n> terminal_tolerance=<bp> link_distance=<bp> label_threshold=<fraction> edit_distance=<n> ultra_fast=<true-or-false> manual_curation=<true-or-false>
#columns	<tab-separated column names>
```

`patterns` is the number of search patterns after expansion and deduplication. `commit` is the short commit checked out when the binary was built, or `unknown` when Git metadata was unavailable. It does not indicate whether that commit was built from a clean worktree. The normal report on standard output is unchanged.

BED files contain BED records only. BEDGraph files begin with their standard browser `track` declaration followed by four-column data. Keeping run metadata in the companion report avoids comment headers that strict coordinate converters reject. Telomere block files still have explicitly documented custom fields, so schema-aware converters must be told they are BED3 plus nine custom fields.

## Telomere block BED files

Both block files use zero-based, half-open coordinates and share one 12-column layout. Only the first three fields are standard BED; the rest are Teloscope-specific, so readers that only understand predefined BED fields should consume the first three columns. UCSC bigBed conversion requires a matching AutoSql definition for the custom fields. See the [UCSC BED specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [bigBed custom-field documentation](https://genome.ucsc.edu/goldenPath/help/bigBed.html).

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `chr` | FASTA record ID |
| 2 | `start` | block start |
| 3 | `end` | block end |
| 4 | `teloLen` | sum of the piece lengths |
| 5 | `teloLabel` | strand: `p` (forward, the `CCCTAA` family) or `q` (reverse, `TTAGGG`); `b` only on interstitial rows |
| 6 | `closestEnd` | for a telomere, the end it belongs to (`p` or `q`); for an interstitial row, the nearer record end, `p` on an exact midpoint |
| 7 | `fwdCan` | exact forward canonical matches |
| 8 | `revCan` | exact reverse canonical matches |
| 9 | `fwdNonCan` | forward-oriented variant matches |
| 10 | `revNonCan` | reverse-oriented variant matches |
| 11 | `chrSize` | full FASTA record length |
| 12 | `teloType` | `scaffold` or `contig` in the terminal file; junction class in the interstitial file: `fusion`, `tail_to_tail`, `fragmentation`, or `single` |

`teloLen` equals `end - start` for every interstitial row, and for a terminal telomere built from a single piece. Only a fragmented terminal telomere — one built from more than one piece — has `teloLen` smaller than `end - start` (`fragmented_p`/`fragmented_q` in the report). A block never contains a gap: blocks are built per contig, and a contig is by definition a run of called bases with no `N` inside.

Canonical coverage is the number of bases covered by exact canonical repeats, with overlaps counted once. A terminal piece must average at least `-y` canonical coverage over its own length to qualify. An interstitial row must average at least `-y` all-repeat coverage — canonical and variant matches together — over its own length. Counts alone do not give you the coverage, since matches can overlap and non-canonical motifs can differ in length.

For an ordinary telomere `teloLabel` and `closestEnd` agree: a `p` row carries the forward motif and sits near the start, a `q` row the reverse motif near the end. They differ on a discordant arm, which the report's `granular` column marks with `*` and the `anomaly` column names as `discordant_p`/`discordant_q`.

“Forward” is a sequence-family convention, not a reference `+` strand annotation. Teloscope orders the canonical motif and its reverse complement lexicographically and calls the smaller string forward. With the default motif pair, forward is `CCCTAA`, normally seen at a chromosome start, and reverse is `TTAGGG`, normally seen at a chromosome end. Each concrete seed after IUPAC expansion is assigned to the closer canonical orientation (ties go to forward), and its edit-distance variants inherit that orientation.

`teloType=scaffold` marks an arm row in the terminal file. `teloType=contig` marks a contig-terminal row in the terminal file; those only appear with `-n`. In the interstitial file `teloType` is always a junction class, judged against the nearest row within `--link-distance` on the same contig: `fusion` for a reverse array then a forward one, `tail_to_tail` for forward then reverse, `fragmentation` for two arrays of the same strand, and `single` when nothing is that close.

## `*_gaps.bed`

Columns:

1. `chr`
2. `start`
3. `end`

Each row marks one contiguous run of `N`, `n`, `X`, or `x`. Blocks are built per contig — the called sequence between gaps — so a block never spans a gap. `-d/--max-block-distance` bounds how much non-telomeric sequence a single piece may bridge; `--link-distance` decides whether two pieces close enough join into one telomere instead of staying separate.

To attach the nearest gap to each terminal and interstitial block with BEDTools, first select their common BED3 prefix:

```sh
cut -f1-3 asm.fa_terminal_telomeres.bed > blocks.bed
cut -f1-3 asm.fa_interstitial_telomeres.bed >> blocks.bed
bedtools closest -a blocks.bed -b asm.fa_gaps.bed -d > blocks_with_nearest_gap.tsv
```

The output contains the three block fields, the three gap fields, and the block-to-gap distance.

## `*_report.tsv`

The report contains two tables.

Path Summary columns in ultra-fast mode:

- `pos`
- `header`
- `telomeres`
- `labels`
- `gaps`
- `type`
- `anomaly`
- `granular`

Full-scan mode adds three more columns; fast mode never reports `its`:

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
- flagged and clean scaffolds, and the discordant and fragmented arms behind them

The summary line `Fragmented arms` replaces the old `Extra terminal blocks` and `Balanced arms` lines. Contig-terminal rows, written to the terminal BED only with `-n`, never count toward `type`, `telomeres`, the anomalies, or the length statistics in either table.

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
