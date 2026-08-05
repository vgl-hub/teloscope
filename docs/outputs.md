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

## Run provenance

The always-written `*_report.tsv` starts with three comment lines:

```text
#teloscope version=0.1.6 commit=<short-commit-or-unknown>
#params canonical=<forward>/<reverse> patterns=<count> window=<bp> step=<bp> terminal_limit=<bp> max_match_dist=<bp> max_block_dist=<bp> min_block_len=<bp> min_block_density=<fraction> edit_distance=<n> ultra_fast=<true-or-false> manual_curation=<true-or-false>
#columns	<tab-separated column names>
```

`patterns` is the number of search patterns after expansion and deduplication. `commit` is the short commit checked out when the binary was built, or `unknown` when Git metadata was unavailable. It does not indicate whether that commit was built from a clean worktree. The normal report on standard output is unchanged.

BED files contain BED records only. BEDGraph files begin with their standard browser `track` declaration followed by four-column data. Keeping run metadata in the companion report avoids comment headers that strict coordinate converters reject. Telomere block files still have explicitly documented custom fields, so schema-aware converters must be told their BED4+ field count.

## Telomere block BED files

Both block files use zero-based, half-open coordinates. Their first four fields form a standard BED4 prefix; column 4 is the BED `name` field and contains Teloscope's `p`, `q`, or `b` label. Later fields are Teloscope-specific, rather than BED `score` or `strand` fields with invented semantics. Readers that only understand predefined BED fields should consume the first four columns. Schema-aware readers should declare these files as BED4+5 or BED4+6 as appropriate. UCSC bigBed conversion additionally requires a matching AutoSql definition for custom fields. See the [UCSC BED specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [bigBed custom-field documentation](https://genome.ucsc.edu/goldenPath/help/bigBed.html).

Schema migration: the v0.1.6 layout moves the label from column 5 to the standard BED `name` field in column 4 and replaces the former length and marginal-count fields with the four joint counts below. Positional consumers must update their indices; `scripts/teloscope_report.py` accepts both layouts during migration.

Common BED4+5 fields:

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `chrom` | FASTA record ID |
| 2 | `chromStart` | block start |
| 3 | `chromEnd` | block end; block length is `chromEnd - chromStart` |
| 4 | `name` | `p`, `q`, or `b` block label |
| 5 | `fwdCan` | exact forward canonical matches |
| 6 | `revCan` | exact reverse canonical matches |
| 7 | `fwdNonCan` | forward-oriented variant matches |
| 8 | `revNonCan` | reverse-oriented variant matches |
| 9 | `chromSize` | full FASTA record length |

`*_terminal_telomeres.bed` is BED4+6: it adds column 10, `blockType`, whose value is `scaffold` or `contig`. The interstitial file needs no type column because every row in that file is interstitial.

The four count cells are mutually exclusive. All useful marginal totals are derived from them:

```text
forward matches       = fwdCan + fwdNonCan
reverse matches       = revCan + revNonCan
canonical matches     = fwdCan + revCan
non-canonical matches = fwdNonCan + revNonCan
all matches           = fwdCan + revCan + fwdNonCan + revNonCan
```

Canonical matched bases are `(fwdCan + revCan) * canonical motif length`, and canonical block density divides that value by `chromEnd - chromStart`. With a fixed-length pattern set, total matched-base coverage can be derived the same way. Counts alone do not determine total matched bases when one run mixes non-canonical motifs of different lengths.

For an interstitial block, `p` means more than 66.6% of all matches are forward-oriented, `q` means less than 33.3%, and `b` is the interval between those thresholds. Terminal blocks retain `p` or `q` from the end-specific scan; their labels are not recalculated from the counts.

“Forward” is a sequence-family convention, not a reference `+` strand annotation. Teloscope orders the canonical motif and its reverse complement lexicographically and calls the smaller string forward. With the default motif pair, forward is `CCCTAA`, normally seen at a chromosome start, and reverse is `TTAGGG`, normally seen at a chromosome end. Each concrete seed after IUPAC expansion is assigned to the closer canonical orientation (ties go to forward), and its edit-distance variants inherit that orientation.

`blockType=scaffold` marks a block within `terminal_limit` of a FASTA record end. `blockType=contig` marks a block found by a terminal scan of an internal ungapped-segment end; these rows are emitted only with `-n/--manual-curation`.

## `*_gaps.bed`

Columns:

1. `chr`
2. `start`
3. `end`

Each row marks one contiguous run of `N`, `n`, `X`, or `x`. Teloscope splits a FASTA record into ungapped segments at every such run and builds blocks within one segment at a time. A terminal or interstitial block therefore cannot span a row in `*_gaps.bed`.

A gap acts as a contig boundary. Teloscope never fuses repeats on opposite sides of an unknown run into one block, even when the match-merging distance is longer than the run.

To attach the nearest gap to each terminal and interstitial block with BEDTools, first select their common BED4 prefix:

```sh
cut -f1-4 asm.fa_terminal_telomeres.bed > blocks.bed
cut -f1-4 asm.fa_interstitial_telomeres.bed >> blocks.bed
bedtools closest -a blocks.bed -b asm.fa_gaps.bed -d > blocks_with_nearest_gap.tsv
```

The output contains the four block fields, the three gap fields, and the block-to-gap distance.

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
