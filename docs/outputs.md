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
#params canonical=<forward>/<reverse> patterns=<count> window=<bp> step=<bp> terminal_limit=<bp> max_match_dist=<bp> max_block_dist=<bp> min_block_len=<bp> min_block_density=<fraction> min_block_counts=<n> min_its_length=<bp> terminal_tolerance=<bp> edit_distance=<n> ultra_fast=<true-or-false> manual_curation=<true-or-false>
#columns	<tab-separated column names>
```

`patterns` is the number of search patterns after expansion and deduplication. `commit` is the short commit checked out when the binary was built, or `unknown` when Git metadata was unavailable. It does not indicate whether that commit was built from a clean worktree. The normal report on standard output is unchanged.

BED files contain BED records only. BEDGraph files begin with their standard browser `track` declaration followed by four-column data. Keeping run metadata in the companion report avoids comment headers that strict coordinate converters reject. Telomere block files still have explicitly documented custom fields, so schema-aware converters must be told they are BED3 plus fifteen custom fields.

## Telomere block BED files

Both block files use zero-based, half-open coordinates and share one 18-column layout. Columns 1 to 11 hold the v0.1.5 fields in their v0.1.5 positions, with one change worth knowing: column 5 now means the strand in both files, where v0.1.5 put the scan direction there for terminal blocks. It can also read `b` on a terminal row, which v0.1.5 never produced, so a consumer that switches on `p` or `q` at column 5 needs a case for it. Columns 12 to 18 are new in v0.1.6. Column 18, `status`, reads `canonical`, `discordant` or `balanced`, so an orientation anomaly has coordinates rather than only a scaffold name. Only the first three fields are standard BED; the rest are Teloscope-specific, so readers that only understand predefined BED fields should consume the first three columns. UCSC bigBed conversion requires a matching AutoSql definition for the custom fields. See the [UCSC BED specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [bigBed custom-field documentation](https://genome.ucsc.edu/goldenPath/help/bigBed.html).

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `chrom` | FASTA record ID |
| 2 | `chromStart` | block start |
| 3 | `chromEnd` | block end |
| 4 | `blockLen` | `chromEnd - chromStart` |
| 5 | `label` | strand composition: `p`, `q`, or `b` |
| 6 | `forwardCount` | forward matches, canonical and variant |
| 7 | `reverseCount` | reverse matches, canonical and variant |
| 8 | `canonicalCount` | exact canonical matches, both strands |
| 9 | `nonCanonicalCount` | variant matches, both strands |
| 10 | `chromSize` | full FASTA record length |
| 11 | `blockType` | `scaffold` in the terminal file, `contig` in the interstitial file |
| 12 | `arm` | position: `p` for the start side, `q` for the end side |
| 13 | `gapStatus` | `gapped` if the block spans an assembly gap, `contiguous` otherwise |
| 14 | `fwdCan` | exact forward canonical matches |
| 15 | `revCan` | exact reverse canonical matches |
| 16 | `fwdNonCan` | forward-oriented variant matches |
| 17 | `revNonCan` | reverse-oriented variant matches |
| 18 | `status` | `canonical`, or `discordant` / `balanced` when the block's strand does not match its arm |

Columns 6 to 9 are the row and column totals of columns 14 to 17. They are kept so that consumers written against v0.1.5 column positions do not have to change, and because the joint cells cannot be recovered from the totals. `scripts/teloscope_report.py` reads this layout and the interim 10-column one.

The four count cells are mutually exclusive. All useful marginal totals are derived from them:

```text
forward matches       = fwdCan + fwdNonCan
reverse matches       = revCan + revNonCan
canonical matches     = fwdCan + revCan
non-canonical matches = fwdNonCan + revNonCan
all matches           = fwdCan + revCan + fwdNonCan + revNonCan
```

Canonical matched bases are `(fwdCan + revCan) * canonical motif length`, and canonical block density divides that value by `chromEnd - chromStart`. With a fixed-length pattern set, total matched-base coverage can be derived the same way. Counts alone do not determine total matched bases when one run mixes non-canonical motifs of different lengths.

Column 5 records strand composition in both files: `p` means more than 66.6% of all matches in the block are forward-oriented, `q` means less than 33.3%, and `b` is the interval between those thresholds. Column 12 records position, which end of the sequence the block sits nearer. In v0.1.5 column 5 held the strand for interstitial blocks and the scan direction for terminal ones; it now means the same thing in both files, and the position it used to carry for terminal blocks moved to column 12. For an ordinary telomere the two agree. They differ on an inverted or balanced array, which the report's granular column marks with `*` and `~` respectively, and column 18 names outright.

“Forward” is a sequence-family convention, not a reference `+` strand annotation. Teloscope orders the canonical motif and its reverse complement lexicographically and calls the smaller string forward. With the default motif pair, forward is `CCCTAA`, normally seen at a chromosome start, and reverse is `TTAGGG`, normally seen at a chromosome end. Each concrete seed after IUPAC expansion is assigned to the closer canonical orientation (ties go to forward), and its edit-distance variants inherit that orientation.

`blockType=scaffold` marks a row in the terminal file and `blockType=contig` a row in the interstitial file. A terminal block has to start within `--terminal-tolerance` called bases of a record end, so every terminal row is scaffold-terminal. Whether a block spans an assembly gap is a separate question and lives in column 13.

## `*_gaps.bed`

Columns:

1. `chr`
2. `start`
3. `end`

Each row marks one contiguous run of `N`, `n`, `X`, or `x`. Blocks are built once per FASTA record, so a block may bridge a run of `N` no longer than `-d/--max-block-distance`. Column 13 of both block files, `gapStatus`, marks the rows where that happened. A longer run still ends the block, and so does a run of called non-telomeric sequence longer than the same limit, which is what keeps two arrays at one end from merging into a single telomere.

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
- flagged and clean scaffolds, and the discordant arms, balanced arms and extra terminal blocks behind them

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
