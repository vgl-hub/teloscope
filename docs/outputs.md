[Back to README](../README.md)

# Outputs

Teloscope writes different outputs in FASTA mode and GFA mode.

## File naming

FASTA outputs keep the input file name as a prefix:

- `asm.fa_terminal_telomeres.bed`
- `asm.fa_gaps.bed`
- `asm.fa_report.tsv`

GFA mode writes one graph file:

- `asm.gfa.telo.annotated.gfa`

## FASTA mode outputs

Always written:

| File | Purpose |
| --- | --- |
| `*_terminal_telomeres.bed` | terminal telomere blocks |
| `*_gaps.bed` | gap intervals from runs of `N` |
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

## `*_terminal_telomeres.bed`

Columns:

1. `chr`
2. `start`
3. `end`
4. `length`
5. `label`
6. `fwdCount`
7. `revCount`
8. `canonCount`
9. `nonCanonCount`
10. `chrSize`
11. `type`

`label` is `p`, `q`, or `b`. `type` is `scaffold` for scaffold-terminal blocks and `contig` for contig-terminal blocks that are only emitted with `-n/--manual-curation`.

## `*_gaps.bed`

Columns:

1. `chr`
2. `start`
3. `end`

Each row marks one contiguous run of `N`.

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

## GFA mode output

GFA mode writes the original graph plus synthetic telomere segments and links:

- segment lines start with `S telomere_...`
- link lines start with `L telomere_...`

Each synthetic telomere segment includes:

- `LN:i:6`
- `RC:i:6000`
- `TL:i:<detected_block_length_bp>`

Each telomere link points from the synthetic node to the assembly segment with `0M` overlap. When paths are present, only path-terminal segment ends are annotated and the link orientation follows the path context. When no paths are present, segments are scanned independently.

No BED, BEDgraph, or TSV files are written in GFA mode.
