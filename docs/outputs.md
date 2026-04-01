[← Back to README](../README.md)

Outputs
============

### GFA mode

When the input is a GFA file, Teloscope scans each segment for terminal telomeric repeats and writes an annotated GFA:

* `assembly.telo.annotated.gfa` The original graph with synthetic telomere nodes appended. Each telomere node is a placeholder segment (`S telomere_<seg>_<start|end> * LN:i:6 RC:i:6000 TL:i:<len>`) linked to the parent assembly segment (`L telomere_... + <seg> <orient> 0M`). The `TL` tag records the actual telomere block length in bp. Designed for visualization in BandageNG.

No BED or BEDgraph files are produced in GFA mode.

### FASTA mode

Teloscope always produces the following files:

* `terminal_telomeres.bed` Telomere block annotations for the assembly.

    Columns: `chr`, `start`, `end`, `length`, `label`, `fwdCount`, `revCount`, `canonCount`, `nonCanonCount`, `chrSize`, `type`.

    The `label` column assigns each block to an arm: `p` (mostly forward-strand matches), `q` (mostly reverse-strand), or `b` (balanced). The `type` column is either `scaffold` (block is near a scaffold end) or `contig` (block is near a contig end but not a scaffold end). By default only scaffold-terminal telomeres are reported; use `-n`/`--manual-curation` to also include contig-terminal blocks.

* `gaps.bed` Assembly gap coordinates (runs of Ns). Columns: `chr`, `start`, `end`.

* `report.tsv` The Path Summary and Assembly Summary tables (also printed to stdout), written as a TSV file for downstream parsing.

Additional optional outputs (disabled in ultra-fast mode):

* `window_repeat_density.bedgraph` Fraction of each window covered by any repeat match, 0 to 1 (`-r`).
* `window_canonical_ratio.bedgraph` Canonical share of repeat density per window, 0 to 1; -1 when no matches (`-r`).
* `window_strand_ratio.bedgraph` Forward-strand share of repeat density per window, 0 to 1; -1 when no matches (`-r`).
* `window_gc.bedgraph` GC content per window (`-g`).
* `window_entropy.bedgraph` Shannon entropy per window (`-e`).
* `canonical_matches.bed` Coordinates of every canonical repeat match in the assembly (`-m`).
* `noncanonical_matches.bed` Coordinates of non-canonical repeat matches in terminal regions (`-m`).
* `interstitial_telomeres.bed` Telomere-like blocks found away from contig ends, i.e., interstitial telomeres or ITSs (`-i`).

All output filenames are prefixed with the input filename (e.g., `asm.fa_terminal_telomeres.bed`).

### Summary report

Teloscope prints two tables to stdout and writes them to `report.tsv`:

**Path Summary Report.** One row per sequence. Columns in ultra-fast mode: `pos`, `header`, `telomeres`, `labels`, `gaps`, `type`, `granular`. In full-scan mode, three columns are added: `its`, `canonical`, `windows`.

**Assembly Summary Report.** Aggregated counts including:
* Total paths, gaps, and telomeres
* Telomere length statistics (mean, median, min, max)
* Chromosome counts by telomere number (two, one, zero)
* Chromosome classification counts (see [Classification](classification.md))
