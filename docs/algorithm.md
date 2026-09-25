[Back to README](index.md)

# Teloscope algorithm

Teloscope has three input modes:

- FASTA mode scans sequence ends, groups telomeric matches into blocks, and classifies each path or scaffold.
- GFA mode scans graph segments and writes an annotated graph with synthetic telomere nodes.
- Reads mode (FASTQ or BAM input) measures and subsets reads in one pass.

## FASTA mode

1. Read the input assembly.
2. Expand the requested repeat patterns, including IUPAC codes and allowed edit-distance variants.
3. Add reverse complements of every pattern.
4. Build a multi-pattern search structure and scan each sequence.
5. Split each scanned region into contigs — runs of called bases with no `N` — and find the canonical matches and cover runs of each strand inside every contig.
6. At the scaffold's two ends (first contig start, last contig end; with `-n`, every contig end), walk the exact canonical repeats inward. The first one inside the start zone (the smaller of `--terminal-tolerance` and `-t`) anchors the telomere. A block runs from the anchor inward to where its strand's canonical coverage still averages `-y`, bridging at most `-d` of non-telomeric sequence, and stops at the outer edge of a real array of the other strand. A block must be strand-pure: at least the `--label-threshold` share of its exact repeats on its own strand. The next exact repeat within `-d` of the block's end starts another block of the same strand, or ends the chain if it begins a real array of the other strand. `teloLen` is the sum of the block lengths; it is smaller than `end - start` only when the telomere is fragmented. A telomere never crosses a run of `N`, and `-t` never bounds it.
7. Outside the telomeres, `-k` chains matches into seeds, and same-strand seeds within `-d` join into one span. The span is trimmed to where all-match coverage averages `-y`. It is kept with at least four exact canonical repeats and coverage of at least `-y`. Its junction class — `fusion`, `tail_to_tail`, `fragmentation` or `single`, defined in [outputs](outputs.md) — is judged against the nearest row within `-d` on the same contig.
8. Label every block with its strand composition and with the arm it sits on.
9. Classify the sequence as `t2t`, `incomplete`, or `none`, and separately record any orientation anomalies.
10. Write BED, TSV, and optional BEDgraph outputs.

## FASTA scanning modes

A scaffold's own ends are the start of its first contig and the end of its last contig; leading and trailing runs of `N` don't count.

By default Teloscope runs in fast mode. It reads only the first contig's head window and the last contig's tail window — the first `-t` bases of each. The window grows inward in steps of `-t` while a repeat lies within reach of its inner edge. This builds the p and q arms and is usually enough for terminal telomere annotation, while keeping whole-genome runs fast.

Any of `-r`, `-g`, `-e`, `-m`, or `-i` forces the full scan: every contig is read whole, and every array outside the arms is an ordinary interstitial row.

`-n/--manual-curation` keeps fast mode but reads both end windows of every contig and builds both chains on each. A telomere at an internal contig end becomes a contig row in the terminal BED instead of an interstitial row.

The interstitial file is always written. In fast mode it holds what the end windows found: by default the first contig's head and the last contig's tail, with `-n` both end windows of every contig. Nothing beyond them is scanned until you add `-i`.

## GFA mode

1. Read the graph header, segments, links, and paths.
2. Decide which segment ends are valid scan targets.
3. Scan the available segment sequence for terminal telomeric repeats.
4. Create one synthetic telomere segment for each detected terminal block.
5. Connect each synthetic node back to the matching assembly segment end with an `L` link at `0M` overlap, the direct adjacency a cap represents.
6. Write the result as `<input>.telo.annotated.gfa`.

When paths are present, Teloscope annotates only path-terminal segment ends. This keeps the graph output aligned with assembly path ends rather than every raw segment end. When no paths are present, each segment is treated independently.

Synthetic telomere nodes are placeholders. They carry tags that preserve the detected telomere length while keeping the graph easier to display in BandageNG.

## Reads mode

Each read is scanned twice with the FASTA-mode patterns: the assembly block rule (`-l`) gives its BED rows, and a whole-read scan with a fixed 42 bp floor decides whether it is kept. A read with a BED row is always kept. BAM reverse-strand records are measured in read orientation; secondary, supplementary and hard-clipped records are not measured but can be kept.

## Pattern handling

`-c/--canonical` and `-p/--patterns` have different roles:

- `-c` sets the reference repeat used for canonical versus non-canonical counting and for `p` or `q` labeling.
- `-p` sets the actual search set.

If `-p` is omitted, Teloscope derives the search set from `-c`. Reverse complements are always added automatically.

## Windowed outputs

Windowed outputs are optional and apply to FASTA mode only:

- repeat density
- canonical ratio
- strand ratio
- GC content
- entropy

Window size comes from `-w`. Step size comes from `-s`. When `-s` equals `-w`, the output is a standard non-overlapping BEDgraph track.
