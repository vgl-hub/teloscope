[Back to README](index.md)

# Teloscope algorithm

Teloscope has two input modes:

- FASTA mode scans sequence ends, groups telomeric matches into blocks, and classifies each path or scaffold.
- GFA mode scans graph segments and writes an annotated graph with synthetic telomere nodes.
- FASTQ subset mode streams reads and writes only records with Teloscope-valid telomeric blocks.
- BAM subset mode streams alignments and writes only records whose stored `SEQ` has a valid block.

## FASTA mode

1. Read the input assembly.
2. Expand the requested repeat patterns, including IUPAC codes and allowed edit-distance variants.
3. Add reverse complements of every pattern.
4. Build a multi-pattern search structure and scan each sequence.
5. Split each scanned region into contigs — runs of called bases with no `N` — and find the canonical matches and cover runs of each strand inside every contig.
6. At the scaffold's own two ends — the first contig's start and the last contig's end — walk the canonical matches inward (with `-n`, every contig end is walked this way). The first one whose piece qualifies, inside the start zone (the smaller of `--terminal-tolerance` and `-t`), anchors the telomere. An anchor is an exact repeat; a variant repeat at the very tip is not part of the telomere. A piece is the anchor trimmed inward to where that strand's canonical coverage still averages `-y`, bridging at most `-d` of called non-telomeric sequence, and stopped at the outer edge of a real array of the other strand (real: an array that would itself reach `-l`). A piece must also be strand-pure: at least the `--label-threshold` share (default `0.667`) of its exact repeats must be on its own strand, so a stretch that alternates strands every repeat is never part of a telomere and becomes an interstitial row labelled `b` instead. Chain onward: a canonical match starting within `--link-distance` of the last piece either extends the chain with another piece of the same strand, or, when it is a real array of the other strand, ends the chain and becomes an interstitial row instead. `teloLen` is the sum of the piece lengths, so it is less than the block's own length exactly when the telomere is fragmented. A telomere never crosses a run of `N`, and `-t` never bounds how far it can reach — it only sets the start zone and the fast-mode window.
7. Score everything outside the telomeres as interstitial rows: `-k` chains matches into seeds, and consecutive same-strand seeds within `-d` join into one span. Trim the span's edges to where all-match coverage still averages `-y`. Keep it when it holds at least four exact canonical repeats and its all-match coverage is at least `-y` times its own length; there is no minimum length. Label its junction against the nearest row within `--link-distance` on the same contig: a reverse array then a forward one is `fusion`, forward then reverse is `tail_to_tail`, two arrays of the same strand are `fragmentation`, and no near neighbour is `single`.
8. Label every block with its strand composition and with the arm it sits on.
9. Classify the sequence as `t2t`, `incomplete`, or `none`, and separately record any orientation anomalies.
10. Write BED, TSV, and optional BEDgraph outputs.

## FASTA scanning modes

A scaffold's own ends are the start of its first contig and the end of its last contig; leading and trailing runs of `N` don't count.

By default Teloscope runs in fast mode. It reads only the first contig's head window and the last contig's tail window — the first `-t` bases of each. The window grows inward in steps of `-t` while a repeat lies within reach of its inner edge, so `-t` never bounds a result. This builds the p and q arms and is usually enough for terminal telomere annotation, while keeping whole-genome runs fast.

Any of `-r`, `-g`, `-e`, `-m`, or `-i` forces the full scan: every contig is read whole, and every array outside the arms is an ordinary interstitial row.

`-n/--manual-curation` also forces the full scan, and additionally builds both chains on every contig. A telomere at an internal contig end becomes a `contig` row in the terminal BED instead of an interstitial row.

The interstitial file is always written. In fast mode it holds only what the head and tail windows found; nothing beyond them is scanned until you add `-i` or `-n`.

## GFA mode

1. Read the graph header, segments, links, and paths.
2. Decide which segment ends are valid scan targets.
3. Scan the available segment sequence for terminal telomeric repeats.
4. Create one synthetic telomere segment for each detected terminal block.
5. Connect each synthetic node back to the matching assembly segment end with an `L` link at `0M` overlap, the direct adjacency a cap represents.
6. Write the result as `<input>.telo.annotated.gfa`.

When paths are present, Teloscope annotates only path-terminal segment ends. This keeps the graph output aligned with assembly path ends rather than every raw segment end. When no paths are present, each segment is treated independently.

Synthetic telomere nodes are placeholders. They carry tags that preserve the detected telomere length while keeping the graph easier to display in BandageNG.

## FASTQ subset mode

`--fastq-subset` reads FASTQ records in bounded batches, scans each read as a whole sequence with the same pattern expansion and block filters used by FASTA mode, and writes unchanged passing FASTQ records to stdout. Read order is preserved. By default records stream to stdout so the output can be piped straight into a mapper; pass `-o` to save them to `<output>/<input>_telomeric.fastq` instead. Diagnostics and final counts are written to stderr. FASTQ subset mode defaults to a 42 bp minimum block length, while assembly annotation keeps the 300 bp default.

## BAM subset mode

`--bam-subset` reads BGZF-compressed BAM directly through `zlib`, without HTSlib or command-line converters. It preserves the BAM header, scans each record's stored `SEQ`, and writes passing records unchanged. Primary, secondary, supplementary, mapped, and unmapped records are evaluated independently. Records without `SEQ` are dropped and counted separately.

Records are processed in bounded byte and record batches. Worker threads score sequences while the main thread writes passing records in input order. The output is valid BGZF with an EOF marker; no index is copied or generated. Missing input EOF markers produce a warning, while malformed BGZF or BAM data is rejected.

FASTQ and BAM use the same read-scoring wrapper and 42 bp default. BAM I/O is isolated from scoring so a future SAM parser or optional CRAM backend can reuse the same filter.

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
