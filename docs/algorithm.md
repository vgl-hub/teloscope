[Back to README](../README.md)

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
5. Merge nearby matches into repeat groups, then merge nearby groups into telomere blocks.
6. Filter blocks by minimum length and minimum repeat density.
7. Label each surviving block as `p`, `q`, or `b` from strand composition.
8. Mark blocks as scaffold-terminal or contig-terminal from `-t/--terminal-limit`.
9. Classify the sequence as `t2t`, `incomplete`, `misassembly`, `discordant`, or `none`.
10. Write BED, TSV, and optional BEDgraph outputs.

## FASTA scanning modes

By default Teloscope runs in ultra-fast mode. It scans only the first and last `-t` base pairs of each sequence. This is usually enough for terminal telomere annotation and keeps whole-genome runs fast.

If any genome-wide output flag is enabled (`-r`, `-g`, `-e`, `-m`, or `-i`), ultra-fast mode is disabled automatically. In that case Teloscope scans the full sequence and can report ITS blocks, genome-wide windows, and individual matches.

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

`--fastq-subset` reads FASTQ records in bounded batches, scans each read as a whole sequence with the same pattern expansion and block filters used by FASTA mode, and writes unchanged passing FASTQ records to stdout. Read order is preserved. By default records stream to stdout so the output can be piped straight into a mapper; pass `-o` to save them to `<output>/<input>_telomeric.fastq` instead. Diagnostics and final counts are written to stderr. FASTQ subset mode defaults to a 60 bp minimum block length, while assembly annotation keeps the 500 bp default.

## BAM subset mode

`--bam-subset` reads BGZF-compressed BAM directly through `zlib`, without HTSlib or command-line converters. It preserves the BAM header, scans each record's stored `SEQ`, and writes passing records unchanged. Primary, secondary, supplementary, mapped, and unmapped records are evaluated independently. Records without `SEQ` are dropped and counted separately.

Records are processed in bounded byte and record batches. Worker threads score sequences while the main thread writes passing records in input order. The output is valid BGZF with an EOF marker; no index is copied or generated. Missing input EOF markers produce a warning, while malformed BGZF or BAM data is rejected.

FASTQ and BAM use the same read-scoring wrapper and 60 bp default. BAM I/O is isolated from scoring so a future SAM parser or optional CRAM backend can reuse the same filter.

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
