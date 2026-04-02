[Back to README](../README.md)

# Teloscope algorithm

Teloscope has two input modes:

- FASTA mode scans sequence ends, groups telomeric matches into blocks, and classifies each path or scaffold.
- GFA mode scans graph segments and writes an annotated graph with synthetic telomere nodes.

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
5. Link each synthetic node back to the matching assembly segment end.
6. Write the result as `<input>.telo.annotated.gfa`.

When paths are present, Teloscope annotates only path-terminal segment ends. This keeps the graph output aligned with assembly path ends rather than every raw segment end. When no paths are present, each segment is treated independently.

Synthetic telomere nodes are placeholders. They carry tags that preserve the detected telomere length while keeping the graph easy to display in BandageNG.

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
