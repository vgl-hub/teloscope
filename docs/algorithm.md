[Back to README](index.md)

# Algorithm

Teloscope has three input modes:

- FASTA mode scans sequence ends, groups telomeric matches into blocks, and classifies each sequence.
- GFA mode scans graph segments and writes the graph back with telomere caps.
- Reads mode (FASTQ or BAM) measures and subsets reads in one pass.

## FASTA mode

1. Read the assembly.
2. Expand the patterns: IUPAC codes, `-x` substitutions, and reverse complements.
3. Split each sequence into contigs at gaps (runs of `N`, `n`, `X`, or `x`) and scan each contig with a multi-pattern trie, collecting each strand's canonical matches and the stretches they cover.
4. Build a telomere at each of the scaffold's two ends, or at every contig end with `-n`. One score runs through every step: a base covered by an exact repeat of the strand in hand earns `1 - y` and any other base costs `y`, for `-y` of `y`, so a stretch scores its exact-repeat fraction minus `y`, times its length: zero or more exactly when exact repeats cover at least `-y` of it. The score is kept in integers (millionths), so ties and thresholds are exact.
   - **Blocks.** `-k` groups one strand's matches into blocks, exact and variant: a match joins while it starts within `-k` of the block's end. A block needs `--min-block-counts` matches.
   - **Pieces.** A block that meets `-y` as a whole is one piece, variant repeats included, so a telomere can end on a variant repeat. A block below `-y` is cut into its maximal scoring segments (Ruzzo and Tompa), each a piece: a dense array is neither stretched into a sparse tail nor dropped because of one. A piece must be strand-pure, more than `--label-threshold` of its exact repeats on its own strand.
   - **Walk.** The first piece of either strand inside the start zone starts the telomere, which walks inward over the pieces of its strand: each piece adds its score and each base between two pieces costs `y`. The walk only reaches a piece whose block lies within `-d` of the last block, and it stops at a real array of the other strand (one whose own walk spans `-l`). The telomere ends where the summed score peaks, the furthest peak on a tie. So every part of the telomere that runs to its inner end meets `-y`, and nothing within `-d` beyond that end would.
   - **Length.** `-l` applies to the span, start to end. A walk that falls short of `-l` is dropped and the search resumes at the next piece in the zone, so a stub cannot hide a telomere. A telomere never crosses `N`, and `-t` never bounds it.
   - `teloLen` sums the pieces; a telomere of more than one piece is `fragmented`. What lies outside the telomere goes to step 5.
5. Outside the telomeres, each taken from its start to its end, `-k` groups the matches of both strands into seeds, and three opposite-strand matches in a row start a new seed. Each seed is cut into its maximal scoring segments on the coverage of all its matches, and a segment is kept with at least four exact canonical repeats. Its [junction class](outputs.md#telomere-block-bed) comes from the rows within `-d` on the same contig.
6. Label every block by strand and by the end it belongs to, then [classify](classification.md) the sequence.
7. Write the BED files, the report, and the optional window tracks.

## FASTA scanning modes

A scaffold's ends are the start of its first contig and the end of its last contig.

- **Fast**, the default: reads the first `-t` bases of the first contig and the last `-t` bases of the last contig, growing inward while the telomere builder still looks that far. This builds both arms, identical to the full scan, and keeps whole-genome runs fast. The interstitial file holds what these windows found.
- **Full**: any of `-r`, `-g`, `-e`, `-m`, or `-i` reads every contig whole, and every other array that passes step 5 becomes an interstitial row.
- **Contig ends** (`-n`): keeps the fast scan but reads both end windows of every contig. A telomere at an internal contig end becomes a `contig` row in the terminal BED instead of an interstitial row.

## GFA mode

1. Read the header, segments, links, and paths.
2. Pick the segment ends to scan: path-terminal ends when there are paths, every segment end otherwise.
3. Scan each end for a terminal telomere.
4. Add one cap segment per telomere and join it to its segment end with an `L` link at `0M` overlap.
5. Write `<input>.telo.annotated.gfa`.

Caps are placeholders: a tag keeps the telomere length, and the segment stays small in BandageNG.

## Reads mode

Each read's ends are measured with the assembly rules (`-l`, `-t 2000`, `--terminal-tolerance 300`), which give its BED rows. A read with no row is scanned once more for any telomeric block of at least 42 bp, which keeps it in the subset. BAM reverse-strand records are measured in read orientation. Secondary, supplementary, and hard-clipped records are not measured but can be kept.
