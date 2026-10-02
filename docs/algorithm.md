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
4. Build a telomere at each of the scaffold's two ends, or at every contig end with `-n`:
   - `-k` chains one strand's matches, exact and variant: a match joins while it starts within `-k` of the chain's end.
   - The first exact canonical repeat inside the start zone anchors a piece. The piece runs inward along its chain: a base covered by an exact repeat of its strand scores 1, any other base costs `y / (1 - y)` for `-y` of `y`, and the piece ends where the running score peaks, the furthest peak on a tie. A dense array is neither stretched into a sparse tail nor dropped because of one, and what is left of the chain goes to step 5.
   - A piece must hold `--min-block-counts` exact repeats and be strand-pure: more than `--label-threshold` of its exact repeats on its own strand. It stops at the outer edge of a real array of the other strand.
   - The next exact repeat of the same strand, in another chain and within `-d` of the piece's end, anchors another piece. It joins only if it outweighs the gap between the two: its length must be at least the gap times `y / (1 - y)`, so a small array far behind a telomere stays an interstitial row. A real array of the other strand ends the telomere; an array is real when pieces joined this way from its edge sum to `-l`.
   - `teloLen` sums the pieces and must reach `-l`. Pieces that fall short are dropped, and the search resumes behind the first of them. A telomere never crosses `N`, and `-t` never bounds it.
5. Outside the telomeres, `-k` chains the matches of both strands into seeds, and three opposite-strand matches in a row start a new seed. Each seed is cut the same way on the coverage of all its matches, starting again behind every peak, and a cut part is kept with at least four exact canonical repeats. Its [junction class](outputs.md#telomere-block-bed) comes from the rows within `-d` on the same contig.
6. Label every block by strand and by the end it belongs to, then [classify](classification.md) the sequence.
7. Write the BED files, the report, and the optional window tracks.

## FASTA scanning modes

A scaffold's ends are the start of its first contig and the end of its last contig.

- **Fast**, the default: reads the first `-t` bases of the first contig and the last `-t` bases of the last contig, growing inward by `-t` while a repeat lies within reach of the inner edge. This builds both arms and keeps whole-genome runs fast. The interstitial file holds what these windows found.
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
