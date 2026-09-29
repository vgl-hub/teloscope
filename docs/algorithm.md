[Back to README](index.md)

# Algorithm

Teloscope has three input modes:

- FASTA mode scans sequence ends, groups telomeric matches into blocks, and classifies each sequence.
- GFA mode scans graph segments and writes the graph back with telomere caps.
- Reads mode (FASTQ or BAM) measures and subsets reads in one pass.

## FASTA mode

1. Read the assembly.
2. Expand the patterns: IUPAC codes, `-x` substitutions, and reverse complements.
3. Build a multi-pattern trie and scan each sequence.
4. Split each scanned region into contigs, runs of called bases with no `N`, and find the canonical matches and cover runs of each strand inside every contig.
5. Build a telomere at each of the scaffold's two ends, or at every contig end with `-n`:
   - The first exact canonical repeat inside the start zone anchors it.
   - A block runs inward while its strand's canonical coverage still averages `-y`, bridges at most `-d` of non-telomeric sequence, and stops at the outer edge of a real array of the other strand.
   - A block must be strand-pure: at least `--label-threshold` of its exact repeats on its own strand.
   - The next exact repeat within `-d` of the block's end starts another block of the same strand, or ends the chain if it begins a real array of the other strand.
   - `teloLen` sums the blocks. A telomere never crosses `N`, and `-t` never bounds it.
6. Outside the telomeres, `-k` chains matches into seeds, and same-strand seeds within `-d` join into one span. The span is trimmed to where all-match coverage averages `-y`, and kept with at least four exact canonical repeats. Its [junction class](outputs.md#telomere-block-bed) comes from the nearest row within `-d` on the same contig.
7. Label every block by strand and by the end it belongs to, then [classify](classification.md) the sequence.
8. Write the BED files, the report, and the optional window tracks.

## FASTA scanning modes

A scaffold's ends are the start of its first contig and the end of its last contig.

- **Fast**, the default: reads the first `-t` bases of the first contig and the last `-t` bases of the last contig, growing inward by `-t` while a repeat lies within reach of the inner edge. This builds both arms and keeps whole-genome runs fast. The interstitial file holds what these windows found.
- **Full**: any of `-r`, `-g`, `-e`, `-m`, or `-i` reads every contig whole, and every array outside the arms becomes an interstitial row.
- **Contig ends** (`-n`): keeps the fast scan but reads both end windows of every contig. A telomere at an internal contig end becomes a `contig` row in the terminal BED instead of an interstitial row.

## GFA mode

1. Read the header, segments, links, and paths.
2. Pick the segment ends to scan: path-terminal ends when there are paths, every segment end otherwise.
3. Scan each end for a terminal telomere.
4. Add one cap segment per telomere and join it to its segment end with an `L` link at `0M` overlap.
5. Write `<input>.telo.annotated.gfa`.

Caps are placeholders: a tag keeps the telomere length, and the segment stays small in BandageNG.

## Reads mode

Each read is scanned twice with the assembly patterns. The block rule (`-l`) gives its BED rows, and a whole-read scan with a fixed 42 bp floor decides whether it is kept; a read with a BED row is always kept. BAM reverse-strand records are measured in read orientation. Secondary, supplementary, and hard-clipped records are not measured but can be kept.
