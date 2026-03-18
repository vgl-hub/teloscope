[← Back to README](../README.md)

How it works
============

Briefly, **Teloscope** reads an assembly and decomposes its parts. Input patterns (with optional IUPAC codes and edit distance) are expanded into all concrete variants along with their reverse complements, and a multi-pattern search tree is built for fast matching. Each contig is scanned: in ultra-fast mode (default), only the regions near each contig end are scanned; otherwise, the full sequence is analyzed.

Nearby repeat matches are merged into groups, and nearby groups are further merged into telomere blocks. Short blocks and blocks with low repeat density are discarded. Each surviving block is labeled by strand orientation (`p`, `q`, or `b`), and the set of blocks for a chromosome determines its classification (T2T, incomplete, misassembled, discordant, or none). Terminal telomere BED annotations are always written. Optional window metrics, match coordinates, and ITS blocks are produced when the corresponding flags are enabled.
