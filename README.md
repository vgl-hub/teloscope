Teloscope
============
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2025.10.14.682431.svg)](https://doi.org/10.1101/2025.10.14.682431)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Introduction
------------
Teloscope is a comprehensive telomere annotation tool. It rapidly matches, counts, and reports telomeric repeats from genome assemblies (.fa) or (.fa.gz). Teloscope reports all these metrics in BED/BEDgraph files and produces a summary report. To install `teloscope`, use: 
```sh
git clone https://github.com/vgl-hub/teloscope.git --recursive;
cd teloscope; 
make -j
```

Usage
------------
    teloscope input.[fa][.gz] [options]
    teloscope -f input.[fa][.gz] -o [output/dir] -j [threads] -c [canonical] -p [patterns] -w [window size] -s [step size] -d [max-block-dist] -l [min-block-len]

The input file can be passed as the first positional argument or with `-f`.

**Note:** Teloscope automatically explores the input repeats and their reverse complements. If none are provided, it will scan for the canonical CCCTAA/TTAGGG repeats. 

Examples
------------
* Example: Minimal single case to run Teloscope.

        teloscope -f "${file}"

* Example: Multiple input case to run Teloscope. 

        teloscope -f "${file}" -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG

* Example: Set window and step sizes. Calculating window metrics.

        teloscope -f "${file}" -o "${out_path}" -j 16 -c TTAGGG -p NNNGGG -w 2000 -s 1000 -g -e -r --verbose

* Example: Allowing all outputs and using all flags.

        teloscope -f "${file}" -o "${out_path}" -j 16 -c TTAGGG -p TBAGGG,TTRGGG,YTAGGG  -w 2000 -s 1000 -d 200 -l 1000 -r -g -e -m -i -t 50000 --verbose --cmd
  
**Note:** Teloscope accepts nucleotides in IUPAC format and generates all possible pattern combinations. 

Parameters 
------------

To check out all options and flags, please use:
`teloscope -h` 

```
Required Parameters:
        '-f'    --input-sequence        Input fasta file (or pass as first positional argument).
        '-o'    --output        Set output route. [Default: Input path]
        '-c'    --canonical     Set canonical pattern. [Default: TTAGGG]
        '-p'    --patterns      Set patterns to explore, separate them by commas [Default: TTAGGG]
        '-j'    --threads       Set maximum number of threads. [Default: max. available]
        '-t'    --terminal-limit        Set terminal limit for exploring telomere variant regions (TVRs). [Default: 50000]
        '-k'    --max-match-distance    Set maximum distance for merging matches. [Default: 50]
        '-d'    --max-block-distance    Set maximum block distance for extension. [Default: 200]
        '-l'    --min-block-length      Set minimum block length. [Default: 500]
        '-y'    --min-block-density     Set minimum block density. [Default: 0.5]

Optional Parameters:
        '-w'    --window        Set sliding window size. [Default: 1000]
        '-s'    --step  Set sliding window step. [Default: 1000 (non-overlapping)]
        '-r'    --out-win-repeats       Output canonical/noncanonical repeats and density by window. [Default: false]
        '-g'    --out-gc        Output GC content for each window. [Default: false]
        '-e'    --out-entropy   Output Shannon entropy for each window. [Default: false]
        '-m'    --out-matches   Output all canonical and terminal non-canonical matches. [Default: false]
        '-i'    --out-its       Output assembly interstitial telomere (ITSs) regions.[Default: false]
        '-u'    --ultra-fast    Ultra-fast mode. Only scans terminal telomeres at contig ends. [Default: true]
        '-n'    --manual-curation    Retain all terminal telomeres (contig + scaffold) in BED output. [Default: scaffold only]
        '-v'    --version       Print current software version.
        '-h'    --help  Print current software options.
        --verbose       Verbose output.
        --cmd   Print command line.
```

Outputs
------------
Teloscope outputs telomere annotations in BED format:

* `terminal_telomeres.bed` Annotation of the full telomere in the assembly. Columns: chr, start, end, length, label, fwdCounts, revCounts, canCounts, nonCanCounts, chrSize, type (scaffold/contig). By default only scaffold-terminal telomeres are reported; use `--manual-curation` to include contig-terminal blocks.

Additional optional outputs (Disabled with --ultra-fast):

* `interstitial_telomeres.bed` Blocks of adjacent canonical repeat matches. Outside of the ends, it represents interstitial telomeres (ITSs).
* `window_fwd.bedgraph` Forward repeat counts per window (-r). IGV-compatible single-value BEDgraph.
* `window_rev.bedgraph` Reverse repeat counts per window (-r).
* `window_canonical.bedgraph` Canonical repeat counts per window (-r).
* `window_noncanonical.bedgraph` NonCanonical repeat counts per window (-r).
* `window_gc.bedgraph` GC content per window (-g).
* `window_entropy.bedgraph` Shannon entropy per window (-e).
* `canonical_matches.bed` Coordinates of canonical repeats throughout the assembly.
* `noncanonical_matches.bed` Coordinates of non-canonical repeats in terminal regions of contigs.

**Note:** Window step defaults to window size (non-overlapping windows). This produces valid single-value BEDgraph files compatible with IGV and UCSC genome browsers. Use `-s` to set a smaller step for overlapping windows if needed.

How it works
------------
Briefly, **Teloscope** reads an assembly and decomposes its parts. It uses prefix trees and sliding windows to efficiently perform multiple string matching and counting of telomeric repeats. It analyzes the informational properties of the sliding windows to find telomeric blocks. These blocks are collected, post-processed, and filtered according to their positional and conformational properties. 

How to cite
------------

**Teloscope** is part of the gfastar tool suite. 
If you use **Teloscope** in your research, please cite:

**The complete genome of a songbird**  
Giulio Formenti, Nivesh Jain, **Jack A. Medico**, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Inês Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceição, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O’Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnéa Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis
bioRxiv 2025.10.14.682431; doi: https://doi.org/10.1101/2025.10.14.682431




