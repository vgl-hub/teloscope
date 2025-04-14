Teloscope
============

Introduction
------------
Teloscope is a universal telomere annotation tool. It comprehensively runs matching, counting, and reporting telomeric repeats from genome assemblies (.fa) or (.fa.gz). Teloscope reports all these metrics in BED/BEDgraph files and produces a summary report. To install `teloscope`, use: 
```sh
git clone https://github.com/vgl-hub/teloscope.git --recursive;
cd teloscope; 
make -j
```

Usage
------------
    teloscope -f input.[fa][.gz] -o [output/dir] -j [threads] -c [canonical] -p [patterns] -w [window size] -s [step size] -d [max-block-dist] -l [min-block-len]

**Note:** Teloscope automatically explores the input repeats and their reverse complements. If none are provided, it will scan for the canonical CCCTAA/TTAGGG repeats. 

Examples
------------
* Example: Minimal single case to run Teloscope.

        teloscope -f "${file}" -o "${out_path}" -u

* Example: Multiple input case to run Teloscope. Set window and step sizes.

        teloscope -f "${file}" -o "${out_path}" -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG -w 2000 -s 1000

* Example: Calculating window metrics.

        teloscope -f "${file}" -o "${out_path}" -j 16 -c TTAGGG -p NNNGGG -w 1000 -s 500 -g -e -r --verbose

* Example: Allowing all outputs and using all flags.

        teloscope -f "${file}" -o "${out_path}" -j 16 -c TTAGGG -p TBAGGG,TTRGGG,YTAGGG  -w 2000 -s 1000 -d 200 -l 1000 -r -g -e -m -i -t 50000 --verbose --cmd
  
**Note:** Teloscope accepts nucleotides in IUPAC format and generates all possible pattern combinations. 

Parameters 
------------

To check out all options and flags, please use:
`teloscope -h` 

```
Required Parameters:
        '-f'    --input-sequence        Initiate tool with fasta file.
        '-o'    --output        Set output route.
        '-c'    --canonical     Set canonical pattern. [Default: TTAGGG]
        '-p'    --patterns      Set patterns to explore, separate them by commas [Default: TTAGGG]
        '-w'    --window        Set sliding window size. [Default: 1000]
        '-s'    --step  Set sliding window step. [Default: 500]
        '-j'    --threads       Set maximum number of threads. [Default: max. available]
        '-l'    --min-block-length      Set minimum block length for merging. [Default: 500]
        '-d'    --max-block-distance    Set maximum block distance for merging. [Default: 200]
        '-t'    --terminal-limit        Set terminal limit for exploring telomere variant regions (TVRs). [Default: 50000]

Optional Parameters:
        '-r'    --out-win-repeats       Output canonical/noncanonical repeats and density by window. [Default: false]
        '-g'    --out-gc        Output GC content for each window. [Default: false]
        '-e'    --out-entropy   Output Shannon entropy for each window. [Default: false]
        '-m'    --out-matches   Output all canonical and terminal non-canonical matches. [Default: false]
        '-i'    --out-its       Output assembly interstitial telomere (ITSs) regions.[Default: false]
        '-u'    --ultra-fast    Ultra-fast mode. Only scans terminal telomeres at contig ends. [Default: false]
        '-v'    --version       Print current software version.
        '-h'    --help  Print current software options.
        --verbose       Verbose output.
        --cmd   Print command line.
```

Outputs
------------
Teloscope outputs telomere annotations in BED format:

* `terminal_telomeres.bed` Annotation of the full telomere in the assembly. This is made of canonical and non-canonical repeats. 

Additional optional outputs (Disabled with --ultra-fast):

* `interstitial_telomeres.bed` Blocks of adjacent canonical repeat matches. Outside of the ends, it represents interstitial telomeres (ITSs).
* `window_metrics.bedgraph` Tabulated file with calculated window metrics such as forward repeats, reverse repeats, canonical repeats, non-canonical repeats (-r), GC% (-g) and Shannon Entropy (-e) by window. 
* `canonical_matches.bed` Coordinates of canonical repeats throughout the assembly. 
* `noncanonical_matches.bed` Coordinates of non-canonical repeats in terminal regions of contigs. 

How it works
------------
Briefly, **Teloscope** reads an assembly and decomposes its parts. It uses prefix trees and sliding windows to efficiently perform multiple string matching and counting of telomeric repeats. It analyzes the informational properties of the sliding windows to find telomeric blocks. These blocks are collected, post-processed, and filtered according to their positional and conformational properties. 

How to cite
------------

**Teloscope** is part of the gfastar tool suite. If you use **Teloscope** in your research, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs.
Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1093/bioinformatics/btac460
