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
    teloscope -f input.[fa][.gz] -o [output/dir] -j [threads] -c [canonical] -p [patterns] -w [window size] -s [step size] -d [max-block-dist] -l [min-block-len] -k

**Note:** Teloscope automatically explores the input repeats and their reverse complements. If none are provided, it will scan for the canonical CCCTAA/TTAGGG repeats. 

Examples
------------
* Example:

        teloscope -f "${file}" -o "${out_path}" -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG -w 2000 -s 1000 -k

* Example:

        teloscope -f "${file}" -o "${out_path}" -j 16 -c TTAGGG -p NNNGGG -w 1000 -s 500 -d 200 -l 1000 -k --verbose

* Example:

        teloscope -f "${file}" -o "${out_path}" -j 16 -c TTAGGG -p TBAGGG,TTRGGG,YTAGGG  -w 2000 -s 1000 -d 200 -l 1000 -k --verbose
  
**Note:** Teloscope accepts nucleotides in IUPAC format and generates all possible pattern combinations. 

Parameters 
------------

To check out all options and flags, please use:
`teloscope -h` 

```
Required Parameters:
        '-f'    --input-sequence        Initiate tool with fasta/fasta.gz file.
        '-o'    --output        Set output route.
        '-c'    --canonical     Set canonical pattern. [Default: TTAGGG]
        '-p'    --patterns      Set patterns to explore, separate them by commas [Default: TTAGGG]
        '-w'    --window        Set sliding window size. [Default: 1000]
        '-s'    --step  Set sliding window step. [Default: 500]
        '-j'    --threads       Set the maximum number of threads. [Default: max. available]
        '-l'    --min-block-length      Set minimum block length for evaluation. [Default: 2000]
        '-d'    --max-block-distance    Set maximum block distance for merging. [Default: 200]

Optional Parameters:
        '-m'    --mode  Set analysis modes, separate them by commas. [Options: all,match,gc,entropy]
        '-k'    --keep-window-data      Keep window data for analysis, memory-intensive. [Default: false]
        '-v'    --version       Print current software version.
        '-h'    --help  Print current software options.
        --verbose       verbose output.
```

Outputs
------------
Teloscope outputs telomere annotations in BED files. All the outputs are: 
* `telomere_blocks_all.bed` Annotation of the full telomere in the assembly. This is made of canonical and non-canonical repeats. 
* `telomere_blocks_canonical.bed` Blocks of adjacent canonical repeat matches. Outside of the ends, it represents interstitial telomeres (ITSs).
* `window_metrics.tsv` Tabulated file with calculated window metrics such as GC% and Shannon Entropy
* `window_repeats.bedgraph` File with canonical repeats, non-canonical repeats, canonical densities, and non-canonical densities by window. 
* `canonical_matches.bed` Coordinates of canonical repeats throughout the assembly. 
* `noncanonical_matches.bed` Coordinates of non-canonical repeats in terminal regions of contigs. 

How it works
------------
Briefly, **Teloscope** reads an assembly and decomposes its parts. It uses prefix trees and sliding windows to efficiently perform multiple string matching and counting of telomeric repeats. It analyzes the informational properties of the sliding windows to find telomeric blocks. These blocks are collected, post-processed, and filtered according to their positional and conformational properties. 

How to cite
------------
If you use **Teloscope** in your research, please, cite this repository. 
https://github.com/vgl-hub/teloscope/
