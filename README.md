# Teloscope
A telomere annotation tool.

Teloscope is a fast and comprehensive tool for matching, counting and reporting telomeric pattens from genome assemblies (.fa) or (.fa.gz). Teloscope allows the calculation of other metrics such as: 

* GC
* Shannon Entropy

Teloscope reports all these metrics and the final telomere coordinates in BED/BEDgraph files and produces a summary report.

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/teloscope.git --recursive` and `make -j` in `teloscope` folder.

## Usage

`teloscope -f input.[fasta][.gz] -p TTAGGG,TTAGGGG -w [window size] -s [step size]`

To check out all options and flags, please use:
`teloscope -h`

**Note:** Teloscope automatically explores for the input patterns and their reverse complements. If none are provided, it will scan for the canonical CCCTAA/TTAGGG repeats. 

## Description

Briefly, **Teloscope** reads any assembly and decompose its parts. It uses prefix trees and sliding windows to match and count telomeric patterns efficiently. It also analyzes the informational properties of the assembly to distinguish canonical and non-canonical telomere repeats.

## How to cite

If you use **Teloscope**, please cite the current repository. 