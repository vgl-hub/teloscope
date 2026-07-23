# Teloscope
[<img alt="github" src="https://img.shields.io/badge/github-vgl--hub/teloscope-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/vgl-hub/teloscope)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-teloscope-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/teloscope/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2025.10.14.682431.svg)](https://doi.org/10.1101/2025.10.14.682431)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/version.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/platforms.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/license.svg)](https://anaconda.org/bioconda/teloscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/teloscope/badges/downloads.svg)](https://anaconda.org/bioconda/teloscope)

Teloscope scans assembly ends for telomeric repeats. It reads `FASTA`, `FASTA.gz`, and `GFA` inputs, merges repeat matches into telomere blocks, classifies scaffolds in FASTA mode, and writes files that are easy to inspect in BED, TSV, BEDgraph, PDF, or GFA form. It can also subset telomeric FASTQ reads and BAM records.

In FASTA mode, Teloscope writes terminal telomere annotations, gap coordinates, and a summary table. In GFA mode, it writes an annotated graph for BandageNG. Synthetic telomere nodes attach to the assembly with `L` links at `0M` overlap, the direct adjacency a cap represents, so BandageNG draws them as caps. `J` jump records stay reserved for real assembly gaps.

## What Teloscope writes

- FASTA mode: `*_terminal_telomeres.bed`, `*_gaps.bed`, `*_report.tsv`, plus optional window tracks, match BED files, ITS BED files, and a PDF report.
- GFA mode: `<input>.telo.annotated.gfa` with telomere placeholder segments linked to the original graph, plus `<input>.telo.annotated.colors.csv` that paints the caps green for BandageNG.
- FASTQ subset mode: unchanged passing FASTQ records on stdout.
- BAM subset mode: a valid BAM stream containing the original header and unchanged passing alignment records.

## Install

From source:

```sh
git clone https://github.com/vgl-hub/teloscope.git --recursive
cd teloscope
make -j
```

To build the helper binaries as well:

```sh
make all -j
```

From Bioconda:

```sh
conda install -c bioconda teloscope
```

Build requirements:

- C++17 compiler
- `zlib`
- `pthread`
- the `gfalibs` submodule

If the submodule is missing, run:

```sh
git submodule update --init --recursive
```

For `--plot-report`, install Python 3 with `matplotlib`, `numpy`, and `pandas`.

## Quick start

| Task | Command |
| --- | --- |
| Scan a vertebrate assembly with the default motif | `teloscope asm.fa` |
| Read compressed FASTA directly | `teloscope asm.fa.gz` |
| Write all optional FASTA outputs | `teloscope asm.fa -o results/ -r -g -e -m -i --plot-report` |
| Keep one inspected submitter naming family | `teloscope asm.fa --include-prefix hap1_chr -o results/` |
| Keep exact sequence IDs from a BED/list file (safest for database assemblies) | `teloscope asm.fa --include-bed chromosomes.ids -o results/` |
| Switch to a plant canonical repeat | `teloscope asm.fa -c CCCTAAA` |
| Search explicit motif variants | `teloscope asm.fa -c TTAGGG -p TTAGGG,TCAGGG,TGAGGG,TTGGGG` |
| Annotate a graph for BandageNG | `teloscope asm.gfa -o results/` |
| Subset telomeric HiFi reads before mapping | `teloscope --fastq-subset reads.fq.gz -j 32 \| minimap2 -ax map-hifi ref.fa -` |
| Subset telomeric records from BAM | `teloscope --bam-subset reads.bam -j 32 > telomeric.bam` |
| Read decompressed stdin | `zcat asm.fa.gz | teloscope -o results/` |

Notes:

- Teloscope always searches both each input pattern and its reverse complement.
- If `-p` is omitted, Teloscope derives the search set from `-c`.
- Any genome-wide output flag (`-r`, `-g`, `-e`, `-m`, `-i`) disables ultra-fast mode automatically.
- GFA mode attaches telomere caps with `L` links at `0M` overlap; `J` records stay reserved for real assembly gaps.
- Gzipped stdin is not supported. Decompress before piping.
- `--fastq-subset` writes FASTQ to stdout and diagnostics to stderr. Pass `-o` to save the reads to a file instead of streaming them.
- `--bam-subset` writes BAM to stdout and diagnostics to stderr. Pass `-o` to write `<input_stem>_telomeric.bam`.
- Read subset modes use a `42` bp default minimum block length; assembly annotation keeps the `300` bp default. Use `-l` to override either mode.
- BAM support has no external bioinformatics runtime dependency: it uses `zlib` directly and does not require HTSlib, `samtools`, or another converter.

## Filter assembly records

Record filtering is off by default, so existing FASTA and GFA behavior is unchanged. When filtering is active, the selection domain is every FASTA record, every supported GFA1 `P` path, or every segment in a pathless GFA1 graph.

| Flag | Value | Default | Effect |
| --- | --- | --- | --- |
| `--include-bed FILE` | BED or one-ID-per-line file | not set | keep records whose complete primary ID occurs in column 1 |
| `--exclude-bed FILE` | BED or one-ID-per-line file | not set | remove records whose complete primary ID occurs in column 1 |
| `--include-prefix LIST` | comma-separated literal prefixes | not set | keep records whose primary ID starts with any listed prefix |
| `--exclude-prefix LIST` | comma-separated literal prefixes | not set | remove records whose primary ID starts with any listed prefix |

The BED pair uses the familiar gfastats flag names for an auditable exact-ID workflow, but Teloscope deliberately treats column 1 as a whole-record selector. This avoids coordinate edits that could change terminal context. The options are long-only because Teloscope already uses `-i` for ITS output and `-e` for entropy output. The prefix pair handles a known accession family without forcing users to generate a list or enabling error-prone regex syntax.

All four flags can be repeated. If no include flag is set, every input record starts selected. If one or more include flags are set, their matches form one union. All exclude matches are then removed from that union, so exclusion wins. Exact IDs and prefixes are case-sensitive. Prefixes are literal text, not regular expressions or globs: use `--include-prefix sample_`, not `sample_*`. Spaces around comma-separated prefixes are trimmed; an empty prefix is an error.

For FASTA, the primary ID is the first whitespace-delimited token after `>`. An NCBI header such as `>NC_000001.11 Homo sapiens chromosome 1` therefore has the ID `NC_000001.11`. Matching includes the accession version. For GFA1 with `P` records, selectors match path names. A pathless GFA1 graph uses segment names. Filtering does not remove supported original topology or excluded paths; it limits the terminal segment ends that Teloscope scans.

Each `--include-bed` or `--exclude-bed` file may mix one-column ID rows and BED3+ rows. Blank lines, `#` comments, `track` lines, and `browser` lines are ignored. BED start and end fields must be unsigned integers with start no greater than end. Coordinates are validated but are not used to crop, splice, or remap a sequence: selection always applies to the whole record.

Every exact ID and every prefix must match at least one ID in the input domain. Teloscope exits on an unmatched selector, duplicate FASTA primary ID, invalid selector file, or an empty final selection. It prints the selected and input counts to stderr. FASTA summary output also records both counts, and all BED, BEDgraph, TSV, and plot-report outputs contain selected records only. Filtering happens after the assembly is loaded, so it reduces scanning and output size but not input parsing or peak loader memory.

These four flags apply only to FASTA and supported GFA1 assembly files. Uncompressed filtered stdin is treated as FASTA; filtered GFA input must use a `.gfa` or `.gfa.gz` filename. Filtered GFA2, GFA1 `C` containment and `W` walk records, and unknown GFA record types are rejected because the current graph writer cannot preserve them safely. Teloscope also rejects filtered FASTQ/BAM input and rejects the flags with `--fastq-subset` or `--bam-subset`.

### Database FASTA names

Do not treat one prefix as a universal chromosome label. In NCBI assembly packages, `GCA_...` and `GCF_...` identify the assembly; they are not the sequence IDs at the start of `*_genomic.fna` headers. Those FASTA IDs are sequence accession.version values. RefSeq families such as `NC_`, `NT_`, `NW_`, and `NZ_` describe broad record classes. GenBank/INSDC WGS accessions use project-specific four- or six-letter series, and NCBI classifies `CM` as a scaffold/CON accession family rather than a universal chromosome class. Ensembl, UCSC, and submitter FASTAs use other schemes.

For a rigorous chromosome-only set, use the role and molecule type in the NCBI [genome sequence report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/) or [assembly report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/), then pass the matching sequence accessions as an exact list. In an assembly report, column 4 distinguishes chromosomes from mitochondria and other assembled molecules; column 5 matches a GenBank/GCA genomic FASTA, while column 7 matches a paired RefSeq/GCF genomic FASTA:

```sh
# GenBank/GCA FASTA IDs
awk -F '\t' '$1 !~ /^#/ && $2 == "assembled-molecule" && $4 == "Chromosome" {print $5}' \
  assembly_report.txt > chromosomes.ids

# RefSeq/GCF FASTA IDs
awk -F '\t' '$1 !~ /^#/ && $2 == "assembled-molecule" && $4 == "Chromosome" && $7 != "na" {print $7}' \
  assembly_report.txt > chromosomes.ids

teloscope assembly_genomic.fna.gz --include-bed chromosomes.ids -o results/
```

Use `--include-prefix` when you have inspected the actual headers and a prefix cleanly represents the set you want. NCBI’s [accession-prefix table](https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/) lists INSDC formats, but the assembly/sequence role remains the safer chromosome classifier.

## Typical output layout

FASTA run:

```text
results/
  asm.fa_terminal_telomeres.bed
  asm.fa_gaps.bed
  asm.fa_report.tsv
  asm.fa_window_repeat_density.bedgraph
  asm.fa_window_canonical_ratio.bedgraph
  asm.fa_window_strand_ratio.bedgraph
  asm.fa_plot_report.pdf
```

GFA run:

```text
results/
  asm.gfa.telo.annotated.gfa
  asm.gfa.telo.annotated.colors.csv
```

## Documentation

| Page | Covers |
| --- | --- |
| [Parameters](docs/parameters.md) | command-line flags, defaults, and flag interactions |
| [Outputs](docs/outputs.md) | every output file, naming rules, and column layouts |
| [Classification](docs/classification.md) | FASTA scaffold classes and the `granular` labels |
| [Algorithm](docs/algorithm.md) | how FASTA mode and GFA mode are processed |
| [Report generation](docs/report.md) | `--plot-report`, ITS plotting, standalone plotting, and report inputs |
| [Simulation](docs/simulation.md) | the synthetic benchmark generator and evaluator |
| [Testing](docs/testing.md) | validator runs, report regression checks, and test regeneration |
| [Troubleshooting](docs/troubleshooting.md) | common build, input, and runtime failures |
| [Release checklist](docs/release.md) | GitHub, Bioconda, and Zenodo release steps |
| [Validation format](validateFiles/README.md) | the `.tst` harness, including directive-mode GFA cases |

## Repo layout

- `src/` and `include/`: main C++ implementation and headers
- `scripts/`: Python and shell helpers, including `teloscope_report.py`, `plot_its.py`, and plotting regression checks
- `docs/`: user-facing documentation
- `testFiles/`: public FASTA and GFA fixtures plus expected outputs
- `validateFiles/`: `.tst` manifests used by `teloscope-validate`
- `gfalibs/`: graph I/O submodule used for GFA parsing and writing

## Validation

Build the validator and run the checked-in suite:

```sh
make validate
build/bin/teloscope-validate validateFiles
```

Run the report layout regression script:

```sh
python3 scripts/test_teloscope_report.py
```

Run the gap BED regression script:

```sh
bash scripts/test_gaps_bed.sh
```

## Citation

Teloscope is part of the gfastar tool suite.

If you use Teloscope, cite:

The complete genome of a songbird  
Giulio Formenti, Nivesh Jain, Jack A. Medico, Marco Sollitto, Dmitry Antipov, Suziane Barcellos, Matthew Biegler, Ines Borges, J King Chang, Ying Chen, Haoyu Cheng, Helena Conceicao, Matthew Davenport, Lorraine De Oliveira, Erick Duarte, Gillian Durham, Jonathan Fenn, Niamh Forde, Pedro A. Galante, Kenji Gerhardt, Alice M. Giani, Simona Giunta, Juhyun Kim, Aleksey Komissarov, Bonhwang Koo, Sergey Koren, Denis Larkin, Chul Lee, Heng Li, Kateryna Makova, Patrick Masterson, Terence Murphy, Kirsty McCaffrey, Rafael L.V. Mercuri, Yeojung Na, Mary J. O'Connell, Shujun Ou, Adam Phillippy, Marina Popova, Arang Rhie, Francisco J. Ruiz-Ruano, Simona Secomandi, Linnea Smeds, Alexander Suh, Tatiana Tilley, Niki Vontzou, Paul D. Waters, Jennifer Balacco, Erich D. Jarvis  
bioRxiv 2025.10.14.682431  
https://doi.org/10.1101/2025.10.14.682431
