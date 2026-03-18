[← Back to README](../README.md)

Simulation and validation
============

Teloscope includes a simulation framework for benchmarking detection accuracy across mutation rates. Build with:
```sh
make simulate
```

This produces `teloscope-simulate`, a standalone tool with two modes:

**Generate mode** produces synthetic multi-contig FASTA files with known telomere positions:
```sh
teloscope-simulate -n 1000 -r 1e-4 -s 42 -o testFiles/simulate/rate_1e-4
```

Each synthetic sequence (~50 kbp) is built from five segments:

1. **p-arm canonical** 12 kbp or 2000 exact copies of CCCTAA.
2. **p-arm TVR** 600 bp or 100 copies of CCCTAA, each with one random substitution.
3. **Central segment** 25 kbp of random nucleotides.
4. **q-arm TVR** 600 bp or 100 copies of TTAGGG, each with one random substitution.
5. **q-arm canonical** 12 kbp or 2000 exact copies of TTAGGG.

After construction, per-base mutations are applied across the entire sequence at the specified rate (`-r`). The TVR regions have exactly one mismatch per repeat, so they are detectable at `-x 1` but not at `-x 0`. This tests both default and edit-distance modes.

Outputs: `sequences.fa` (synthetic FASTA) and `ground_truth.tsv` (per-contig telomere coordinates with exact boundaries for each segment).

**Evaluate mode** compares teloscope BED output against the ground truth:
```sh
teloscope-simulate --evaluate -g ground_truth.tsv -b terminal_telomeres.bed
```

Output columns: `total_tips`, `detected`, `sensitivity`, `mean_bias_bp`, `mean_abs_err_bp`, `tvr_rate`, `fp_blocks`.

### Simulation parameters

| Flag | Description | Default |
|------|-------------|---------|
| `-n` | Number of sequences to generate | `1000` |
| `-R` | Canonical repeats per tip (each repeat is 6 bp) | `2000` |
| `-T` | TVR (variant) repeats per tip (each repeat is 6 bp with 1 substitution) | `100` |
| `-L` | Central segment length in bp | `25000` |
| `-r` | Per-base mutation rate applied to the full sequence | `0.0` |
| `-s` | Random seed for reproducibility | `42` |
| `-o` | Output directory | `testFiles/simulate` |

### Pipeline script

A sweep across mutation rates (1e-6 through 1e-2) can be run with:
```sh
bash .github/workflows/val-simulate.sh
```

Configurable via environment: `SIM_N=1000000 SIM_SEED=42 bash .github/workflows/val-simulate.sh`
