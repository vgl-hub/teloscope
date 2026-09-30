[Back to README](index.md)

# Simulation

`teloscope-simulate` builds synthetic assemblies with known telomeres and scores Teloscope's calls against them.

```sh
make simulate    # writes build/bin/teloscope-simulate
```

## Generate

```sh
build/bin/teloscope-simulate -n 1000 -r 1e-4 -s 42 -o testFiles/simulate/rate_1e-4
```

This writes `sequences.fa` and `ground_truth.tsv`. Each sequence runs canonical telomere, TVR block, random sequence, TVR block, canonical telomere; the mutation rate then applies to all of it.

| Flag | Meaning | Default |
| --- | --- | --- |
| `-n` | number of sequences | `1000` |
| `-R` | canonical repeats per telomere | `2000` |
| `-T` | TVR repeats per telomere | `100` |
| `-L` | random internal length in bp | `25000` |
| `-r` | per-base mutation rate | `0.0` |
| `-s` | random seed | `42` |
| `-o` | output directory | `testFiles/simulate` |

## Evaluate

```sh
build/bin/teloscope -f testFiles/simulate/rate_1e-4/sequences.fa -o testFiles/simulate/rate_1e-4/teloscope_out
build/bin/teloscope-simulate --evaluate \
  -g testFiles/simulate/rate_1e-4/ground_truth.tsv \
  -b testFiles/simulate/rate_1e-4/teloscope_out/sequences.fa_terminal_telomeres.bed
```

The evaluator prints the seven values on one tab-separated line, without a header:

| Field | Meaning |
| --- | --- |
| `total_tips` | two per sequence |
| `detected` | tips overlapped by a row with the right label: `p` or `b` at the start, `q` or `b` at the end |
| `sensitivity` | `detected / total_tips` |
| `mean_bias_bp` | mean of called minus true canonical length |
| `mean_abs_err_bp` | mean absolute length error |
| `tvr_rate` | share of detected tips whose row reaches into the TVR |
| `fp_blocks` | rows that overlap no true telomere or TVR |

## Sweep

`.github/workflows/val-simulate.sh` runs rates from `1e-6` to `1e-2` and prints one line per rate. `SIM_N` and `SIM_SEED` override its 1000 sequences and seed 42:

```sh
SIM_N=10000 SIM_SEED=7 bash .github/workflows/val-simulate.sh
```
