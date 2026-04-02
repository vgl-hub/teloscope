[Back to README](../README.md)

# Simulation and validation

Teloscope includes a synthetic benchmark generator and evaluator in `teloscope-simulate`.

Build it with:

```sh
make simulate
```

The binary is written to `build/bin/teloscope-simulate`.

## Generate synthetic assemblies

Example:

```sh
build/bin/teloscope-simulate -n 1000 -r 1e-4 -s 42 -o testFiles/simulate/rate_1e-4
```

This writes:

- `sequences.fa`
- `ground_truth.tsv`

Each synthetic sequence contains:

1. left canonical telomere
2. left TVR block
3. random internal sequence
4. right TVR block
5. right canonical telomere

After the sequence is built, the requested mutation rate is applied across the whole sequence.

## Evaluate Teloscope output

Run Teloscope on the synthetic FASTA, then compare the called terminal BED file to the ground truth:

```sh
build/bin/teloscope -f testFiles/simulate/rate_1e-4/sequences.fa -o testFiles/simulate/rate_1e-4/teloscope_out
build/bin/teloscope-simulate --evaluate \
  -g testFiles/simulate/rate_1e-4/ground_truth.tsv \
  -b testFiles/simulate/rate_1e-4/teloscope_out/sequences.fa_terminal_telomeres.bed
```

The evaluator prints:

- `total_tips`
- `detected`
- `sensitivity`
- `mean_bias_bp`
- `mean_abs_err_bp`
- `tvr_rate`
- `fp_blocks`

## Generator parameters

| Flag | Meaning | Default |
| --- | --- | --- |
| `-n` | number of sequences | `1000` |
| `-R` | canonical repeats per telomere end | `2000` |
| `-T` | TVR repeats per telomere end | `100` |
| `-L` | internal random segment length in bp | `25000` |
| `-r` | per-base mutation rate | `0.0` |
| `-s` | random seed | `42` |
| `-o` | output directory | `testFiles/simulate` |

## Sweep script

The repo includes a shell wrapper for mutation-rate sweeps:

```sh
bash .github/workflows/val-simulate.sh
```

Optional environment overrides:

```sh
SIM_N=1000000 SIM_SEED=42 bash .github/workflows/val-simulate.sh
```
