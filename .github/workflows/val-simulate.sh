#!/bin/bash
set -e

SIMULATE=build/bin/teloscope-simulate
TELOSCOPE=build/bin/teloscope
OUTDIR=testFiles/simulate
N=${SIM_N:-1000}
SEED=${SIM_SEED:-42}

# Header
echo -e "rate\ttips\tdetected\tsensitivity\tbias_bp\tabs_err_bp\ttvr_rate\tfp_blocks"

for rate in 1e-6 1e-5 1e-4 1e-3 1e-2; do
    DIR="$OUTDIR/rate_${rate}"
    TELOUT="$DIR/teloscope_out"
    mkdir -p "$TELOUT"

    $SIMULATE -n "$N" -r "$rate" -s "$SEED" -o "$DIR" 2>/dev/null
    $TELOSCOPE -f "$DIR/sequences.fa" -o "$TELOUT" 2>/dev/null >/dev/null

    BED="$TELOUT/sequences.fa_terminal_telomeres.bed"
    RESULT=$($SIMULATE --evaluate -g "$DIR/ground_truth.tsv" -b "$BED")
    echo -e "${rate}\t${RESULT}"
done
