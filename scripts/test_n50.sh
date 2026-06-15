#!/bin/bash
# Test Scaffold N50 / Contig N50 in the Assembly Summary Report.
# Crafts a FASTA with known scaffold/contig structure and checks the reported N50s.
#
#   scaf1: 1000 bp + (50 N gap) + 600 bp   -> pathSize 1650, contigs {1000, 600}
#   scaf2: 800 bp                          -> pathSize 800,  contig  {800}
#
#   scaffolds {1650, 800} (total 2450)  -> Scaffold N50 = 1650
#   contigs   {1000, 800, 600} (total 2400) -> Contig N50 = 800
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TELO="$REPO_DIR/build/bin/teloscope"
TMP_DIR="$REPO_DIR/testFiles/tmp_n50"

red()   { printf "\033[0;31m%s\033[0m" "$1"; }
green() { printf "\033[0;32m%s\033[0m" "$1"; }

if [ ! -x "$TELO" ]; then
    echo "Error: teloscope binary not found at $TELO (run 'make' first)"; exit 1
fi

rm -rf "$TMP_DIR"; mkdir -p "$TMP_DIR"
FA="$TMP_DIR/n50_basic.fa"
{
    echo ">scaf1"
    awk 'BEGIN{for(i=0;i<1000;i++)printf "A"; for(i=0;i<50;i++)printf "N"; for(i=0;i<600;i++)printf "A"; print ""}'
    echo ">scaf2"
    awk 'BEGIN{for(i=0;i<800;i++)printf "A"; print ""}'
} > "$FA"

OUT=$("$TELO" -f "$FA" -o "$TMP_DIR/out" 2>/dev/null) || true

scaf=$(printf '%s\n' "$OUT" | awk -F'\t' '/^Scaffold N50:/{print $2; exit}')
ctg=$(printf '%s\n' "$OUT"  | awk -F'\t' '/^Contig N50:/{print $2; exit}')

EXP_SCAF=1650
EXP_CTG=800
fail=0

if [ "${scaf:-}" = "$EXP_SCAF" ]; then
    echo "  $(green PASS) Scaffold N50 = $scaf"
else
    echo "  $(red FAIL) Scaffold N50 = '${scaf:-<missing>}' (expected $EXP_SCAF)"; fail=1
fi
if [ "${ctg:-}" = "$EXP_CTG" ]; then
    echo "  $(green PASS) Contig N50 = $ctg"
else
    echo "  $(red FAIL) Contig N50 = '${ctg:-<missing>}' (expected $EXP_CTG)"; fail=1
fi

rm -rf "$TMP_DIR"
[ "$fail" -eq 0 ] && { echo "test_n50: all passed"; exit 0; } || { echo "test_n50: FAILED"; exit 1; }
