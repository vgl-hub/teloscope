#!/bin/bash
# Test gap BED output: verifies *_gaps.bed against expected files.
# All tests should FAIL before gap output is implemented,
# and PASS after.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TELO="$REPO_DIR/build/bin/teloscope"
TMP_DIR="$REPO_DIR/testFiles/tmp"
EXPECTED_DIR="$REPO_DIR/testFiles/expected"

PASS=0
FAIL=0
TOTAL=0

red()   { printf "\033[0;31m%s\033[0m" "$1"; }
green() { printf "\033[0;32m%s\033[0m" "$1"; }

run_test() {
    local fa_file="$1"       # e.g. gapped_t2t.fa
    local expected="$2"      # e.g. testFiles/expected/gapped_t2t.fa_gaps.bed
    local desc="$3"
    TOTAL=$((TOTAL + 1))

    rm -rf "$TMP_DIR"
    mkdir -p "$TMP_DIR"

    # Run teloscope
    "$TELO" -f "$REPO_DIR/testFiles/$fa_file" -o "$TMP_DIR" >/dev/null 2>&1 || true

    # Find the gaps.bed file
    local gaps_file
    gaps_file=$(find "$TMP_DIR" -name "*_gaps.bed" 2>/dev/null | head -1)

    if [ -z "$gaps_file" ]; then
        FAIL=$((FAIL + 1))
        echo "  $(red FAIL) $desc — _gaps.bed not found"
        return
    fi

    # Compare content (ignore trailing newlines)
    if diff <(sed '/^$/d' "$gaps_file") <(sed '/^$/d' "$expected") >/dev/null 2>&1; then
        PASS=$((PASS + 1))
        echo "  $(green PASS) $desc"
    else
        FAIL=$((FAIL + 1))
        echo "  $(red FAIL) $desc — content mismatch:"
        diff <(sed '/^$/d' "$gaps_file") <(sed '/^$/d' "$expected") || true
    fi
}

# Check binary exists
if [ ! -x "$TELO" ]; then
    echo "Error: teloscope binary not found at $TELO"
    echo "Run 'make' first."
    exit 1
fi

echo "Running gap BED output tests..."
echo ""

# --- Gapped files: should produce non-empty _gaps.bed ---
echo "Single-gap test cases:"
run_test "gapped_t2t.fa"            "$EXPECTED_DIR/gapped_t2t.fa_gaps.bed"            "gapped_t2t (gap@1500-1600)"
run_test "gapped_none.fa"           "$EXPECTED_DIR/gapped_none.fa_gaps.bed"           "gapped_none (gap@1000-1100)"
run_test "gapped_incomplete.fa"     "$EXPECTED_DIR/gapped_incomplete.fa_gaps.bed"     "gapped_incomplete (gap@1600-1700)"
run_test "gapped_incomplete_q.fa"   "$EXPECTED_DIR/gapped_incomplete_q.fa_gaps.bed"   "gapped_incomplete_q (gap@2000-2100)"
run_test "gapped_discordant.fa"     "$EXPECTED_DIR/gapped_discordant.fa_gaps.bed"     "gapped_discordant (gap@2000-2100)"
run_test "gapped_discordant_q.fa"   "$EXPECTED_DIR/gapped_discordant_q.fa_gaps.bed"   "gapped_discordant_q (gap@2600-2700)"
run_test "gapped_misassembly.fa"    "$EXPECTED_DIR/gapped_misassembly.fa_gaps.bed"    "gapped_misassembly (gap@3700-3800)"
run_test "gapped_misassembly_qq.fa" "$EXPECTED_DIR/gapped_misassembly_qq.fa_gaps.bed" "gapped_misassembly_qq (gap@500-600)"

echo ""
echo "Multi-gap test case:"
run_test "multi_gap_t2t.fa"         "$EXPECTED_DIR/multi_gap_t2t.fa_gaps.bed"         "multi_gap_t2t (gaps@1000-1050,1450-1500)"

echo ""
echo "No-gap test cases (expect empty _gaps.bed):"
run_test "t2t.fa"                   "$EXPECTED_DIR/t2t.fa_gaps.bed"                   "t2t (no gaps)"
run_test "no_telo.fa"               "$EXPECTED_DIR/no_telo.fa_gaps.bed"               "no_telo (no gaps)"
run_test "multi.fa"                 "$EXPECTED_DIR/multi.fa_gaps.bed"                 "multi (no gaps)"

echo ""
echo "Results: $PASS passed, $FAIL failed (out of $TOTAL)"

rm -rf "$TMP_DIR"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
