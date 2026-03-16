#!/bin/bash
set -e

# Run standard .tst validation
build/bin/teloscope-validate validateFiles

# File count tests
PASS="\033[0;32mPASS\033[0m"
FAIL="\033[0;31mFAIL\033[0m"
EXIT=0

check_file_count() {
    local desc="$1"
    local expected="$2"
    local dir="$3"
    local actual
    actual=$(ls -1 "$dir" | wc -l)
    if [ "$actual" -eq "$expected" ]; then
        echo -e "$PASS $desc (expected $expected files, got $actual)"
    else
        echo -e "$FAIL $desc (expected $expected files, got $actual)"
        ls -1 "$dir"
        EXIT=1
    fi
}

TMPDIR="testFiles/tmp"
mkdir -p "$TMPDIR"

# Test: default mode produces 1 file (terminal_telomeres.bed)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/bTaeGut7_chr33_mat.fa.gz -o "$TMPDIR" 2>/dev/null >/dev/null
check_file_count "default mode file count" 1 "$TMPDIR"

# Test: -r produces 4 files (terminal + density + canonical_ratio + strand_ratio bedgraphs)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/bTaeGut7_chr33_mat.fa.gz -o "$TMPDIR" -r 2>/dev/null >/dev/null
check_file_count "-r flag file count" 4 "$TMPDIR"

# Test: -g produces 2 files (terminal + gc bedgraph)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/bTaeGut7_chr33_mat.fa.gz -o "$TMPDIR" -g 2>/dev/null >/dev/null
check_file_count "-g flag file count" 2 "$TMPDIR"

# Test: -e produces 2 files (terminal + entropy bedgraph)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/bTaeGut7_chr33_mat.fa.gz -o "$TMPDIR" -e 2>/dev/null >/dev/null
check_file_count "-e flag file count" 2 "$TMPDIR"

# Test: -r -g -e produces 6 files (terminal + 3 repeat + gc + entropy)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/bTaeGut7_chr33_mat.fa.gz -o "$TMPDIR" -r -g -e 2>/dev/null >/dev/null
check_file_count "-r -g -e flag file count" 6 "$TMPDIR"

# Test: positional argument works like -f
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope testFiles/bTaeGut7_chr33_mat.fa.gz -o "$TMPDIR" 2>/dev/null >/dev/null
check_file_count "positional arg file count" 1 "$TMPDIR"

rm -rf "$TMPDIR"/* 2>/dev/null || true
exit $EXIT
