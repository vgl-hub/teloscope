#!/bin/bash
set -e

# Run standard .tst validation (53 tests)
build/bin/teloscope-validate validateFiles

# File count and output pattern tests
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

check_output_contains() {
    local desc="$1"
    local pattern="$2"
    local output="$3"
    if echo "$output" | grep -qF "$pattern"; then
        echo -e "$PASS $desc"
    else
        echo -e "$FAIL $desc (pattern '$pattern' not found)"
        EXIT=1
    fi
}

check_output_not_contains() {
    local desc="$1"
    local pattern="$2"
    local output="$3"
    if echo "$output" | grep -qF "$pattern"; then
        echo -e "$FAIL $desc (pattern '$pattern' should not appear)"
        EXIT=1
    else
        echo -e "$PASS $desc"
    fi
}

TMPDIR="testFiles/tmp"
mkdir -p "$TMPDIR"

# =====================================================
# File count tests (output flag combinations)
# =====================================================

# Test: default mode produces 2 files (terminal_telomeres.bed + report.tsv)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" 2>/dev/null >/dev/null
check_file_count "default mode file count" 3 "$TMPDIR"

# Test: -r produces 5 files (terminal + report + density + canonical_ratio + strand_ratio)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -r 2>/dev/null >/dev/null
check_file_count "-r flag file count" 6 "$TMPDIR"

# Test: -g produces 3 files (terminal + report + gc)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -g 2>/dev/null >/dev/null
check_file_count "-g flag file count" 4 "$TMPDIR"

# Test: -e produces 3 files (terminal + report + entropy)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -e 2>/dev/null >/dev/null
check_file_count "-e flag file count" 4 "$TMPDIR"

# Test: -m produces 4 files (terminal + report + canonical_matches + noncanonical_matches)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -m 2>/dev/null >/dev/null
check_file_count "-m flag file count" 5 "$TMPDIR"

# Test: -i produces 3 files (terminal + report + interstitial_telomeres)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -i 2>/dev/null >/dev/null
check_file_count "-i flag file count" 4 "$TMPDIR"

# Test: -r -g -e produces 7 files
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -r -g -e 2>/dev/null >/dev/null
check_file_count "-r -g -e file count" 8 "$TMPDIR"

# Test: -r -m -g -e -i produces 10 files (terminal + report + 3 repeat + 2 match + its + gc + entropy)
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/t2t.fa -o "$TMPDIR" -r -m -g -e -i 2>/dev/null >/dev/null
check_file_count "-r -m -g -e -i file count" 11 "$TMPDIR"

# Test: positional argument works like -f
rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope testFiles/t2t.fa -o "$TMPDIR" 2>/dev/null >/dev/null
check_file_count "positional arg file count" 3 "$TMPDIR"

# =====================================================
# Classification consistency: default vs -r/-m/-i
# The window loop fix ensures these all agree
# =====================================================

T2T_DEFAULT=$(build/bin/teloscope -f testFiles/t2t.fa 2>/dev/null)
T2T_R=$(build/bin/teloscope -f testFiles/t2t.fa -r -o "$TMPDIR" 2>/dev/null)
T2T_M=$(build/bin/teloscope -f testFiles/t2t.fa -m -o "$TMPDIR" 2>/dev/null)

check_output_contains "t2t default: type=t2t" "t2t	PQ" "$T2T_DEFAULT"
check_output_contains "t2t -r: type=t2t" "t2t	PQ" "$T2T_R"
check_output_contains "t2t -m: type=t2t" "t2t	PQ" "$T2T_M"

# =====================================================
# Classification tests on synthetic data
# =====================================================

# All 10 scaffold types
OUT=$(build/bin/teloscope -f testFiles/t2t.fa 2>/dev/null)
check_output_contains "t2t classification" "t2t	PQ" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/incomplete_p.fa 2>/dev/null)
check_output_contains "incomplete_p classification" "incomplete	P" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/incomplete_q.fa 2>/dev/null)
check_output_contains "incomplete_q classification" "incomplete	Q" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/no_telo.fa 2>/dev/null)
check_output_contains "no_telo classification" "none" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/misassembly.fa 2>/dev/null)
check_output_contains "misassembly Pp classification" "misassembly	Pp" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/misassembly_qq.fa 2>/dev/null)
check_output_contains "misassembly Qq classification" "misassembly	Qq" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/discordant.fa 2>/dev/null)
check_output_contains "discordant classification" "discordant	P*" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/gapped_t2t.fa 2>/dev/null)
check_output_contains "gapped_t2t classification" "gapped_t2t	PQ" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/gapped_misassembly.fa 2>/dev/null)
check_output_contains "gapped_misassembly classification" "gapped_misassembly	Pp" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/gapped_incomplete.fa 2>/dev/null)
check_output_contains "gapped_incomplete classification" "gapped_incomplete	P" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/gapped_none.fa 2>/dev/null)
check_output_contains "gapped_none classification" "gapped_none" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/gapped_discordant.fa 2>/dev/null)
check_output_contains "gapped_discordant classification" "gapped_discordant	P*" "$OUT"

# =====================================================
# Spelling consistency
# =====================================================

OUT=$(build/bin/teloscope -f testFiles/t2t.fa 2>/dev/null)
check_output_contains "summary has Misassembled" "Misassembled:" "$OUT"
check_output_contains "summary has Gapped misassembled" "Gapped misassembled:" "$OUT"
check_output_not_contains "no Missassembled typo" "Missassembled:" "$OUT"

# =====================================================
# Flag behavior tests
# =====================================================

# -l flag: raising min block length kills short blocks
OUT=$(build/bin/teloscope -f testFiles/t2t.fa -l 1000 2>/dev/null)
check_output_contains "-l 1000 kills 600bp blocks" "none" "$OUT"

# -t flag: small terminal limit misses telomeres
OUT=$(build/bin/teloscope -f testFiles/t2t.fa -t 100 2>/dev/null)
check_output_contains "-t 100 misses telomeres" "none" "$OUT"

# -y flag: density threshold controls block survival
OUT=$(build/bin/teloscope -f testFiles/density_edge.fa -y 0.8 2>/dev/null)
check_output_contains "-y 0.8 kills 50% density block" "incomplete	Q" "$OUT"

OUT=$(build/bin/teloscope -f testFiles/density_edge.fa -y 0.3 2>/dev/null)
check_output_contains "-y 0.3 keeps 50% density block" "t2t	PQ" "$OUT"

# -k flag: small merge distance prevents block formation
OUT=$(build/bin/teloscope -f testFiles/t2t.fa -k 1 2>/dev/null)
check_output_contains "-k 1 prevents block formation" "none" "$OUT"

# Multi-contig counts
OUT=$(build/bin/teloscope -f testFiles/multi.fa 2>/dev/null)
check_output_contains "multi: 3 paths" "Total paths:	3" "$OUT"
check_output_contains "multi: 3 telomeres" "Total telomeres:	3" "$OUT"
check_output_contains "multi: T2T count" "T2T:	1" "$OUT"
check_output_contains "multi: Incomplete count" "Incomplete:	1" "$OUT"
check_output_contains "multi: No telomeres count" "No telomeres:	1" "$OUT"

# Short contig: no crash, classified as none
OUT=$(build/bin/teloscope -f testFiles/short_contig.fa 2>/dev/null)
check_output_contains "short contig: no crash, type=none" "none" "$OUT"

# =====================================================
# GFA annotation tests
# =====================================================

rm -rf "$TMPDIR"/* 2>/dev/null || true
build/bin/teloscope -f testFiles/gfa_telo.gfa -o "$TMPDIR" 2>/dev/null
GFA_OUT="$TMPDIR/gfa_telo.gfa.telo.annotated.gfa"

# Check annotated GFA was produced
if [ -f "$GFA_OUT" ]; then
    echo -e "$PASS GFA annotation: output file created"
else
    echo -e "$FAIL GFA annotation: output file missing"
    EXIT=1
fi

# Check telomere nodes exist for segments with telomeres
check_output_contains "GFA: seg_t2t start node" "telomere_seg_t2t_start" "$(cat "$GFA_OUT")"
check_output_contains "GFA: seg_t2t end node" "telomere_seg_t2t_end" "$(cat "$GFA_OUT")"
check_output_contains "GFA: seg_ponly start node" "telomere_seg_ponly_start" "$(cat "$GFA_OUT")"
check_output_contains "GFA: seg_qonly end node" "telomere_seg_qonly_end" "$(cat "$GFA_OUT")"

# Check no telomere nodes for seg_none
check_output_not_contains "GFA: no telomere for seg_none" "telomere_seg_none" "$(cat "$GFA_OUT")"

# Check edge orientations: start → L telo + seg +, end → L telo + seg -
GFA_CONTENT=$(cat "$GFA_OUT")
check_output_contains "GFA: start edge orient (+)" "telomere_seg_t2t_start	+	seg_t2t	+" "$GFA_CONTENT"
check_output_contains "GFA: end edge orient (-)" "telomere_seg_t2t_end	+	seg_t2t	-" "$GFA_CONTENT"
check_output_contains "GFA: ponly start edge" "telomere_seg_ponly_start	+	seg_ponly	+" "$GFA_CONTENT"
check_output_contains "GFA: qonly end edge" "telomere_seg_qonly_end	+	seg_qonly	-" "$GFA_CONTENT"

# Check original edges are preserved
check_output_contains "GFA: original edge preserved" "seg_t2t	+	seg_ponly	+	0M" "$GFA_CONTENT"

# Check telomere node tags (LN, RC, TL)
check_output_contains "GFA: telomere LN tag" "LN:i:6" "$GFA_CONTENT"
check_output_contains "GFA: telomere RC tag" "RC:i:6000" "$GFA_CONTENT"
check_output_contains "GFA: telomere TL tag" "TL:i:600" "$GFA_CONTENT"

# =====================================================
# Error handling
# =====================================================

# Invalid -k 0
if build/bin/teloscope -f testFiles/t2t.fa -k 0 2>/dev/null; then
    echo -e "$FAIL -k 0 should fail"
    EXIT=1
else
    echo -e "$PASS -k 0 rejected"
fi

# Invalid -d 0
if build/bin/teloscope -f testFiles/t2t.fa -d 0 2>/dev/null; then
    echo -e "$FAIL -d 0 should fail"
    EXIT=1
else
    echo -e "$PASS -d 0 rejected"
fi

rm -rf "$TMPDIR"/* 2>/dev/null || true
exit $EXIT
