#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD=$(mktemp -d "${TMPDIR:-/tmp}/teloscope_bam_sanitize.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT
cp -a "$ROOT/gfalibs" "$BUILD/gfalibs"
find "$BUILD/gfalibs" -maxdepth 1 -type f \
    \( -name '*.o' -o -name '*.gcda' -o -name '*.gcno' \) -delete

SANITIZERS="-O1 -fno-omit-frame-pointer -fsanitize=address,undefined"
make -C "$ROOT" BUILD="$BUILD/bin" GFALIBS_DIR="$BUILD/gfalibs" CFLAGS="$SANITIZERS" \
    LDFLAGS="-pthread -fsanitize=address,undefined" head -j2

ASAN_OPTIONS="${ASAN_OPTIONS:-abort_on_error=1:detect_leaks=0}" \
UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}" \
TELOSCOPE="$BUILD/bin/teloscope" BAM_MUTATION_CASES="${BAM_MUTATION_CASES:-512}" \
    python3 "$ROOT/scripts/test_bam_subset.py"

ASAN_OPTIONS="${ASAN_OPTIONS:-abort_on_error=1:detect_leaks=0}" \
UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}" \
BUILD_DIR="$BUILD/bin" FAULT_CXXFLAGS="$SANITIZERS" \
    FAULT_LDFLAGS="-fsanitize=address,undefined" \
    bash "$ROOT/scripts/test_bgzf_faults.sh"
