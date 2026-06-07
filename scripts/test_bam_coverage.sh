#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD=$(mktemp -d "${TMPDIR:-/tmp}/teloscope_bam_coverage.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT
cp -a "$ROOT/gfalibs" "$BUILD/gfalibs"
find "$BUILD/gfalibs" -maxdepth 1 -type f \
    \( -name '*.o' -o -name '*.gcda' -o -name '*.gcno' \) -delete

make -C "$ROOT" BUILD="$BUILD/bin" GFALIBS_DIR="$BUILD/gfalibs" \
    CFLAGS="-O0 --coverage" LDFLAGS="-pthread --coverage" head -j2

TELOSCOPE="$BUILD/bin/teloscope" BAM_MUTATION_CASES="${BAM_MUTATION_CASES:-512}" \
    python3 "$ROOT/scripts/test_bam_subset.py"
BUILD_DIR="$BUILD/bin" FAULT_CXXFLAGS="-O0 -g --coverage" \
    FAULT_LDFLAGS="--coverage" bash "$ROOT/scripts/test_bgzf_faults.sh"
python3 "$ROOT/scripts/check_bam_coverage.py" --object-dir "$BUILD/bin/.o"
