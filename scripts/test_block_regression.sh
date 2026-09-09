#!/usr/bin/env bash
# Compare block calls between the current tree and a reference commit.
# Builds the reference in a throwaway worktree so the current build is untouched.
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

OLD_REF=${OLD_REF:-$(git merge-base origin/main HEAD 2>/dev/null || echo main)}
ROSTER=${ROSTER:-validateFiles/repair_roster.tsv}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/teloscope_regress.XXXXXX")
trap 'rm -rf "$WORK"; git worktree remove --force "$WORK/old" 2>/dev/null || true' EXIT

echo "reference: $OLD_REF"

echo "== building current =="
make head -j"$(nproc 2>/dev/null || echo 4)" >/dev/null

echo "== building reference in a worktree =="
git worktree add -q --detach "$WORK/old" "$OLD_REF"
git -C "$WORK/old" submodule update --init --recursive -q
make -C "$WORK/old" head -j"$(nproc 2>/dev/null || echo 4)" >/dev/null

echo "== extracting =="
python3 scripts/check_block_regression.py extract \
    --binary "$WORK/old/build/bin/teloscope" --out "$WORK/old_rows.tsv"
python3 scripts/check_block_regression.py extract \
    --binary "$ROOT/build/bin/teloscope" --out "$WORK/new_rows.tsv"

cp "$WORK/old_rows.tsv" "$WORK/new_rows.tsv" "$ROOT/" 2>/dev/null || true

echo "== comparing =="
python3 scripts/check_block_regression.py compare \
    --old "$WORK/old_rows.tsv" --new "$WORK/new_rows.tsv" --roster "$ROSTER"
