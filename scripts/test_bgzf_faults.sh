#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT/build/bin"}
OUTPUT=$(mktemp "${TMPDIR:-/tmp}/teloscope_bgzf_faults.XXXXXX")
trap 'rm -f "$OUTPUT"' EXIT
read -r -a COMPILE_FLAGS <<< "${FAULT_CXXFLAGS:--O0 -g}"
read -r -a LINK_FLAGS <<< "${FAULT_LDFLAGS:-}"

"${CXX:-g++}" -std=gnu++17 "${COMPILE_FLAGS[@]}" -I"$ROOT/include" \
    "$ROOT/tests/test_bgzf_faults.cpp" "$BUILD_DIR/.o/bgzf" \
    -Wl,--wrap=inflateInit2_ -Wl,--wrap=inflateEnd \
    -Wl,--wrap=deflateInit2_ -Wl,--wrap=deflate \
    -Wl,--wrap=deflateEnd -Wl,--wrap=deflateBound \
    "${LINK_FLAGS[@]}" -lz -o "$OUTPUT"
"$OUTPUT"
