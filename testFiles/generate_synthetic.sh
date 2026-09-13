#!/bin/bash
# Generate synthetic FASTA/GFA fixtures and testFiles/synthetic/manifest.tsv from declared intent.
# Byte-identical regeneration of existing fixtures is a hard constraint, enforced by --check.
#
# Usage:
#   generate_synthetic.sh [-o|--output OUTDIR]   write fixtures and the manifest
#   generate_synthetic.sh --check   regenerate into a temp dir and diff; nonzero on drift
#   generate_synthetic.sh --list   print fixture ids
#   generate_synthetic.sh --owned   print every path this script owns

set -euo pipefail
# Pathname expansion off: an S:<literal> escape hatch could otherwise glob on '*'.
set -f

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

readonly CANON_FWD="CCCTAA"
readonly CANON_REV="TTAGGG"
readonly PLANT_FWD="CCCTAAA"
readonly PLANT_REV="TTTAGGG"

readonly FILLER="ACGATCGATCGACTGACTGACGATCGATCGACTGACTGACGATCGATCGACTGACTGACGATCGATCGACTGACTGACGATCGATCGACTGACTG"

die() { printf 'generate_synthetic.sh: %s\n' "$*" >&2; exit 1; }

# Fixture expectations are stated against these defaults; fail loudly if one drifts.
assert_default() {
    local name=$1 expected=$2 file=$3 pattern=$4 actual
    [ -f "$REPO_ROOT/$file" ] || die "drift guard: $file not found"
    actual=$(sed -n "s/.*${pattern}.*/\1/p" "$REPO_ROOT/$file" | head -1)
    [ -n "$actual" ] || die "drift guard: could not read $name from $file"
    [ "$actual" = "$expected" ] || die \
        "drift guard: $name is $actual in $file but fixtures are written for $expected.
    Revisit the fixture expectations before changing this constant."
}

check_defaults() {
    assert_default maxMatchDist   50    include/input.h     'maxMatchDist *= *\([0-9]*\)'
    assert_default maxBlockDist   500   include/input.h     'maxBlockDist *= *\([0-9]*\)'
    assert_default minBlockLen    300   include/input.h     'minBlockLen *= *\([0-9]*\)'
    assert_default minBlockCounts 2     include/input.h     'minBlockCounts *= *\([0-9]*\)'
    assert_default terminalLimit  50000 include/input.h     'terminalLimit *= *\([0-9]*\)'
    assert_default termTolerance  3000  include/input.h     'terminalTolerance *= *\([0-9]*\)'
    assert_default linkDistance   1000  include/input.h     'linkDistance *= *\([0-9]*\)'
    assert_default editDistance   1     include/input.h     'editDistance *= *\([0-9]*\)'
    # hardcoded, not reachable from the CLI -- see docs/conflicts.md REG-003
    assert_default minCanonicalCount 4  src/teloscope.cpp   'minCanonicalCount *= *\([0-9]*\)'
    assert_default labelThreshold 0.667f include/input.h    'labelThreshold *= *\([0-9.]*f\)'
}

repeat_motif() {
    local unit=$1 n=$2 out="" cur=$1
    [ "$n" -ge 0 ] || die "repeat_motif: negative count"
    while [ "$n" -gt 0 ]; do
        if [ $((n & 1)) -eq 1 ]; then out="${out}${cur}"; fi
        cur="${cur}${cur}"
        n=$((n >> 1))
    done
    printf '%s' "$out"
}

make_filler() {
    local len=$1 out=""
    [ "$len" -ge 0 ] || die "make_filler: negative length"
    while [ ${#out} -lt "$len" ]; do out="${out}${FILLER}"; done
    printf '%s' "${out:0:$len}"
}

make_gap() {
    local len=$1 ch=${2:-N}
    repeat_motif "$ch" "$len"
}

hamming() {
    local a=$1 b=$2 i d=0
    [ ${#a} -eq ${#b} ] || { printf '99'; return; }
    for ((i = 0; i < ${#a}; i++)); do
        [ "${a:i:1}" = "${b:i:1}" ] || d=$((d + 1))
    done
    printf '%s' "$d"
}

# Assert a variant motif is 1-2 substitutions from some canonical motif.
assert_variant() {
    local motif=$1 best=99 d
    for canon in "$CANON_FWD" "$CANON_REV" "$PLANT_FWD" "$PLANT_REV"; do
        d=$(hamming "$motif" "$canon")
        [ "$d" -lt "$best" ] && best=$d
    done
    [ "$best" -ge 1 ] && [ "$best" -le 2 ] || die \
        "V:$motif is Hamming $best from every canonical motif; expected 1 or 2"
}

# DSL: F/R/V/M (repeats), L/l (filler), N/X (gaps), S (literal); '+' joins tokens.
build_tokens() {
    local spec=$1 out="" tok kind body unit count
    local IFS='+'
    for tok in $spec; do
        kind=${tok%%:*}
        body=${tok#*:}
        case "$kind" in
            F|R|V|M)
                unit=${body%x*}
                count=${body##*x}
                [ "$unit" != "$body" ] || die "token $tok: expected <unit>x<count>"
                case "$kind" in
                    F) [ "$unit" = "$CANON_FWD" ] || [ "$unit" = "$PLANT_FWD" ] || die \
                           "F:$unit is not a forward canonical motif" ;;
                    R) [ "$unit" = "$CANON_REV" ] || [ "$unit" = "$PLANT_REV" ] || die \
                           "R:$unit is not a reverse canonical motif" ;;
                    V) assert_variant "$unit" ;;
                esac
                out="${out}$(repeat_motif "$unit" "$count")"
                ;;
            L) out="${out}$(make_filler "$body")" ;;
            l) out="${out}$(make_filler "$body" | tr 'ACGT' 'acgt')" ;;
            N) out="${out}$(make_gap "$body" N)" ;;
            X) out="${out}$(make_gap "$body" X)" ;;
            S) out="${out}${body}" ;;
            *) die "unknown token kind '$kind' in '$tok'" ;;
        esac
    done
    printf '%s' "$out"
}

FX_ID=(); FX_PATH=(); FX_SPEC=(); FX_FLAGS=(); FX_EXPECT=(); FX_INTENT=()

fx() { # id path record_spec flags expect intent
    FX_ID+=("$1"); FX_PATH+=("$2"); FX_SPEC+=("$3")
    FX_FLAGS+=("$4"); FX_EXPECT+=("$5"); FX_INTENT+=("$6")
}

GFX_ID=(); GFX_PATH=(); GFX_BODY=(); GFX_INTENT=()

gfx() { # id path body intent
    GFX_ID+=("$1"); GFX_PATH+=("$2"); GFX_BODY+=("$3"); GFX_INTENT+=("$4")
}

# A file this script does not own: not written, not checked, only carried in the manifest.
XFX_ID=(); XFX_PATH=(); XFX_SCAFFOLD=(); XFX_FLAGS=(); XFX_EXPECT=(); XFX_INTENT=()

xfx() { # id path scaffold flags expect intent
    XFX_ID+=("$1"); XFX_PATH+=("$2"); XFX_SCAFFOLD+=("$3")
    XFX_FLAGS+=("$4"); XFX_EXPECT+=("$5"); XFX_INTENT+=("$6")
}

write_fasta() {
    local dest=$1 spec=$2 rec header body
    local IFS=';'
    : > "$dest"
    for rec in $spec; do
        header=${rec%%=*}
        body=${rec#*=}
        [ "$header" != "$rec" ] || die "record '$rec' is missing '<header>='"
        {
            printf '>%s\n' "$header"
            build_tokens "$body"
            printf '\n'
        } >> "$dest"
    done
}

source "$SCRIPT_DIR/synthetic_fixtures.sh"

MANIFEST_COLUMNS="id	path	scaffold	flags	expect_type	expect_anomaly	expect_granular	expect_telomeres	expect_labels	expect_gaps	expect_its	expect_telolen	intent"

# Absent key = default; key present but empty means the field really is empty.
field() { # expect_string key default
    case ";$1;" in
        *";$2="*) printf '%s' "$(printf '%s' "$1" | tr ';' '\n' | sed -n "s/^$2=//p" | head -1)" ;;
        *)        printf '%s' "$3" ;;
    esac
}

write_manifest() {
    local dir=$1 i rec header spec out="$1/synthetic/manifest.tsv"
    mkdir -p "$dir/synthetic"
    {
        printf '# Expected classification per fixture, derived from docs/classification.md and the\n'
        printf '# fixture construction -- never read back from the binary. Regenerate with\n'
        printf '# testFiles/generate_synthetic.sh; CI checks that regeneration is a no-op.\n'
        printf '# A "?" in any expect_ column means the outcome is blocked on a docs/conflicts.md\n'
        printf '# decision; scripts/test_synthetic_intent.py fails on it rather than passing silently.\n'
        printf '%s\n' "$MANIFEST_COLUMNS"
        local -a records expects
        for i in "${!FX_ID[@]}"; do
            IFS=';' read -r -a records <<< "${FX_SPEC[$i]}"
            IFS='|'  read -r -a expects <<< "${FX_EXPECT[$i]}"
            [ "${#expects[@]}" -eq 1 ] || [ "${#expects[@]}" -eq "${#records[@]}" ] || die \
                "${FX_ID[$i]}: ${#records[@]} records but ${#expects[@]} expectations"
            local n=0 e
            for rec in "${records[@]}"; do
                header=${rec%%=*}
                if [ "${#expects[@]}" -eq 1 ]; then e=${expects[0]}; else e=${expects[$n]}; fi
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "${FX_ID[$i]}" "${FX_PATH[$i]}" "$header" "${FX_FLAGS[$i]}" \
                    "$(field "$e" type '?')" \
                    "$(field "$e" anom '?')" \
                    "$(field "$e" gran '?')" \
                    "$(field "$e" telo '?')" \
                    "$(field "$e" labels '?')" \
                    "$(field "$e" gaps '?')" \
                    "$(field "$e" its '-')" \
                    "$(field "$e" telolen '-')" \
                    "${FX_INTENT[$i]}"
                n=$((n + 1))
            done
        done
        local j e
        for j in "${!XFX_ID[@]}"; do
            e=${XFX_EXPECT[$j]}
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${XFX_ID[$j]}" "${XFX_PATH[$j]}" "${XFX_SCAFFOLD[$j]}" "${XFX_FLAGS[$j]}" \
                "$(field "$e" type '?')" \
                "$(field "$e" anom '?')" \
                "$(field "$e" gran '?')" \
                "$(field "$e" telo '?')" \
                "$(field "$e" labels '?')" \
                "$(field "$e" gaps '?')" \
                "$(field "$e" its '-')" \
                "$(field "$e" telolen '-')" \
                "${XFX_INTENT[$j]}"
        done
    } > "$out"
}

emit_all() {
    local dir=$1 i
    for i in "${!FX_ID[@]}"; do
        mkdir -p "$(dirname "$dir/${FX_PATH[$i]}")"
        write_fasta "$dir/${FX_PATH[$i]}" "${FX_SPEC[$i]}"
    done
    for i in "${!GFX_ID[@]}"; do
        mkdir -p "$(dirname "$dir/${GFX_PATH[$i]}")"
        printf '%s' "${GFX_BODY[$i]}" > "$dir/${GFX_PATH[$i]}"
    done
    write_manifest "$dir"
}

owned_paths() {
    local i
    {
        for i in "${!FX_ID[@]}"; do printf '%s\n' "${FX_PATH[$i]}"; done
        for i in "${!GFX_ID[@]}"; do printf '%s\n' "${GFX_PATH[$i]}"; done
        printf 'synthetic/manifest.tsv\n'
    } | sort -u
}

self_test() {
    local motif
    for motif in "$CANON_FWD" "$CANON_REV" "$PLANT_FWD" "$PLANT_REV"; do
        case "$FILLER" in
            *"$motif"*) die "FILLER contains the telomeric motif $motif" ;;
        esac
    done
    [ "${#FILLER}" -eq 95 ] || die "FILLER length changed; every filler offset depends on it"
}

OUTDIR="$SCRIPT_DIR"
MODE=write
while [ $# -gt 0 ]; do
    case "$1" in
        -o|--output) OUTDIR=$2; shift 2 ;;
        --check)     MODE=check; shift ;;
        --list)      MODE=list; shift ;;
        --owned)     MODE=owned; shift ;;
        -h|--help)   sed -n '2,9p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)           die "unknown argument '$1'" ;;
    esac
done

case "$MODE" in
    list)  for i in "${!FX_ID[@]}"; do printf '%s\n' "${FX_ID[$i]}"; done
           for i in "${!GFX_ID[@]}"; do printf '%s\n' "${GFX_ID[$i]}"; done
           for i in "${!XFX_ID[@]}"; do printf '%s\n' "${XFX_ID[$i]}"; done ;;
    owned) owned_paths ;;
    check)
        self_test; check_defaults
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        emit_all "$tmp"
        drift=0
        while read -r rel; do
            if ! cmp -s "$tmp/$rel" "$OUTDIR/$rel" 2>/dev/null; then
                printf 'DRIFT %s\n' "$rel"; drift=$((drift + 1))
            fi
        done < <(owned_paths)
        total=$(owned_paths | wc -l)
        if [ "$drift" -eq 0 ]; then
            printf 'checked %s owned files against %s: 0 differ\n' "$total" "$OUTDIR"
            printf 'defaults verified against source: OK\n'
        else
            printf '%s of %s owned files differ from %s\n' "$drift" "$total" "$OUTDIR" >&2
            exit 1
        fi
        ;;
    write)
        self_test; check_defaults
        emit_all "$OUTDIR"
        printf 'Generated %s fixtures and 1 manifest in %s/\n' \
            "$(( ${#FX_ID[@]} + ${#GFX_ID[@]} ))" "$OUTDIR"
        ;;
esac
