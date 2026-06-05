#!/usr/bin/env bash
# Diff teloscope telomere caps against the reference annotator on shared >10 kb inputs.
# Reference script is not vendored; set REF=/path if it moved.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/.." && pwd)"
telo="$repo/build/bin/teloscope"
ref="${REF:-$repo/gfa/remove_nodes_add_telomere.py}"

[ -x "$telo" ] || { echo "build teloscope first (make)"; exit 2; }
[ -f "$ref" ]  || { echo "missing reference script: $ref (set REF=/path/to/remove_nodes_add_telomere.py)"; exit 2; }

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# Generate the two scenarios: pathless '+' carriers (p/q/t2t) and a path-aware '-' carrier.
python3 - "$work" <<'PY'
import sys, os
w = sys.argv[1]
P = "CCCTAA"*120          # 720 bp start motif
Q = "TTAGGG"*120          # 720 bp end motif
MID = "ACGT"*2500         # 10 kb filler -> contigs ~11.4 kb (> 10 kb threshold)

def seg(name, s): return "S\t%s\t%s\n" % (name, s)

# Scenario 1: no paths -> teloscope path-less; reference uses the single-node shortcut.
with open(os.path.join(w,"s1.gfa"),"w") as f:
    f.write("H\tVN:Z:1.2\n")
    f.write(seg("utig_p", P+MID))        # start only
    f.write(seg("utig_q", MID+Q))        # end only
    f.write(seg("utig_t", P+MID+Q))      # both
with open(os.path.join(w,"s1.scfmap"),"w") as f:
    for c in ("p","q","t"): f.write("path contig_%s utig_%s\n" % (c,c))
with open(os.path.join(w,"s1.bed"),"w") as f:
    f.write("contig_p\t0\t720\t720\tp\n")
    f.write("contig_q\t10720\t11440\t720\tq\n")
    f.write("contig_t\t0\t720\t720\tp\n")
    f.write("contig_t\t10720\t11440\t720\tq\n")
open(os.path.join(w,"empty.paths"),"w").close()

# Scenario 2: a unitig carried reverse ('-'); telomere at its physical END = contig START.
with open(os.path.join(w,"s2.gfa"),"w") as f:
    f.write("H\tVN:Z:1.2\n")
    f.write(seg("utig_r", MID+Q))                 # telomere at physical end
    f.write("P\tcontig_r\tutig_r-\t*\n")          # carried reverse
with open(os.path.join(w,"s2.scfmap"),"w") as f:
    f.write("path contig_r pathR\n")              # not 'utig...' -> use paths.tsv
with open(os.path.join(w,"s2.paths"),"w") as f:
    f.write("pathR\tutig_r-\n")
with open(os.path.join(w,"s2.bed"),"w") as f:
    f.write("contig_r\t0\t720\t720\tp\n")         # telomere at contig start
PY

# Telomere S/L lines minus tool-specific tags; real tab since grep/sed escapes vary.
TAB=$'\t'
norm() { grep -E "^[SL]${TAB}telomere_" "$1" 2>/dev/null | sed -E "s/${TAB}(TL|RC):i:[0-9]+//g" | sort; }

status=0
run_case() {
    local tag="$1" gfa="$2" scfmap="$3" bed="$4" paths="$5"
    python3 "$ref" --telo "$bed" --gfa "$gfa" --scfmap "$scfmap" --paths "$paths" \
        --output "$work/${tag}.ref.gfa" --colors "$work/${tag}.ref.csv" 2>/dev/null
    "$telo" -f "$gfa" -o "$work" >/dev/null 2>&1
    local mine="$work/$(basename "$gfa").telo.annotated.gfa"
    local caps; caps=$(norm "$mine" | grep -c '^L' || true)
    if diff <(norm "$work/${tag}.ref.gfa") <(norm "$mine") >"$work/${tag}.diff"; then
        echo "PASS  $tag  ($caps caps match the reference)"
    else
        echo "FAIL  $tag  (< reference, > teloscope)"; cat "$work/${tag}.diff"; status=1
    fi
}

run_case pathless "$work/s1.gfa" "$work/s1.scfmap" "$work/s1.bed" "$work/empty.paths"
run_case pathrev  "$work/s2.gfa" "$work/s2.scfmap" "$work/s2.bed" "$work/s2.paths"

exit $status
