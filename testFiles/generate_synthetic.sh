#!/bin/bash
# Generate synthetic FASTA test files with known telomere placements.
# Each file is small and deterministic for CI.

DIR="$(dirname "$0")"

# Helper: repeat a motif N times
repeat_motif() { printf "%0.s$1" $(seq 1 "$2"); }

# Filler that contains no telomeric 6-mers (TTAGGG, CCCTAA, or common variants)
FILLER="ACGATCGATCGACTGACTGACGATCGATCGACTGACTGACGATCGATCGACTGACTGACGATCGATCGACTGACTGACGATCGATCGACTGACTG"

# Repeat filler to desired length (approximate)
make_filler() {
    local len=$1
    local out=""
    while [ ${#out} -lt "$len" ]; do
        out="${out}${FILLER}"
    done
    echo "${out:0:$len}"
}

# ============================================================
# 1. t2t.fa — Both arms, correct orientation → t2t
# ============================================================
{
    echo ">chr_t2t"
    parm=$(repeat_motif "CCCTAA" 100)   # 600bp p-arm
    qarm=$(repeat_motif "TTAGGG" 100)   # 600bp q-arm
    mid=$(make_filler 2000)
    echo "${parm}${mid}${qarm}"
} > "$DIR/t2t.fa"

# ============================================================
# 2. incomplete_p.fa — Only p-arm telomere → incomplete
# ============================================================
{
    echo ">chr_incomplete_p"
    parm=$(repeat_motif "CCCTAA" 100)
    tail=$(make_filler 2400)
    echo "${parm}${tail}"
} > "$DIR/incomplete_p.fa"

# ============================================================
# 3. incomplete_q.fa — Only q-arm telomere → incomplete
# ============================================================
{
    echo ">chr_incomplete_q"
    head=$(make_filler 2400)
    qarm=$(repeat_motif "TTAGGG" 100)
    echo "${head}${qarm}"
} > "$DIR/incomplete_q.fa"

# ============================================================
# 4. no_telo.fa — No telomeres → none
# ============================================================
{
    echo ">chr_none"
    make_filler 3000
} > "$DIR/no_telo.fa"

# ============================================================
# 5. misassembly.fa — Two p-arms (Pp) → misassembly
# Two CCCTAA blocks both in valid positions (left half).
# The Pp check detects a second p-block → misassembly.
# ============================================================
{
    echo ">chr_misassembly"
    parm1=$(repeat_motif "CCCTAA" 100)  # 600bp p-arm at start
    mid1=$(make_filler 500)
    parm2=$(repeat_motif "CCCTAA" 100)  # 600bp second p-arm
    tail=$(make_filler 5000)
    echo "${parm1}${mid1}${parm2}${tail}"
} > "$DIR/misassembly.fa"

# ============================================================
# 6. discordant.fa — p-arm at q-side position → discordant
# ============================================================
{
    echo ">chr_discordant"
    head=$(make_filler 2800)
    parm=$(repeat_motif "CCCTAA" 100)
    echo "${head}${parm}"
} > "$DIR/discordant.fa"

# ============================================================
# 7. gapped_t2t.fa — T2T with internal N-gap → gapped_t2t
# ============================================================
{
    echo ">chr_gapped_t2t"
    parm=$(repeat_motif "CCCTAA" 100)
    mid1=$(make_filler 900)
    gap=$(printf 'N%.0s' $(seq 1 100))
    mid2=$(make_filler 900)
    qarm=$(repeat_motif "TTAGGG" 100)
    echo "${parm}${mid1}${gap}${mid2}${qarm}"
} > "$DIR/gapped_t2t.fa"

# ============================================================
# 8. multi.fa — Multiple contigs: t2t + none + incomplete
# ============================================================
{
    echo ">contig_t2t"
    parm=$(repeat_motif "CCCTAA" 100)
    qarm=$(repeat_motif "TTAGGG" 100)
    mid=$(make_filler 2000)
    echo "${parm}${mid}${qarm}"

    echo ">contig_none"
    make_filler 2000

    echo ">contig_incomplete"
    qarm=$(repeat_motif "TTAGGG" 100)
    head=$(make_filler 2400)
    echo "${head}${qarm}"
} > "$DIR/multi.fa"

# ============================================================
# 9. plant.fa — 7-mer canonical (CCCTAAA/TTTAGGG)
# ============================================================
{
    echo ">chr_plant"
    parm=$(repeat_motif "CCCTAAA" 86)  # 602bp
    qarm=$(repeat_motif "TTTAGGG" 86)  # 602bp
    mid=$(make_filler 2000)
    echo "${parm}${mid}${qarm}"
} > "$DIR/plant.fa"

# ============================================================
# 10. edit_test.fa — Mutated repeats (1 sub per repeat)
# Each repeat has exactly 1 substitution: TTAGAG instead of TTAGGG
# Detectable at -x 1 but not -x 0
# ============================================================
{
    echo ">chr_edit_p"
    parm=$(repeat_motif "CTCTAA" 100)   # 1-sub variant of CCCTAA
    mid=$(make_filler 2000)
    qarm=$(repeat_motif "TTAGGG" 100)   # exact canonical
    echo "${parm}${mid}${qarm}"
} > "$DIR/edit_test.fa"

# ============================================================
# 11. short_contig.fa — Very short contig (shorter than one window)
# ============================================================
{
    echo ">chr_tiny"
    make_filler 100
} > "$DIR/short_contig.fa"

# ============================================================
# 12. balanced.fa — Interleaved fwd/rev in one region → 'b' label
# Alternating CCCTAA and TTAGGG repeats create a balanced block
# ============================================================
{
    echo ">chr_balanced"
    head=$(make_filler 1200)
    # 100 interleaved pairs = 1200bp
    bal=""
    for i in $(seq 1 100); do
        bal="${bal}CCCTAATTAGGG"
    done
    tail=$(make_filler 1200)
    echo "${head}${bal}${tail}"
} > "$DIR/balanced.fa"

# ============================================================
# 13. its.fa — Interstitial telomere blocks (middle of sequence)
# Telomere at both ends + telomeric block in the middle
# The middle block must be far enough from ends to not be terminal
# ============================================================
{
    echo ">chr_its"
    parm=$(repeat_motif "CCCTAA" 100)    # 600bp p-arm
    mid1=$(make_filler 3000)
    its=$(repeat_motif "TTAGGG" 100)     # 600bp interstitial
    mid2=$(make_filler 3000)
    qarm=$(repeat_motif "TTAGGG" 100)    # 600bp q-arm
    echo "${parm}${mid1}${its}${mid2}${qarm}"
} > "$DIR/its.fa"

# ============================================================
# 14. density_edge.fa — Block right at density threshold
# Mix canonical repeats with filler to get ~50% density
# ============================================================
{
    echo ">chr_density"
    # Alternating: 6bp match + 6bp filler = 50% density
    parm=""
    for i in $(seq 1 100); do
        parm="${parm}CCCTAAACGATC"
    done
    mid=$(make_filler 2000)
    qarm=$(repeat_motif "TTAGGG" 100)
    echo "${parm}${mid}${qarm}"
} > "$DIR/density_edge.fa"

# ============================================================
# 15. misassembly_qq.fa — Two q-arms (Qq) → misassembly
# Two TTAGGG blocks both in valid positions (right half).
# ============================================================
{
    echo ">chr_misassembly_qq"
    head=$(make_filler 5000)
    qarm1=$(repeat_motif "TTAGGG" 100)  # 600bp first q-arm
    mid=$(make_filler 500)
    qarm2=$(repeat_motif "TTAGGG" 100)  # 600bp second q-arm at end
    echo "${head}${mid}${qarm1}${mid}${qarm2}"
} > "$DIR/misassembly_qq.fa"

# ============================================================
# 16. gapped_misassembly.fa — Pp with internal N-gap → gapped_misassembly
# Gap placed far from both p-blocks to avoid segment boundary issues
# ============================================================
{
    echo ">chr_gapped_misassembly"
    parm1=$(repeat_motif "CCCTAA" 100)  # 600bp p-arm at start
    mid1=$(make_filler 500)
    parm2=$(repeat_motif "CCCTAA" 100)  # 600bp second p-arm
    mid2=$(make_filler 2000)
    gap=$(printf 'N%.0s' $(seq 1 100))  # 100bp N-gap in tail region
    tail=$(make_filler 3000)
    echo "${parm1}${mid1}${parm2}${mid2}${gap}${tail}"
} > "$DIR/gapped_misassembly.fa"

# ============================================================
# 17. gapped_incomplete.fa — Single p-arm with N-gap
# ============================================================
{
    echo ">chr_gapped_incomplete"
    parm=$(repeat_motif "CCCTAA" 100)
    mid=$(make_filler 1000)
    gap=$(printf 'N%.0s' $(seq 1 100))
    tail=$(make_filler 2000)
    echo "${parm}${mid}${gap}${tail}"
} > "$DIR/gapped_incomplete.fa"

# ============================================================
# 18. gapped_none.fa — No telomeres but has N-gap
# ============================================================
{
    echo ">chr_gapped_none"
    head=$(make_filler 1000)
    gap=$(printf 'N%.0s' $(seq 1 100))
    tail=$(make_filler 2000)
    echo "${head}${gap}${tail}"
} > "$DIR/gapped_none.fa"

# ============================================================
# 19. gapped_discordant.fa — Discordant p at q-side with N-gap
# ============================================================
{
    echo ">chr_gapped_discordant"
    head=$(make_filler 2000)
    gap=$(printf 'N%.0s' $(seq 1 100))
    mid=$(make_filler 500)
    parm=$(repeat_motif "CCCTAA" 100)
    echo "${head}${gap}${mid}${parm}"
} > "$DIR/gapped_discordant.fa"

# ============================================================
# 20. gfa_telo.gfa — GFA with segments having telomeric tips
# Four segments: t2t, p-only, q-only, no-telo
# Plus edges between them (simulating a linear assembly graph)
# ============================================================
{
    parm=$(repeat_motif "CCCTAA" 100)   # 600bp p-arm
    qarm=$(repeat_motif "TTAGGG" 100)   # 600bp q-arm
    filler=$(make_filler 2000)

    echo "H	VN:Z:1.2"
    echo "S	seg_t2t	${parm}${filler}${qarm}"
    echo "S	seg_ponly	${parm}${filler}"
    echo "S	seg_qonly	${filler}${qarm}"
    echo "S	seg_none	${filler}"
    echo "L	seg_t2t	+	seg_ponly	+	0M"
    echo "L	seg_ponly	+	seg_qonly	+	0M"
    echo "L	seg_qonly	+	seg_none	+	0M"
} > "$DIR/gfa_telo.gfa"

echo "Generated $(ls -1 "$DIR"/*.fa "$DIR"/*.gfa 2>/dev/null | wc -l) synthetic test files in $DIR/"
