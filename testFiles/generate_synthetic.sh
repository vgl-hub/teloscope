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

# ============================================================
# 21. mirror_inverted_short.fa — Swapped arms on short contig
# TTAGGG (rev/q-pattern) at start + CCCTAA (fwd/p-pattern) at end.
# Short contig: entire sequence within terminalLimit → both found.
# Both at wrong end → both hasValidOr=false → discordant.
# ============================================================
{
    echo ">chr_mirror_inv_short"
    rev_at_start=$(repeat_motif "TTAGGG" 100)   # 600bp q-pattern at p-side
    mid=$(make_filler 2000)
    fwd_at_end=$(repeat_motif "CCCTAA" 100)     # 600bp p-pattern at q-side
    echo "${rev_at_start}${mid}${fwd_at_end}"
} > "$DIR/mirror_inverted_short.fa"

# ============================================================
# 22. mirror_rev_start.fa — q-pattern at p-side (mirror of discordant.fa)
# discordant.fa places CCCTAA (p) at q-side; this places
# TTAGGG (q) at p-side. Verifies symmetric discordant detection.
# ============================================================
{
    echo ">chr_mirror_rev_start"
    qarm=$(repeat_motif "TTAGGG" 100)   # 600bp q-pattern at start
    tail=$(make_filler 2400)
    echo "${qarm}${tail}"
} > "$DIR/mirror_rev_start.fa"

# ============================================================
# 23. mirror_inverted_long.fa — Swapped arms on long contig
# Same as mirror_inverted_short but contig exceeds 2*terminalLimit.
# Directional scan: fwd scans from start (finds nothing),
# rev scans from end (finds nothing) → none.
# Uses -t 1000 to keep file small.
# ============================================================
{
    echo ">chr_mirror_inv_long"
    rev_at_start=$(repeat_motif "TTAGGG" 100)   # 600bp q-pattern at p-side
    mid=$(make_filler 4000)
    fwd_at_end=$(repeat_motif "CCCTAA" 100)     # 600bp p-pattern at q-side
    echo "${rev_at_start}${mid}${fwd_at_end}"
} > "$DIR/mirror_inverted_long.fa"

# ============================================================
# 24. mirror_fwd_end_long.fa — p-pattern only at q-side, long contig
# CCCTAA at far end, nothing at start. Directional scan from start
# finds no fwd matches in zone → no blocks → none.
# Uses -t 1000.
# ============================================================
{
    echo ">chr_mirror_fwd_end_long"
    head=$(make_filler 4600)
    fwd_at_end=$(repeat_motif "CCCTAA" 100)     # 600bp at end
    echo "${head}${fwd_at_end}"
} > "$DIR/mirror_fwd_end_long.fa"

# ============================================================
# 25. mirror_rev_start_long.fa — q-pattern only at p-side, long contig
# TTAGGG at start, nothing at end. Directional scan from end
# finds no rev matches in zone → no blocks → none.
# Uses -t 1000.
# ============================================================
{
    echo ">chr_mirror_rev_start_long"
    rev_at_start=$(repeat_motif "TTAGGG" 100)   # 600bp at start
    tail=$(make_filler 4600)
    echo "${rev_at_start}${tail}"
} > "$DIR/mirror_rev_start_long.fa"

# ============================================================
# 26. mirror_both_start.fa — Both strands at p-side
# CCCTAA (p) + TTAGGG (q) both at start. p-block valid at start,
# q-block invalid at start (hasValidOr=false) → discordant.
# ============================================================
{
    echo ">chr_mirror_both_start"
    parm=$(repeat_motif "CCCTAA" 100)    # 600bp fwd at start
    qarm=$(repeat_motif "TTAGGG" 100)    # 600bp rev right after
    tail=$(make_filler 2000)
    echo "${parm}${qarm}${tail}"
} > "$DIR/mirror_both_start.fa"

# ============================================================
# 27. mirror_both_end.fa — Both strands at q-side
# filler + CCCTAA (p) + TTAGGG (q) at end. q-block valid at end,
# p-block invalid at end (hasValidOr=false) → discordant.
# ============================================================
{
    echo ">chr_mirror_both_end"
    head=$(make_filler 2000)
    parm=$(repeat_motif "CCCTAA" 100)    # 600bp fwd near end
    qarm=$(repeat_motif "TTAGGG" 100)    # 600bp rev at end
    echo "${head}${parm}${qarm}"
} > "$DIR/mirror_both_end.fa"

# ============================================================
# 28. mirror_extend.fa — Terminal block extending past zone
# Dense CCCTAA from start, length > terminalLimit.
# Block starts in zone and extends past it.
# Uses -t 300.
# ============================================================
{
    echo ">chr_mirror_extend"
    parm=$(repeat_motif "CCCTAA" 100)    # 600bp p-arm (exceeds -t 300 zone)
    tail=$(make_filler 2400)
    echo "${parm}${tail}"
} > "$DIR/mirror_extend.fa"

# ============================================================
# 29. mirror_merge.fa — Sub-block merging via blockDist
# Two CCCTAA clusters separated by 100bp filler (> matchDist=50,
# < blockDist=200). Phase 1 creates 2 sub-blocks, Phase 2 merges.
# ============================================================
{
    echo ">chr_mirror_merge"
    clust1=$(repeat_motif "CCCTAA" 50)   # 300bp cluster
    gap=$(make_filler 100)               # 100bp gap (50 < 100 < 200)
    clust2=$(repeat_motif "CCCTAA" 50)   # 300bp cluster
    tail=$(make_filler 2000)
    echo "${clust1}${gap}${clust2}${tail}"
} > "$DIR/mirror_merge.fa"

# ============================================================
# 30. mirror_no_merge.fa — Sub-blocks too far apart to merge
# Two CCCTAA clusters separated by 300bp filler (> blockDist=200).
# Each sub-block has blockLen=240 < minBlockLen=500 → both fail.
# ============================================================
{
    echo ">chr_mirror_no_merge"
    clust1=$(repeat_motif "CCCTAA" 40)   # 240bp cluster
    gap=$(make_filler 300)               # 300bp gap (> 200)
    clust2=$(repeat_motif "CCCTAA" 40)   # 240bp cluster
    tail=$(make_filler 2000)
    echo "${clust1}${gap}${clust2}${tail}"
} > "$DIR/mirror_no_merge.fa"

# ============================================================
# 31. boundary_no_terminal.fa — Matches outside terminal zones
# With -t 500: no matches in [0,500) or [6700,7200)
# → no terminal blocks → fwdBoundary=0 revBoundary=7200
# → ITS covers everything → 2 ITS blocks
# ============================================================
{
    echo ">chr_boundary_no_term"
    head=$(make_filler 2000)
    fwd_mid=$(repeat_motif "CCCTAA" 100)    # 600bp at pos 2000
    gap=$(make_filler 2000)
    rev_mid=$(repeat_motif "TTAGGG" 100)    # 600bp at pos 4600
    tail=$(make_filler 2000)
    echo "${head}${fwd_mid}${gap}${rev_mid}${tail}"
} > "$DIR/boundary_no_terminal.fa"

# ============================================================
# 32. boundary_fail_filter.fa — Terminal block too short
# 20 CCCTAA at start = 120bp block. Passes Phase 1 (counts≥2,
# canonical>0) but fails Phase 2 (blockLen=120 < minBlockLen=500).
# Boundary stays at 0 → matches available for ITS.
# ============================================================
{
    echo ">chr_boundary_fail"
    short_telo=$(repeat_motif "CCCTAA" 20)  # 120bp at start
    tail=$(make_filler 3000)
    echo "${short_telo}${tail}"
} > "$DIR/boundary_fail_filter.fa"

# ============================================================
# 33. boundary_cross.fa — Boundaries meet, no ITS
# Dense adjacent p/q blocks: fwdBoundary=1200, revBoundary=1200.
# 1200 < 1200 is false → ITS condition fails → 0 ITS blocks.
# ============================================================
{
    echo ">chr_boundary_cross"
    parm=$(repeat_motif "CCCTAA" 200)       # 1200bp p-arm
    qarm=$(repeat_motif "TTAGGG" 200)       # 1200bp q-arm
    echo "${parm}${qarm}"
} > "$DIR/boundary_cross.fa"

# ============================================================
# 34. boundary_zone_shift.fa — Zone size controls terminal vs ITS
# Telomere blocks at offset 500 and offset 3100.
# -t 1200: both in zone → terminal (t2t)
# -t 400: both outside zone → ITS only (none)
# ============================================================
{
    echo ">chr_boundary_zone"
    head=$(make_filler 500)
    parm=$(repeat_motif "CCCTAA" 100)       # 600bp at pos 500
    mid=$(make_filler 2000)
    qarm=$(repeat_motif "TTAGGG" 100)       # 600bp at pos 3100
    tail=$(make_filler 500)
    echo "${head}${parm}${mid}${qarm}${tail}"
} > "$DIR/boundary_zone_shift.fa"

# ============================================================
# 35. boundary_its_at_edge.fa — ITS starts exactly at fwdBoundary
# p-block at 0-600 → fwdBoundary=600.
# Small TTAGGG cluster at 600-660 right at boundary → ITS picks it up.
# Large TTAGGG at end → q-block, revBoundary well before mid.
# ============================================================
{
    echo ">chr_boundary_edge"
    parm=$(repeat_motif "CCCTAA" 100)       # 600bp p-arm
    its_cluster=$(repeat_motif "TTAGGG" 10) # 60bp at pos 600
    mid=$(make_filler 2200)
    qarm=$(repeat_motif "TTAGGG" 100)       # 600bp q-arm
    echo "${parm}${its_cluster}${mid}${qarm}"
} > "$DIR/boundary_its_at_edge.fa"

# ============================================================
# 36. boundary_multiple_p.fa — Two separate p-blocks
# Gap=300 > blockDist=200 → no Phase 2 merge → 2 terminal blocks.
# fwdBoundary moves to end of outermost block (1500).
# ITS starts at 1500, nothing there → 0 ITS.
# ============================================================
{
    echo ">chr_boundary_multi_p"
    clust1=$(repeat_motif "CCCTAA" 100)     # 600bp at pos 0
    gap=$(make_filler 300)                  # 300bp gap > blockDist
    clust2=$(repeat_motif "CCCTAA" 100)     # 600bp at pos 900
    tail=$(make_filler 2000)
    echo "${clust1}${gap}${clust2}${tail}"
} > "$DIR/boundary_multiple_p.fa"

# ============================================================
# 37. boundary_extend_its.fa — Terminal extends past zone + ITS
# With -t 300: p-block starts in [0,300) zone, extends to 600.
# q-block starts in [4200,4500) zone, extends left to 3900.
# ITS region [600,3900) contains TTAGGG cluster at 2600.
# ============================================================
{
    echo ">chr_boundary_ext_its"
    parm=$(repeat_motif "CCCTAA" 100)       # 600bp p-arm
    mid1=$(make_filler 2000)
    its_telo=$(repeat_motif "TTAGGG" 50)    # 300bp ITS at pos 2600
    mid2=$(make_filler 1000)
    qarm=$(repeat_motif "TTAGGG" 100)       # 600bp q-arm
    echo "${parm}${mid1}${its_telo}${mid2}${qarm}"
} > "$DIR/boundary_extend_its.fa"

echo "Generated $(ls -1 "$DIR"/*.fa "$DIR"/*.gfa 2>/dev/null | wc -l) synthetic test files in $DIR/"
