# shellcheck shell=bash

# Fixture declarations for generate_synthetic.sh. Sourced, never run directly.

# Declaration format: docs/testing.md, section Fixtures.

# <flags>: "-" means defaults; otherwise the exact flags the expectation was checked under.

# Legacy fixtures: byte-identical to what shipped; consumed by val.sh and test_gaps_bed.sh.

fx t2t t2t.fa \
   'chr_t2t=F:CCCTAAx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Both arms, correct orientation'

fx incomplete_p incomplete_p.fa \
   'chr_incomplete_p=F:CCCTAAx100+L:6600' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Only p-arm telomere'

fx incomplete_q incomplete_q.fa \
   'chr_incomplete_q=L:6600+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'Only q-arm telomere'

fx no_telo no_telo.fa \
   'chr_none=L:3000' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'No telomeres'

# -d 200 keeps density-trim from bridging; --link-distance (1000) chains them anyway.
fx misassembly misassembly.fa \
   'chr_misassembly=F:CCCTAAx100+L:500+F:CCCTAAx100+L:5000' '-d 200' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;telolen=1200' \
   'Two p arrays 500 bp apart exceed -d 200 for density-trim but chain within --link-distance into one fragmented arm'

# Filler > tolerance 3000, so the p end never reaches this array too and reclaims it (R3).
fx discordant discordant.fa \
   'chr_discordant=L:7200+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=0' \
   'Forward motif sitting at the q end: arm and strand disagree'

fx gapped_t2t gapped_t2t.fa \
   'chr_gapped_t2t=F:CCCTAAx100+L:900+N:100+L:900+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=1' \
   'T2T with an internal N-gap; gappedness is a separate column, not a type'

fx multi multi.fa \
   'contig_t2t=F:CCCTAAx100+L:2000+R:TTAGGGx100;contig_none=L:2000;contig_incomplete=L:6600+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'Three records: t2t, none, incomplete'

fx plant plant.fa \
   'chr_plant=F:CCCTAAAx86+L:2000+R:TTTAGGGx86' '-c CCCTAAA' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Plant 7-mer canonical repeat'

# Terminal density is canonical-only (D10): a variant array never anchors a piece, at any -x.
fx edit_test edit_test.fa \
   'chr_edit_p=V:CTCTAAx100+L:2400' '-x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'A one-substitution-per-repeat variant array is never a terminal anchor, at -x 0 or -x 1'

fx short_contig short_contig.fa \
   'chr_tiny=L:100' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'Contig shorter than one window'

# Interleaved orientations never qualify as a chain piece (R1/R3): moved to the interior.
fx balanced balanced.fa \
   'chr_balanced=L:3500+M:CCCTAATTAGGGx100+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   'Interleaved forward and reverse repeats: no qualifying arm, one interior interstitial b row'

fx its its.fa \
   'chr_its=F:CCCTAAx100+L:3000+R:TTAGGGx100+L:3000+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'Telomeres at both ends plus one interstitial array in the middle'

fx density_edge density_edge.fa \
   'chr_density=M:CCCTAAACGATCx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'p-arm alternating 6 bp match and 6 bp filler: density exactly at -y, kept on the tie' \

fx misassembly_qq misassembly_qq.fa \
   'chr_misassembly_qq=L:5000+L:500+R:TTAGGGx100+L:500+R:TTAGGGx100' '-d 200' \
   'type=incomplete;anom=fragmented_q;gran=Q;telo=1;labels=q;gaps=0;telolen=1200' \
   'Two q arrays 500 bp apart exceed -d 200 for density-trim but chain within --link-distance into one fragmented arm'

fx gapped_misassembly gapped_misassembly.fa \
   'chr_gapped_misassembly=F:CCCTAAx100+L:500+F:CCCTAAx100+L:2000+N:100+L:3000' '-d 200' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=1;telolen=1200' \
   'Two p arrays 500 bp apart chain into one fragmented arm; an unrelated N-gap sits in the tail contig'

fx gapped_incomplete gapped_incomplete.fa \
   'chr_gapped_incomplete=F:CCCTAAx100+L:1000+N:100+L:2000' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1' \
   'Single p-arm with an N-gap'

fx gapped_none gapped_none.fa \
   'chr_gapped_none=L:1000+N:100+L:2000' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1' \
   'No telomeres but one N-gap'

fx gapped_discordant gapped_discordant.fa \
   'chr_gapped_discordant=L:2000+N:100+L:500+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=1' \
   'Forward motif at the q end, with an N-gap'

fx mirror_inverted_short mirror_inverted_short.fa \
   'chr_mirror_inv_short=R:TTAGGGx100+L:2000+F:CCCTAAx100' '-' \
   'type=t2t;anom=discordant_p,discordant_q;gran=P*Q*;telo=2;labels=pq;gaps=0' \
   'Both arms present but each carries the other end motif: an inverted terminal repeat'

# Filler > tolerance 3000, so the q end never reaches this array too and reclaims it (R3).
fx mirror_rev_start mirror_rev_start.fa \
   'chr_mirror_rev_start=R:TTAGGGx100+L:7200' '-' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=0' \
   'Reverse motif at the p end: the mirror of discordant.fa'

fx mirror_inverted_long mirror_inverted_long.fa \
   'chr_mirror_inv_long=R:TTAGGGx100+L:4000+F:CCCTAAx100' '-t 1000' \
   'type=t2t;anom=discordant_p,discordant_q;gran=P*Q*;telo=2;labels=pq;gaps=0' \
   'Inverted arms on a contig longer than twice the terminal zone'

fx mirror_fwd_end_long mirror_fwd_end_long.fa \
   'chr_mirror_fwd_end_long=L:4600+F:CCCTAAx100' '-t 1000' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=0' \
   'Forward motif only, at the far end'

fx mirror_rev_start_long mirror_rev_start_long.fa \
   'chr_mirror_rev_start_long=R:TTAGGGx100+L:4600' '-t 1000' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=0' \
   'Reverse motif only, at the start'

# A real reverse array right behind the p arm ends its chain (R1); filler > tolerance keeps
# it out of reach from the q end, so it lands as a tail_to_tail interstitial row (D11).
fx mirror_both_start mirror_both_start.fa \
   'chr_mirror_both_start=F:CCCTAAx100+R:TTAGGGx100+L:7000' '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0;its=1' \
   'An abutting reverse array ends the p chain at once and is reported as an interstitial row'

# The mirror: a real forward array ahead of the q arm ends its chain, out of reach from p.
fx mirror_both_end mirror_both_end.fa \
   'chr_mirror_both_end=L:7000+F:CCCTAAx100+R:TTAGGGx100' '-i' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0;its=1' \
   'An abutting forward array ends the q chain at once and is reported as an interstitial row'

fx mirror_extend mirror_extend.fa \
   'chr_mirror_extend=F:CCCTAAx100+L:2400' '-t 300' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Terminal block starts inside the zone and extends past it'

fx mirror_merge mirror_merge.fa \
   'chr_mirror_merge=F:CCCTAAx50+L:100+F:CCCTAAx50+L:6500' '-d 200' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Two clusters closer than -d merge into one block'

# -d 200 keeps the two 240 bp clusters separate; each is below -l on its own.
fx mirror_no_merge mirror_no_merge.fa \
   'chr_mirror_no_merge=F:CCCTAAx40+L:300+F:CCCTAAx40+L:2000' '-d 200' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'Two clusters further apart than -d, each below -l: neither survives'

fx boundary_no_terminal boundary_no_terminal.fa \
   'chr_boundary_no_term=L:2000+F:CCCTAAx100+L:2000+R:TTAGGGx100+L:2000' '-i -t 500' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'Both arrays outside the terminal zone: no terminal blocks, two interstitial ones'

fx boundary_fail_filter boundary_fail_filter.fa \
   'chr_boundary_fail=F:CCCTAAx20+L:3000' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'Terminal array of 120 bp is below -l and is dropped'

fx boundary_cross boundary_cross.fa \
   'chr_boundary_cross=F:CCCTAAx200+R:TTAGGGx200' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=0' \
   'Dense abutting p and q arrays leave no room between them for an interstitial block'

fx boundary_zone_shift boundary_zone_shift.fa \
   'chr_boundary_zone=L:500+F:CCCTAAx100+L:2000+R:TTAGGGx100+L:500' '-t 1200' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Zone size decides terminal versus interstitial: at -t 1200 both arrays are terminal'

fx boundary_zone_shift_narrow boundary_zone_shift.fa \
   'chr_boundary_zone=L:500+F:CCCTAAx100+L:2000+R:TTAGGGx100+L:500' '-t 400' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'The same file at -t 400: both arrays fall outside the zone'

fx boundary_its_at_edge boundary_its_at_edge.fa \
   'chr_boundary_edge=F:CCCTAAx100+R:TTAGGGx10+L:2200+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'A 60 bp array abutting the p arm is too short to be a real array (R1) and does not clamp the chain, but it clears the interstitial gate: no length floor there any more'

fx boundary_multiple_p boundary_multiple_p.fa \
   'chr_boundary_multi_p=F:CCCTAAx100+L:300+F:CCCTAAx100+L:5700' '-d 200' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;telolen=1200' \
   'Two p arrays 300 bp apart exceed -d 200 for density-trim but chain within --link-distance into one fragmented arm'

# The 1000 bp gap is exactly --link-distance (inclusive), so the two q pieces chain into
# one fragmented arm regardless of -t; nothing is left over for the interstitial file.
fx boundary_extend_its boundary_extend_its.fa \
   'chr_boundary_ext_its=F:CCCTAAx100+L:2000+R:TTAGGGx50+L:1000+R:TTAGGGx100' '-i -t 300' \
   'type=t2t;anom=fragmented_q;gran=PQ;telo=2;labels=pq;gaps=0;its=0;telolen=900' \
   'Terminal blocks extend past the -t zone: the middle array chains into the q arm, not interstitial'

# Nothing bounds the extent (R2): a 1200 bp array is reported whole even at -t 300, and the
# fast-scan window extension (R6) makes fast and full scan agree.
fx uncapped_extent_fast uncapped_extent.fa \
   'chr_uncapped_extent=F:CCCTAAx200+L:5000' '-t 300' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0;telolen=1200' \
   'A 1200 bp array is not clipped at the -t 300 window: fast mode extends past it'

fx uncapped_extent_full uncapped_extent.fa \
   'chr_uncapped_extent=F:CCCTAAx200+L:5000' '-i -t 300' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0;telolen=1200;its=0' \
   'The same file under a full scan: identical to the fast-mode result'

# -t 50 shrinks the zone to almost nothing; with -k 200 -x 1 both halves stay interstitial.
# -n is dropped here: at -t 50 each half also anchors its own contig's internal end, which
# would move it to the terminal BED instead of testing the plain interstitial builder.
fx its_gap_split its_gap_split.fa \
   'chr_its_gap_split=L:100+F:CCCTAAx100+N:100+R:TTAGGGx100+L:100' '-i -t 50 -k 200 -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=2' \
   'One array split by a gap into two contigs; each half is its own interstitial row, classed single'

fx its_gap_headtohead its_gap_headtohead.fa \
   'chr_its_gap_headtohead=L:100+R:TTAGGGx100+N:100+F:CCCTAAx100+L:100' '-i -t 50 -k 200 -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=2' \
   'Head-to-head arrays split by a gap: two independent contigs, each an interstitial row classed single (a gap breaks any junction link, D11)'

fx its_strand_pure its_strand_pure.fa \
   'chr_its_strand_pure=L:100+V:TTAGGAx7+F:CCCTAAx13+L:100' '-i -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   'Reverse variant halo fails the exact-canonical-count gate on its own (inversion cut splits it from the core); only the pure forward core is reported'

fx its_headtohead its_headtohead.fa \
   'chr_its_headtohead=L:100+F:CCCTAAx20+V:CTCTAAx7+V:TTAGGAx7+R:TTAGGGx20+L:100' '-i -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'F and R canonical runs joined by a variant halo flip orientation partway; the inversion cut still splits them into two rows, a head-to-head pair (D11)'

# Shipped before this generator (used by generate-tests.cpp, test_gaps_bed.sh); bytes are fixed.

fx discordant_pp discordant_pp.fa \
   'chr_discordant_pp=L:2996+F:CCCTAAx84+L:900+F:CCCTAAx100' '-d 200' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;telolen=1104' \
   'Both forward arrays are reachable from either end; strand decides (R3): one fragmented p arm, no q telomere'

fx extra_invalid_p extra_invalid_p.fa \
   'chr_extra_invalid_p=F:CCCTAAx100+L:3396+F:CCCTAAx84+L:500' '-d 200' \
   'type=t2t;anom=discordant_q;gran=PQ*;telo=2;labels=pq;gaps=0' \
   'A p arm plus a forward array at the q end, which reads as a discordant q arm'

fx gapped_discordant_q gapped_discordant_q.fa \
   'chr_gapped_discordant_q=R:TTAGGGx100+L:2000+N:100+L:2300' '-' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=1' \
   'Reverse motif at the p end, with an N-gap'

fx gapped_incomplete_q gapped_incomplete_q.fa \
   'chr_gapped_incomplete_q=L:2000+N:100+L:2300+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=1' \
   'Single q-arm with an N-gap'

fx gapped_misassembly_qq gapped_misassembly_qq.fa \
   'chr_gapped_misassembly_qq=L:500+N:100+L:2200+R:TTAGGGx100+L:496+R:TTAGGGx84+L:600' '-d 200' \
   'type=incomplete;anom=fragmented_q;gran=Q;telo=1;labels=q;gaps=1;telolen=1104' \
   'Two q arrays 496 bp apart exceed -d 200 for density-trim but chain within --link-distance into one fragmented arm; an unrelated leading N-gap'

fx multi_gap_t2t multi_gap_t2t.fa \
   'chr_multi_gap_t2t=F:CCCTAAx100+L:400+N:50+L:400+N:50+L:400+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=2' \
   'T2T with two internal N-gaps'

_gfa_parm=$(repeat_motif "$CANON_FWD" 100)
_gfa_qarm=$(repeat_motif "$CANON_REV" 100)
_gfa_fill=$(make_filler 2000)

gfx gfa_telo gfa_telo.gfa \
"H	VN:Z:1.2
S	seg_t2t	${_gfa_parm}${_gfa_fill}${_gfa_qarm}
S	seg_ponly	${_gfa_parm}${_gfa_fill}
S	seg_qonly	${_gfa_fill}${_gfa_qarm}
S	seg_none	${_gfa_fill}
L	seg_t2t	+	seg_ponly	+	0M
L	seg_ponly	+	seg_qonly	+	0M
L	seg_qonly	+	seg_none	+	0M
" \
   'Four segments: t2t, p-only, q-only, no-telo, chained by links'

unset _gfa_parm _gfa_qarm _gfa_fill

# Real genomes: not owned, never written or checked; carried in the manifest for the intent harness.

xfx bTaeGut7_mat bTaeGut7_chr33_mat.fa.gz chr33_mat '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Real zebra finch chr33, maternal haplotype'

xfx bTaeGut7_pat bTaeGut7_chr33_pat.fa.gz chr33_pat '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Real zebra finch chr33, paternal haplotype'

# VGP excerpts from the 26.09.12 edge-case panel (notebook section in brackets).
xfx vgp_probe_mega vgp_probe.fa.gz mega_OZ124247.1_p_0-600000 '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0;its=0;telolen=538251' \
   'Pochard OZ124247.1 first 600 kb: a single uncapped p arm, 4-538255 [S4.1, S6.2 C01, REG-018]'

xfx vgp_probe_frag32 vgp_probe.fa.gz frag32_OZ221982.1_q_302227075-302307075 '-i' \
   'type=incomplete;anom=fragmented_q;gran=Q;telo=1;labels=q;gaps=0;its=1;telolen=44993' \
   'Frog OZ221982.1 last 80 kb: pieces chain into one fragmented q arm, 0-79163 [S4.2, S7.1]'

xfx vgp_probe_mito vgp_probe.fa.gz mito_CM010492.2_whole '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=0' \
   'A 16 kb mitochondrial genome: no telomere'

xfx vgp_turtle vgp_turtle.fa.gz giant_CM098529.1_84917045-85927174 '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=3' \
   'Turtle CM098529.1: a near-1 Mb reverse array on the first contig, too far from either scaffold end to be an arm, reported as three interstitial rows; never bridged across the 200 bp gap [S4.3]'

xfx vgp_turtle_n vgp_turtle.fa.gz giant_CM098529.1_84917045-85927174 '-i -n' \
   'type=none;anom=.;gran=q;telo=0;labels=none;gaps=1;its=0' \
   'Same file with -n: the array anchors its own contig'"'"'s internal end and becomes a lowercase contig row, leaving the interstitial BED'

# Combinatorial axes
source "$SCRIPT_DIR/synthetic_axis_anomaly.sh"
source "$SCRIPT_DIR/synthetic_axis_threshold.sh"
source "$SCRIPT_DIR/synthetic_axis_structure.sh"
source "$SCRIPT_DIR/synthetic_axis_more.sh"
source "$SCRIPT_DIR/synthetic_axis_junction.sh"
