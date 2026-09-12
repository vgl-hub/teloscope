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
   'chr_incomplete_p=F:CCCTAAx100+L:2400' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Only p-arm telomere'

fx incomplete_q incomplete_q.fa \
   'chr_incomplete_q=L:2400+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'Only q-arm telomere'

fx no_telo no_telo.fa \
   'chr_none=L:3000' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'No telomeres'

# -d 200 keeps the two same-end blocks separate; the default -d 500 would merge them (val.sh:133).
fx misassembly misassembly.fa \
   'chr_misassembly=F:CCCTAAx100+L:500+F:CCCTAAx100+L:5000' '-d 200' \
   'type=incomplete;anom=misassembly;gran=Pp;telo=1;labels=p;gaps=0' \
   'Two p-arms at the same end: one extra terminal block'

fx discordant discordant.fa \
   'chr_discordant=L:2800+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=0' \
   'Forward motif sitting at the q end: arm and strand disagree'

fx gapped_t2t gapped_t2t.fa \
   'chr_gapped_t2t=F:CCCTAAx100+L:900+N:100+L:900+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=1' \
   'T2T with an internal N-gap; gappedness is a separate column, not a type'

fx multi multi.fa \
   'contig_t2t=F:CCCTAAx100+L:2000+R:TTAGGGx100;contig_none=L:2000;contig_incomplete=L:2400+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'Three records: t2t, none, incomplete'

fx plant plant.fa \
   'chr_plant=F:CCCTAAAx86+L:2000+R:TTTAGGGx86' '-c CCCTAAA' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Plant 7-mer canonical repeat'

# Density gate (teloscope.cpp:403) counts canonical coverage only; uncallable at any -y, REG-013.
fx edit_test edit_test.fa \
   'chr_edit_p=V:CTCTAAx100+L:2000+R:TTAGGGx100' '-x 1' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'One substitution per p-arm repeat: found at -x 1, not at -x 0'

fx short_contig short_contig.fa \
   'chr_tiny=L:100' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'Contig shorter than one window'

fx balanced balanced.fa \
   'chr_balanced=L:1200+M:CCCTAATTAGGGx100+L:1200' '-' \
   'type=incomplete;anom=balanced_p;gran=P~;telo=1;labels=p;gaps=0' \
   'Interleaved forward and reverse repeats give one array a balanced strand label'

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
   'type=incomplete;anom=misassembly;gran=qQ;telo=1;labels=q;gaps=0' \
   'Two q-arms at the same end: one extra terminal block'

fx gapped_misassembly gapped_misassembly.fa \
   'chr_gapped_misassembly=F:CCCTAAx100+L:500+F:CCCTAAx100+L:2000+N:100+L:3000' '-d 200' \
   'type=incomplete;anom=misassembly;gran=Pp;telo=1;labels=p;gaps=1' \
   'Two p-arms plus an unrelated N-gap in the tail'

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

fx mirror_rev_start mirror_rev_start.fa \
   'chr_mirror_rev_start=R:TTAGGGx100+L:2400' '-' \
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

fx mirror_both_start mirror_both_start.fa \
   'chr_mirror_both_start=F:CCCTAAx100+R:TTAGGGx100+L:2000' '-' \
   'type=incomplete;anom=misassembly;gran=Pp*;telo=1;labels=p;gaps=0' \
   'Abutting opposite orientations are cut apart by inversionRun, giving two p blocks'

fx mirror_both_end mirror_both_end.fa \
   'chr_mirror_both_end=L:2000+F:CCCTAAx100+R:TTAGGGx100' '-' \
   'type=incomplete;anom=misassembly;gran=q*Q;telo=1;labels=q;gaps=0' \
   'Same at the q end; the >= tie-break elects the outermost block as the q arm'

fx mirror_extend mirror_extend.fa \
   'chr_mirror_extend=F:CCCTAAx100+L:2400' '-t 300' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Terminal block starts inside the zone and extends past it'

fx mirror_merge mirror_merge.fa \
   'chr_mirror_merge=F:CCCTAAx50+L:100+F:CCCTAAx50+L:2000' '-d 200' \
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
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=0' \
   'A 60 bp array abutting the p arm is below --min-its-length and is not reported'

fx boundary_multiple_p boundary_multiple_p.fa \
   'chr_boundary_multi_p=F:CCCTAAx100+L:300+F:CCCTAAx100+L:2000' '-d 200' \
   'type=incomplete;anom=misassembly;gran=Pp;telo=1;labels=p;gaps=0' \
   'Two p arrays further apart than -d stay separate: the inner one is an extra block'

fx boundary_extend_its boundary_extend_its.fa \
   'chr_boundary_ext_its=F:CCCTAAx100+L:2000+R:TTAGGGx50+L:1000+R:TTAGGGx100' '-i -t 300' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'Terminal blocks extend past the zone; the middle array stays interstitial'

fx its_gap_split its_gap_split.fa \
   'chr_its_gap_split=L:100+F:CCCTAAx100+N:100+R:TTAGGGx100+L:100' '-i -n -t 50 -k 200 -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=?' \
   'One array split by a gap; both halves are too far from an end to be terminal'

fx its_gap_headtohead its_gap_headtohead.fa \
   'chr_its_gap_headtohead=L:100+R:TTAGGGx100+N:100+F:CCCTAAx100+L:100' '-i -n -t 50 -k 200 -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=?' \
   'Head-to-head arrays across a gap: see docs/conflicts.md REG-004'

fx its_strand_pure its_strand_pure.fa \
   'chr_its_strand_pure=L:100+V:TTAGGAx7+F:CCCTAAx13+L:100' '-i -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=?' \
   'Reverse variant halo around a pure canonical forward core: pooled label is balanced'

fx its_headtohead its_headtohead.fa \
   'chr_its_headtohead=L:100+F:CCCTAAx20+V:CTCTAAx7+V:TTAGGAx7+R:TTAGGGx20+L:100' '-i -x 1' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=?' \
   'All four canonicality and orientation count cells are non-zero'

# Shipped before this generator (used by generate-tests.cpp, test_gaps_bed.sh); bytes are fixed.

fx discordant_pp discordant_pp.fa \
   'chr_discordant_pp=L:2996+F:CCCTAAx84+L:900+F:CCCTAAx100' '-d 200' \
   'type=incomplete;anom=discordant_q,misassembly;gran=q*Q*;telo=1;labels=q;gaps=0' \
   'Two forward arrays both at the q end: discordant and one extra block'

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
   'type=incomplete;anom=misassembly;gran=Qq;telo=1;labels=q;gaps=1' \
   'Two q arrays plus an unrelated N-gap'

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

# Combinatorial axes
source "$SCRIPT_DIR/synthetic_axis_anomaly.sh"
source "$SCRIPT_DIR/synthetic_axis_threshold.sh"
source "$SCRIPT_DIR/synthetic_axis_structure.sh"
source "$SCRIPT_DIR/synthetic_axis_more.sh"
source "$SCRIPT_DIR/synthetic_axis_junction.sh"
