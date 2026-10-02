# shellcheck shell=bash

# Combinatorial tests: incomplete, seed boundary, window invariance.

# ---- incomplete x anomaly x gappedness ----
fx mo_inc_p_clean synthetic/mo_inc_p_clean.fa \
   'chr_mo_inc_p_clean=F:CCCTAAx100+L:6600' '-' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0' \
   'One concordant p arm'

fx mo_inc_p_clean_gapped synthetic/mo_inc_p_clean_gapped.fa \
   'chr_mo_inc_p_clean_gapped=F:CCCTAAx100+L:1300+N:200+L:1300' '-' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=1' \
   'One concordant p arm on a gapped scaffold'

fx mo_inc_q_clean synthetic/mo_inc_q_clean.fa \
   'chr_mo_inc_q_clean=L:6600+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;telo=1;labels=q;gaps=0' \
   'One concordant q arm'

# Filler > tolerance 3000, so the p end never reaches this array too and reclaims it.
fx mo_inc_q_disc synthetic/mo_inc_q_disc.fa \
   'chr_mo_inc_q_disc=L:7200+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;telo=1;labels=q;gaps=0' \
   'One inverted q arm'

# Interleaved orientations never qualify as a chain piece: no telomere at either end.
fx mo_inc_p_bal_gapped synthetic/mo_inc_p_bal_gapped.fa \
   'chr_mo_inc_p_bal_gapped=M:CCCTAATTAGGGx50+L:1300+N:200+L:1300' '-' \
   'type=none;anom=.;telo=0;labels=none;gaps=1' \
   'An interleaved array on the first contig and a pure-filler last contig: no telomere at either end'

fx mo_inc_q_disc_gapped synthetic/mo_inc_q_disc_gapped.fa \
   'chr_mo_inc_q_disc_gapped=L:1300+N:200+L:1300+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;telo=1;labels=q;gaps=1' \
   'One inverted q arm on a gapped scaffold'

# Two same-strand arrays within -d are two pieces of one fragmented arm.
fx mo_inc_p_extra synthetic/mo_inc_p_extra.fa \
   'chr_mo_inc_p_extra=F:CCCTAAx100+L:600+F:CCCTAAx100+L:5400' '-' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;telolen=1200' \
   'Two p arrays 600 bp apart are two pieces of one fragmented p arm (teloLen 1200); no q arm'

fx mo_inc_q_extra synthetic/mo_inc_q_extra.fa \
   'chr_mo_inc_q_extra=L:5400+R:TTAGGGx100+L:600+R:TTAGGGx100' '-' \
   'type=incomplete;anom=fragmented_q;telo=1;labels=q;gaps=0;telolen=1200' \
   'Two q arrays 600 bp apart are two pieces of one fragmented q arm (teloLen 1200); no p arm'

# Each 600 bp array outweighs the 500 bp gap before it, so all three join.
fx mo_inc_p_three_piece synthetic/mo_inc_p_three_piece.fa \
   'chr_mo_inc_p_three_piece=F:CCCTAAx100+L:500+F:CCCTAAx100+L:500+F:CCCTAAx100+L:4400' '-' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;telolen=1800' \
   'Three 600 bp p arrays 500 bp apart each outweigh the gap before them and join: one fragmented p arm, teloLen 1800, hull 0-2800'

# Filler > tolerance 3000 keeps the q end from reaching and reclaiming this structure.
fx mo_inc_p_extra_disc synthetic/mo_inc_p_extra_disc.fa \
   'chr_mo_inc_p_extra_disc=R:TTAGGGx100+L:600+R:TTAGGGx100+L:6600' '-' \
   'type=incomplete;anom=discordant_p,fragmented_p;telo=1;labels=p;gaps=0;telolen=1200' \
   'Two inverted arrays 600 bp apart are two pieces of one arm that is both discordant and fragmented'

fx mo_inc_p_extra_bal synthetic/mo_inc_p_extra_bal.fa \
   'chr_mo_inc_p_extra_bal=M:CCCTAATTAGGGx50+L:1200+F:CCCTAAx100+L:4800' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=1' \
   'Interleaved orientations at the very start fail every anchor attempt and are skipped; the plain forward array 1200 bp in, beyond -d, is the p arm, and the interleaved run is a single interstitial b row'

# The plain array outweighs the 506 bp between the first piece's last exact repeat and its start.
fx mo_elect_longer_loses synthetic/mo_elect_longer_loses.fa \
   'chr_mo_elect_longer_loses=M:CCCTAAACGATCx60+L:500+F:CCCTAAx100+L:5380' '-' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;telolen=1314' \
   'A dense alternating array and a plain array 500 bp behind it are two pieces of one fragmented arm, the plain 600 bp piece outweighing the 506 bp gap; the first piece ends at its last exact repeat (714 bp), not the trailing filler unit'

fx mo_none_gapped_two synthetic/mo_none_gapped_two.fa \
   'chr_mo_none_gapped_two=L:900+N:150+L:900+N:150+L:900' '-' \
   'type=none;anom=.;telo=0;labels=none;gaps=2' \
   'No telomeres and two gaps'

fx mo_t2t_gapped_three synthetic/mo_t2t_gapped_three.fa \
   'chr_mo_t2t_gapped_three=F:CCCTAAx100+L:600+N:100+L:600+N:100+L:600+N:100+L:600+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;telo=2;labels=pq;gaps=3' \
   'A t2t scaffold broken into four contigs'

# Each array is a -k seed of its own; -d only sets the junction class of the two rows.
fx mo_seed_gap_exact synthetic/mo_seed_gap_exact.fa \
   'chr_mo_seed_gap_exact=L:3500+F:CCCTAAx200+L:1000+F:CCCTAAx200+L:3500' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=2' \
   'Two interior arrays 1000 bp apart are separate -k seeds: two interstitial rows, the second exactly -d from the first and so classed fragmentation'

fx mo_seed_gap_beyond synthetic/mo_seed_gap_beyond.fa \
   'chr_mo_seed_gap_beyond=L:3500+F:CCCTAAx200+L:1001+F:CCCTAAx200+L:3500' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=2' \
   'A 1001 bp run exceeds -d: the two interior arrays stay two interstitial rows, classed single'

fx mo_window_overlapping synthetic/mo_window_overlapping.fa \
   'chr_mo_window_overlapping=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-r -w 100 -s 50' \
   'type=t2t;anom=.;telo=2;labels=pq;gaps=0' \
   'Overlapping windows change the tracks, never the call'

fx mo_window_large synthetic/mo_window_large.fa \
   'chr_mo_window_large=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-r -w 2000 -s 2000' \
   'type=t2t;anom=.;telo=2;labels=pq;gaps=0' \
   'A window wider than either arm changes the tracks, never the call'

fx mo_long_scaffold_default synthetic/mo_long_scaffold.fa \
   'chr_mo_long_scaffold=F:CCCTAAx100+L:20000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;telo=2;labels=pq;gaps=0' \
   'A 21 kb scaffold is well inside the default terminal limit'

fx mo_long_scaffold_narrow synthetic/mo_long_scaffold_narrow.fa \
   'chr_mo_long_scaffold_narrow=L:3000+F:CCCTAAx100+L:20000+R:TTAGGGx100+L:3000' '-t 2000' \
   'type=none;anom=.;telo=0;labels=none;gaps=0' \
   'Both arrays sit 3 kb in, outside a 2 kb terminal zone'

fx mo_order_forward synthetic/mo_order_forward.fa \
   'rec_a=F:CCCTAAx100+L:2800+R:TTAGGGx100;rec_b=L:3000;rec_c=L:6600+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;telo=2;labels=pq;gaps=0|type=none;anom=.;telo=0;labels=none;gaps=0|type=incomplete;anom=.;telo=1;labels=q;gaps=0' \
   'Three records in one order'

fx mo_order_reversed synthetic/mo_order_reversed.fa \
   'rec_c=L:6600+R:TTAGGGx100;rec_b=L:3000;rec_a=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;telo=1;labels=q;gaps=0|type=none;anom=.;telo=0;labels=none;gaps=0|type=t2t;anom=.;telo=2;labels=pq;gaps=0' \
   'The same three records in the opposite order get the same per-record answers'

fx mo_sparse_long synthetic/mo_sparse_long.fa \
   'chr_mo_sparse_long=M:CCCTAAACGATCx100+L:6000' '-' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0' \
   'A long array at exactly -y density survives on both length and density'

fx mo_dense_short synthetic/mo_dense_short.fa \
   'chr_mo_dense_short=F:CCCTAAx40+L:2760' '-' \
   'type=none;anom=.;telo=0;labels=none;gaps=0' \
   'A fully dense 240 bp array still fails -l: density does not buy length'

# ---- piece rule: -k chains, exact-repeat score, pieces joined within -d, -l summed over the pieces ----

# The score only counts exact repeats, so a variant-only tail chained to the array costs more than its few exact repeats earn; the tail is left to the interstitial rows.
fx pc_chain_tail synthetic/pc_chain_tail.fa \
   'chr_pc_chain_tail=F:CCCTAAx100+V:CTCTAAx40+F:CCCTAAx1+V:CTCTAAx40+F:CCCTAAx1+V:CTCTAAx40+F:CCCTAAx1+V:CTCTAAx40+F:CCCTAAx1+V:CTCTAAx40+F:CCCTAAx1+L:5000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=1;telolen=600' \
   'A 600 bp exact array chained to a 1230 bp tail of one-mismatch repeats with five exact ones: the arm is the exact array alone, the tail is an interstitial row'

# Fewer than --min-block-counts exact repeats make no piece, however regular.
fx pc_single_every200 synthetic/pc_single_every200.fa \
   'chr_pc_single_every200=F:CCCTAAx100+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:194+F:CCCTAAx1+L:5000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=0;telolen=600' \
   'Single exact repeats every 200 bp behind a 600 bp array are no pieces (one exact repeat each): the arm is the array alone, not fragmented, and no row is left'

fx pc_pairs_every200 synthetic/pc_pairs_every200.fa \
   'chr_pc_pairs_every200=F:CCCTAAx100+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:188+F:CCCTAAx2+L:5000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=0;telolen=600' \
   'Pairs of exact repeats every 200 bp behind a 600 bp array are pieces, but each 12 bp pair is shorter than the gap in front of it and none joins: the arm is the array alone, and a pair holds too few repeats for an interstitial row'

# Pieces add up to -l: four 60 bp arrays 60 bp apart sum to 240, five to 300.
fx pc_tiny_below synthetic/pc_tiny_below.fa \
   'chr_pc_tiny_below=F:CCCTAAx10+L:60+F:CCCTAAx10+L:60+F:CCCTAAx10+L:60+F:CCCTAAx10+L:5000' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=4' \
   'Four 60 bp arrays, each as long as the 60 bp gap before it, join into 240 bp, short of -l: no telomere, and each array is an interstitial row'

fx pc_tiny_sum synthetic/pc_tiny_sum.fa \
   'chr_pc_tiny_sum=F:CCCTAAx10+L:60+F:CCCTAAx10+L:60+F:CCCTAAx10+L:60+F:CCCTAAx10+L:60+F:CCCTAAx10+L:5000' '-i' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;its=0;telolen=300' \
   'Five 60 bp arrays, each as long as the 60 bp gap before it, all join and sum to 300 bp, reaching -l: one fragmented p arm with hull 0-540'

# A variant lead-in is scored as non-repeat: the exact array behind it must outweigh it for the piece to reach the array.
fx pc_leadin_heavy synthetic/pc_leadin_heavy.fa \
   'chr_pc_leadin_heavy=F:CCCTAAx2+V:CTCTAAx500+F:CCCTAAx600+L:6000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=0;telolen=6612' \
   'Two exact repeats at the tip, 3000 bp of one-mismatch repeats, then a 3600 bp exact array: the array outweighs the lead-in and the piece starts at the tip and runs through it (6612 bp)'

fx pc_leadin_light synthetic/pc_leadin_light.fa \
   'chr_pc_leadin_light=F:CCCTAAx2+V:CTCTAAx500+F:CCCTAAx100+L:6000' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=1' \
   'The same lead-in with a 600 bp exact array: it does not outweigh the variants, the tip piece is 12 bp and dropped, and the array starts at 3012, beyond the 3000 bp start zone: no arm, one interstitial row'

# -k is inclusive: a match starting exactly -k past the chain end joins it.
fx pc_k_50 synthetic/pc_k_50.fa \
   'chr_pc_k_50=F:CCCTAAx50+L:50+F:CCCTAAx50+L:6000' '-' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;telolen=650' \
   'Two 300 bp arrays 50 bp apart are one -k chain and one piece (300 + 300 - 50 filler = 550 at its peak, hull 0-650): not fragmented, teloLen equals the hull'

fx pc_k_51 synthetic/pc_k_51.fa \
   'chr_pc_k_51=F:CCCTAAx50+L:51+F:CCCTAAx50+L:6000' '-' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;telolen=600' \
   'Two 300 bp arrays 51 bp apart are two chains: two pieces within -d, so a fragmented p arm with teloLen 600 over hull 0-651'

# -l is compared on the sum of the pieces.
fx pc_sum_300 synthetic/pc_sum_300.fa \
   'chr_pc_sum_300=F:CCCTAAx25+L:150+F:CCCTAAx25+L:6050' '-i' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;its=0;telolen=300' \
   'Two 150 bp arrays 150 bp apart join and sum to exactly -l 300: a fragmented p arm'

fx pc_sum_294 synthetic/pc_sum_294.fa \
   'chr_pc_sum_294=F:CCCTAAx25+L:144+F:CCCTAAx24+L:6056' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=2' \
   'Two arrays of 150 and 144 bp, the second as long as the 144 bp gap before it, join into 294 bp, short of -l: no telomere, two interstitial rows'

# A chain that falls short of -l is dropped and the search resumes behind its first piece, still inside the start zone.
fx pc_restart_stub synthetic/pc_restart_stub.fa \
   'chr_pc_restart_stub=R:TTAGGGx8+L:1500+F:CCCTAAx500+L:6000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=1;telolen=3000' \
   'A 48 bp reverse stub at the tip falls short of -l and is dropped; the 3 kb forward array starting at 1548, inside the start zone, becomes a concordant p arm, and the stub is an interstitial row'

# Interior head-to-head junction: two rows per side, 4 rows in all; their junction classes are checked by check_invariants.py DER-20.
fx pc_h2h_pairs synthetic/pc_h2h_pairs.fa \
   'chr_pc_h2h_pairs=L:4000+R:TTAGGGx60+L:100+R:TTAGGGx60+L:20+F:CCCTAAx60+L:100+F:CCCTAAx60+L:4000' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=4' \
   'Two reverse arrays 100 bp apart, then two forward arrays 20 bp after them and 100 bp apart, all far from both ends: no arm and four interstitial rows (the inversion cut splits the 20 bp junction)'

# ---- piece rule: a later piece joins only when it is at least as long as the gap behind the last joined piece, times y/(1-y) ----

# At -y 0.5 the weight is 1: a 300 bp array exactly 300 bp behind the 600 bp array joins.
fx pc_gap_300_joins synthetic/pc_gap_300_joins.fa \
   'chr_pc_gap_300_joins=F:CCCTAAx100+L:300+F:CCCTAAx50+L:6000' '-i' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;its=0;telolen=900' \
   'A 300 bp array exactly as long as the 300 bp gap before it joins the 600 bp array: one fragmented p arm of teloLen 900, hull 0-1200, no row left over'

fx pc_gap_301_skips synthetic/pc_gap_301_skips.fa \
   'chr_pc_gap_301_skips=F:CCCTAAx100+L:301+F:CCCTAAx50+L:6000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=1;telolen=600' \
   'A 300 bp array one base shorter than the 301 bp gap before it does not join: the p arm is the 600 bp array alone, and the 300 bp array is an interstitial row classed fragmentation by the p arm within -d'

# The skipped 48 bp piece does not end the search: the next piece is measured from the last joined piece, 448 bp away.
fx pc_gap_skip_continue synthetic/pc_gap_skip_continue.fa \
   'chr_pc_gap_skip_continue=F:CCCTAAx100+L:200+F:CCCTAAx8+L:200+F:CCCTAAx100+L:5000' '-i' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;its=0;telolen=1200' \
   'A 48 bp array 200 bp behind the 600 bp array is too small for its gap and is skipped, but the 600 bp array 448 bp from the last joined piece outweighs that gap and joins: one fragmented p arm of teloLen 1200, hull 0-1648, with the 48 bp array inside the hull and no interstitial row'

# Ten 180 bp arrays 420 bp apart: each is shorter than the gap before it, and from the second on the gap exceeds -d.
fx pc_gap_train synthetic/pc_gap_train.fa \
   'chr_pc_gap_train=F:CCCTAAx100+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:5000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=10;telolen=600' \
   'A 600 bp array followed by a train of ten 180 bp arrays 420 bp apart: no member outweighs its gap, the arm is the 600 bp array alone and not fragmented, and the ten arrays are interstitial rows'

fx pc_gap_train_tip synthetic/pc_gap_train_tip.fa \
   'chr_pc_gap_train_tip=F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:420+F:CCCTAAx30+L:5000' '-i' \
   'type=none;anom=.;telo=0;labels=none;gaps=0;its=10' \
   'The same train of ten 180 bp arrays starting at the tip: the first is a piece but nothing outweighs its gap and 180 bp is short of -l, so there is no telomere and the ten arrays are interstitial rows'

# At -y 0.25 a gap weighs 1/3: 300 >= 600/3, while at the default -y 0.5 it weighs 1 and 300 < 600.
fx pc_gap_y_025 synthetic/pc_gap_y.fa \
   'chr_pc_gap_y=F:CCCTAAx100+L:600+F:CCCTAAx50+L:5000' '-i -y 0.25' \
   'type=incomplete;anom=fragmented_p;telo=1;labels=p;gaps=0;its=0;telolen=900' \
   'At -y 0.25 the 600 bp gap weighs 200 and the 300 bp array behind it outweighs it: one fragmented p arm of teloLen 900'

fx pc_gap_y_default synthetic/pc_gap_y.fa \
   'chr_pc_gap_y=F:CCCTAAx100+L:600+F:CCCTAAx50+L:5000' '-i' \
   'type=incomplete;anom=.;telo=1;labels=p;gaps=0;its=1;telolen=600' \
   'At -y 0.5 the 600 bp gap weighs 600 and the 300 bp array behind it does not outweigh it: the p arm is the 600 bp array alone and the 300 bp array is an interstitial row'

# The middle array is 300 bp from the q piece and as long as that gap, so it joins even though it lies past the -t zone.
fx pc_gap_extend_past_t synthetic/pc_gap_extend_past_t.fa \
   'chr_pc_gap_extend_past_t=F:CCCTAAx100+L:2700+R:TTAGGGx50+L:300+R:TTAGGGx100' '-i -t 300' \
   'type=t2t;anom=fragmented_q;telo=2;labels=pq;gaps=0;its=0;telolen=900' \
   'A 300 bp array as long as the 300 bp gap before the q piece joins it past the -t zone: the q arm is fragmented (teloLen 900) and nothing is interstitial'
