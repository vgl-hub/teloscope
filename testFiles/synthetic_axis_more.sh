# shellcheck shell=bash

# Combinatorial tests: incomplete, seed boundary, window invariance.

# ---- incomplete x anomaly x gappedness ----
fx mo_inc_p_clean synthetic/mo_inc_p_clean.fa \
   'chr_mo_inc_p_clean=F:CCCTAAx100+L:6600' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'One concordant p arm'

fx mo_inc_p_clean_gapped synthetic/mo_inc_p_clean_gapped.fa \
   'chr_mo_inc_p_clean_gapped=F:CCCTAAx100+L:1300+N:200+L:1300' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1' \
   'One concordant p arm on a gapped scaffold'

fx mo_inc_q_clean synthetic/mo_inc_q_clean.fa \
   'chr_mo_inc_q_clean=L:6600+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'One concordant q arm'

# Filler > tolerance 3000, so the p end never reaches this array too and reclaims it (R3).
fx mo_inc_q_disc synthetic/mo_inc_q_disc.fa \
   'chr_mo_inc_q_disc=L:7200+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=0' \
   'One inverted q arm'

# Interleaved orientations never qualify as a chain piece (R1/R3): no telomere at either end.
fx mo_inc_p_bal_gapped synthetic/mo_inc_p_bal_gapped.fa \
   'chr_mo_inc_p_bal_gapped=M:CCCTAATTAGGGx50+L:1300+N:200+L:1300' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1' \
   'An interleaved array on the first contig and a pure-filler last contig: no telomere at either end'

fx mo_inc_q_disc_gapped synthetic/mo_inc_q_disc_gapped.fa \
   'chr_mo_inc_q_disc_gapped=L:1300+N:200+L:1300+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=1' \
   'One inverted q arm on a gapped scaffold'

# Two same-strand arrays within --link-distance chain into one fragmented arm (R2); the
# canonical-count election they used to test (REG-005) is retired.
fx mo_inc_p_extra synthetic/mo_inc_p_extra.fa \
   'chr_mo_inc_p_extra=F:CCCTAAx100+L:600+F:CCCTAAx100+L:5400' '-' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0' \
   'Two p arrays 600 bp apart chain into one fragmented arm; no q arm'

fx mo_inc_q_extra synthetic/mo_inc_q_extra.fa \
   'chr_mo_inc_q_extra=L:5400+R:TTAGGGx100+L:600+R:TTAGGGx100' '-' \
   'type=incomplete;anom=fragmented_q;gran=Q;telo=1;labels=q;gaps=0' \
   'Two q arrays 600 bp apart chain into one fragmented arm; no p arm'

# Three pieces 700 bp apart: past -d 500 for density-trim, each still within --link-distance.
fx mo_inc_p_three_piece synthetic/mo_inc_p_three_piece.fa \
   'chr_mo_inc_p_three_piece=F:CCCTAAx100+L:700+F:CCCTAAx100+L:700+F:CCCTAAx100+L:4000' '-' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;telolen=1800' \
   'Three p arrays, each 700 bp from the next, chain into one fragmented arm; teloLen is the sum of the three pieces'

# Filler > tolerance 3000 beyond the second piece, so the q end never reaches this
# structure too and reclaims it (R3).
fx mo_inc_p_extra_disc synthetic/mo_inc_p_extra_disc.fa \
   'chr_mo_inc_p_extra_disc=R:TTAGGGx100+L:600+R:TTAGGGx100+L:6600' '-' \
   'type=incomplete;anom=discordant_p,fragmented_p;gran=P*;telo=1;labels=p;gaps=0' \
   'Two inverted p arrays chain into one fragmented, discordant p arm'

fx mo_inc_p_extra_bal synthetic/mo_inc_p_extra_bal.fa \
   'chr_mo_inc_p_extra_bal=M:CCCTAATTAGGGx50+L:600+F:CCCTAAx100+L:5400' '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0;its=1' \
   'Interleaved orientations at the very start fail every anchor attempt (R3) and are skipped; the plain forward array 600 bp in is the p arm, and the interleaved run is a single interstitial b row'

# Chaining by R2 replaces canonical-count election outright (REG-005 retired).
fx mo_elect_longer_loses synthetic/mo_elect_longer_loses.fa \
   'chr_mo_elect_longer_loses=M:CCCTAAACGATCx60+L:600+F:CCCTAAx100+L:5280' '-' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;telolen=1314' \
   'A dense alternating array and a plain array 600 bp apart chain into one fragmented arm; the first piece ends at its last match (714 bp), not the trailing filler unit'

fx mo_none_gapped_two synthetic/mo_none_gapped_two.fa \
   'chr_mo_none_gapped_two=L:900+N:150+L:900+N:150+L:900' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=2' \
   'No telomeres and two gaps'

fx mo_t2t_gapped_three synthetic/mo_t2t_gapped_three.fa \
   'chr_mo_t2t_gapped_three=F:CCCTAAx100+L:600+N:100+L:600+N:100+L:600+N:100+L:600+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=3' \
   'A t2t scaffold broken into four contigs'

# Interior seed grouping (R5): two same-orientation arrays group into one interstitial span
# when the gap between them is within -d, and split apart again just beyond it. -k (seed
# formation from raw matches) does not affect this at all -- both arrays are single seeds
# either way, so it never comes into play here.
fx mo_seed_gap_exact synthetic/mo_seed_gap_exact.fa \
   'chr_mo_seed_gap_exact=L:3500+F:CCCTAAx100+L:500+F:CCCTAAx100+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   'A 500 bp run is exactly -d: the two interior arrays group into one interstitial row'

fx mo_seed_gap_beyond synthetic/mo_seed_gap_beyond.fa \
   'chr_mo_seed_gap_beyond=L:3500+F:CCCTAAx100+L:501+F:CCCTAAx100+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'A 501 bp run exceeds -d: the two interior arrays stay two interstitial rows, classed fragmentation'

fx mo_window_overlapping synthetic/mo_window_overlapping.fa \
   'chr_mo_window_overlapping=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-r -w 100 -s 50' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Overlapping windows change the tracks, never the call'

fx mo_window_large synthetic/mo_window_large.fa \
   'chr_mo_window_large=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-r -w 2000 -s 2000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'A window wider than either arm changes the tracks, never the call'

fx mo_long_scaffold_default synthetic/mo_long_scaffold.fa \
   'chr_mo_long_scaffold=F:CCCTAAx100+L:20000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'A 21 kb scaffold is well inside the default terminal limit'

fx mo_long_scaffold_narrow synthetic/mo_long_scaffold_narrow.fa \
   'chr_mo_long_scaffold_narrow=L:3000+F:CCCTAAx100+L:20000+R:TTAGGGx100+L:3000' '-t 2000' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'Both arrays sit 3 kb in, outside a 2 kb terminal zone'

fx mo_order_forward synthetic/mo_order_forward.fa \
   'rec_a=F:CCCTAAx100+L:2800+R:TTAGGGx100;rec_b=L:3000;rec_c=L:6600+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'Three records in one order'

fx mo_order_reversed synthetic/mo_order_reversed.fa \
   'rec_c=L:6600+R:TTAGGGx100;rec_b=L:3000;rec_a=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'The same three records in the opposite order get the same per-record answers'

fx mo_sparse_long synthetic/mo_sparse_long.fa \
   'chr_mo_sparse_long=M:CCCTAAACGATCx100+L:6000' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A long array at exactly -y density survives on both length and density'

fx mo_dense_short synthetic/mo_dense_short.fa \
   'chr_mo_dense_short=F:CCCTAAx40+L:2760' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'A fully dense 240 bp array still fails -l: density does not buy length'
