# shellcheck shell=bash

# Combinatorial tests: incomplete, seed boundary, window invariance.

# ---- incomplete x anomaly x gappedness ----
fx mo_inc_p_clean synthetic/mo_inc_p_clean.fa \
   'chr_mo_inc_p_clean=F:CCCTAAx100+L:2800' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'One concordant p arm'

fx mo_inc_p_clean_gapped synthetic/mo_inc_p_clean_gapped.fa \
   'chr_mo_inc_p_clean_gapped=F:CCCTAAx100+L:1300+N:200+L:1300' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1' \
   'One concordant p arm on a gapped scaffold'

fx mo_inc_q_clean synthetic/mo_inc_q_clean.fa \
   'chr_mo_inc_q_clean=L:2800+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'One concordant q arm'

fx mo_inc_q_disc synthetic/mo_inc_q_disc.fa \
   'chr_mo_inc_q_disc=L:2800+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=0' \
   'One inverted q arm'

fx mo_inc_p_bal_gapped synthetic/mo_inc_p_bal_gapped.fa \
   'chr_mo_inc_p_bal_gapped=M:CCCTAATTAGGGx50+L:1300+N:200+L:1300' '-' \
   'type=incomplete;anom=balanced_p;gran=P~;telo=1;labels=p;gaps=1' \
   'One balanced p arm on a gapped scaffold'

fx mo_inc_q_disc_gapped synthetic/mo_inc_q_disc_gapped.fa \
   'chr_mo_inc_q_disc_gapped=L:1300+N:200+L:1300+F:CCCTAAx100' '-' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=1' \
   'One inverted q arm on a gapped scaffold'

fx mo_inc_p_extra synthetic/mo_inc_p_extra.fa \
   'chr_mo_inc_p_extra=F:CCCTAAx100+L:600+F:CCCTAAx100+L:2200' '-' \
   'type=incomplete;anom=misassembly;gran=Pp;telo=1;labels=p;gaps=0' \
   'Two p arrays and no q arm'

fx mo_inc_q_extra synthetic/mo_inc_q_extra.fa \
   'chr_mo_inc_q_extra=L:2200+R:TTAGGGx100+L:600+R:TTAGGGx100' '-' \
   'type=incomplete;anom=misassembly;gran=qQ;telo=1;labels=q;gaps=0' \
   'Two q arrays and no p arm'

fx mo_inc_p_extra_disc synthetic/mo_inc_p_extra_disc.fa \
   'chr_mo_inc_p_extra_disc=R:TTAGGGx100+L:600+R:TTAGGGx100+L:2200' '-' \
   'type=incomplete;anom=discordant_p,misassembly;gran=P*p*;telo=1;labels=p;gaps=0' \
   'Two inverted p arrays: the elected one is flagged and the other is extra'

fx mo_inc_p_extra_bal synthetic/mo_inc_p_extra_bal.fa \
   'chr_mo_inc_p_extra_bal=M:CCCTAATTAGGGx50+L:600+F:CCCTAAx100+L:2200' '-' \
   'type=incomplete;anom=balanced_p,misassembly;gran=P~p;telo=1;labels=p;gaps=0' \
   'A balanced array and a plain forward array both at the p end, tied on canonical count'

# Docs elect the longest block (expect Pp); the binary elects by canonical count (pP), REG-005.
fx mo_elect_longer_loses synthetic/mo_elect_longer_loses.fa \
   'chr_mo_elect_longer_loses=M:CCCTAAACGATCx60+L:600+F:CCCTAAx100+L:2200' '-' \
   'type=incomplete;anom=misassembly;gran=Pp;telo=1;labels=p;gaps=0' \
   'The longer array loses the arm election to a shorter, more canonical one: REG-005'

fx mo_none_gapped_two synthetic/mo_none_gapped_two.fa \
   'chr_mo_none_gapped_two=L:900+N:150+L:900+N:150+L:900' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=2' \
   'No telomeres and two gaps'

fx mo_t2t_gapped_three synthetic/mo_t2t_gapped_three.fa \
   'chr_mo_t2t_gapped_three=F:CCCTAAx100+L:600+N:100+L:600+N:100+L:600+N:100+L:600+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=3' \
   'A t2t scaffold broken into four contigs'

fx mo_seed_gap_exact synthetic/mo_seed_gap_exact.fa \
   'chr_mo_seed_gap_exact=F:CCCTAAx50+L:50+F:CCCTAAx50+L:2400' '-d 1' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A 50 bp run between matches is exactly -k, so they chain into one seed'

fx mo_seed_gap_beyond synthetic/mo_seed_gap_beyond.fa \
   'chr_mo_seed_gap_beyond=F:CCCTAAx50+L:56+F:CCCTAAx50+L:2400' '-d 1' \
   'type=incomplete;anom=misassembly;gran=Pp;telo=1;labels=p;gaps=0' \
   'A 56 bp run exceeds -k, so the matches form two seeds that -d 1 cannot rejoin'

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
   'rec_a=F:CCCTAAx100+L:2800+R:TTAGGGx100;rec_b=L:3000;rec_c=L:2800+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0' \
   'Three records in one order'

fx mo_order_reversed synthetic/mo_order_reversed.fa \
   'rec_c=L:2800+R:TTAGGGx100;rec_b=L:3000;rec_a=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'The same three records in the opposite order get the same per-record answers'

fx mo_sparse_long synthetic/mo_sparse_long.fa \
   'chr_mo_sparse_long=M:CCCTAAACGATCx100+L:2400' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A long array at exactly -y density survives on both length and density'

fx mo_dense_short synthetic/mo_dense_short.fa \
   'chr_mo_dense_short=F:CCCTAAx40+L:2760' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'A fully dense 240 bp array still fails -l: density does not buy length'
