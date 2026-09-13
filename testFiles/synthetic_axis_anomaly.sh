# shellcheck shell=bash

# Axis B: anomaly cross-product on fixed geometry.

# Anomalies emit in fixed order: discordant_p, discordant_q, fragmented_p, fragmented_q.

fx ax_t2t_clean synthetic/ax_t2t_clean.fa \
   'chr_ax_t2t_clean=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Both arms concordant: the control for this axis'

fx ax_t2t_disc_p synthetic/ax_t2t_disc_p.fa \
   'chr_ax_t2t_disc_p=R:TTAGGGx100+L:2800+R:TTAGGGx100' '-' \
   'type=t2t;anom=discordant_p;gran=P*Q;telo=2;labels=pq;gaps=0' \
   'Reverse motif at the p end: the p arm points the wrong way'

fx ax_t2t_disc_q synthetic/ax_t2t_disc_q.fa \
   'chr_ax_t2t_disc_q=F:CCCTAAx100+L:2800+F:CCCTAAx100' '-' \
   'type=t2t;anom=discordant_q;gran=PQ*;telo=2;labels=pq;gaps=0' \
   'Forward motif at the q end: the q arm points the wrong way'

fx ax_t2t_disc_pq synthetic/ax_t2t_disc_pq.fa \
   'chr_ax_t2t_disc_pq=R:TTAGGGx100+L:2800+F:CCCTAAx100' '-' \
   'type=t2t;anom=discordant_p,discordant_q;gran=P*Q*;telo=2;labels=pq;gaps=0' \
   'Both arms inverted: an inverted terminal repeat at each end'

# Interleaved orientations never qualify as a chain piece (R1): they read as one
# interstitial 'b' row, never a terminal arm. Moved to the interior, checked under -i.
fx ax_t2t_bal_p synthetic/ax_t2t_bal_p.fa \
   'chr_ax_t2t_bal_p=L:3500+M:CCCTAATTAGGGx50+L:3500+R:TTAGGGx100' '-i' \
   'type=incomplete;anom=.;gran=Q;telo=1;labels=q;gaps=0;its=1' \
   'Interleaved orientations near the p end: an interior b row, not a p arm; --label-threshold'

fx ax_t2t_bal_q synthetic/ax_t2t_bal_q.fa \
   'chr_ax_t2t_bal_q=F:CCCTAAx100+L:3500+M:CCCTAATTAGGGx50+L:3500' '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0;its=1' \
   'Interleaved orientations near the q end: an interior b row, not a q arm'

fx ax_t2t_bal_pq synthetic/ax_t2t_bal_pq.fa \
   'chr_ax_t2t_bal_pq=L:3500+M:CCCTAATTAGGGx50+L:2000+M:CCCTAATTAGGGx50+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'Two interior interleaved arrays, neither end has a telomere'

fx ax_t2t_disc_p_bal_q synthetic/ax_t2t_disc_p_bal_q.fa \
   'chr_ax_t2t_disc_p_bal_q=R:TTAGGGx100+L:3500+M:CCCTAATTAGGGx50+L:3500' '-i' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=0;its=1' \
   'An inverted p arm and an interior interleaved array'

fx ax_t2t_bal_p_disc_q synthetic/ax_t2t_bal_p_disc_q.fa \
   'chr_ax_t2t_bal_p_disc_q=L:3500+M:CCCTAATTAGGGx50+L:3500+F:CCCTAAx100' '-i' \
   'type=incomplete;anom=discordant_q;gran=Q*;telo=1;labels=q;gaps=0;its=1' \
   'The mirror: an interior interleaved array and an inverted q arm'

# Same-strand pieces within --link-distance chain into one row (R2): the extra
# block becomes a second piece of the same arm, one row, one granular letter.
fx ax_t2t_extra_p synthetic/ax_t2t_extra_p.fa \
   'chr_ax_t2t_extra_p=F:CCCTAAx100+L:600+F:CCCTAAx100+L:2200+R:TTAGGGx100' '-' \
   'type=t2t;anom=fragmented_p;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Two p arrays 600 bp apart chain into one fragmented arm'

fx ax_t2t_extra_q synthetic/ax_t2t_extra_q.fa \
   'chr_ax_t2t_extra_q=F:CCCTAAx100+L:2200+R:TTAGGGx100+L:600+R:TTAGGGx100' '-' \
   'type=t2t;anom=fragmented_q;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Two q arrays 600 bp apart chain into one fragmented arm'

fx ax_t2t_extra_p_disc_q synthetic/ax_t2t_extra_p_disc_q.fa \
   'chr_ax_t2t_extra_p_disc_q=F:CCCTAAx100+L:600+F:CCCTAAx100+L:2200+F:CCCTAAx100' '-' \
   'type=t2t;anom=discordant_q,fragmented_p;gran=PQ*;telo=2;labels=pq;gaps=0' \
   'A fragmented p arm and an unrelated forward array at the q end'

# Filler > tolerance 3000, so the q end never reaches this array too and reclaims it (R3).
fx ax_incomplete_disc_p synthetic/ax_incomplete_disc_p.fa \
   'chr_ax_incomplete_disc_p=R:TTAGGGx100+L:7200' '-' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=0' \
   'One inverted arm and nothing at the other end'

fx ax_incomplete_bal_q synthetic/ax_incomplete_bal_q.fa \
   'chr_ax_incomplete_bal_q=L:3500+M:CCCTAATTAGGGx50+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   'An interior interleaved array and nothing else: no telomere at either end'

fx ax_t2t_disc_p_gapped synthetic/ax_t2t_disc_p_gapped.fa \
   'chr_ax_t2t_disc_p_gapped=R:TTAGGGx100+L:1300+N:200+L:1300+R:TTAGGGx100' '-' \
   'type=t2t;anom=discordant_p;gran=P*Q;telo=2;labels=pq;gaps=1' \
   'An inverted p arm on a gapped scaffold: gappedness changes neither column'

fx ax_t2t_bal_pq_gapped synthetic/ax_t2t_bal_pq_gapped.fa \
   'chr_ax_t2t_bal_pq_gapped=L:3500+M:CCCTAATTAGGGx50+L:1300+N:200+L:1300+M:CCCTAATTAGGGx50+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=2' \
   'Two interior interleaved arrays on a gapped scaffold, one per contig'
