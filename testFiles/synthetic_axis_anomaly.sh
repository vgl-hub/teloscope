# shellcheck shell=bash

# Axis B: anomaly cross-product on fixed geometry.

# Anomalies emit in fixed order: discordant_p, discordant_q, balanced_p, balanced_q, misassembly.

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

fx ax_t2t_bal_p synthetic/ax_t2t_bal_p.fa \
   'chr_ax_t2t_bal_p=M:CCCTAATTAGGGx50+L:2800+R:TTAGGGx100' '-' \
   'type=t2t;anom=balanced_p;gran=P~Q;telo=2;labels=pq;gaps=0' \
   'Interleaved orientations at the p end give a balanced strand label'

fx ax_t2t_bal_q synthetic/ax_t2t_bal_q.fa \
   'chr_ax_t2t_bal_q=F:CCCTAAx100+L:2800+M:CCCTAATTAGGGx50' '-' \
   'type=t2t;anom=balanced_q;gran=PQ~;telo=2;labels=pq;gaps=0' \
   'Interleaved orientations at the q end'

fx ax_t2t_bal_pq synthetic/ax_t2t_bal_pq.fa \
   'chr_ax_t2t_bal_pq=M:CCCTAATTAGGGx50+L:2800+M:CCCTAATTAGGGx50' '-' \
   'type=t2t;anom=balanced_p,balanced_q;gran=P~Q~;telo=2;labels=pq;gaps=0' \
   'Both arms balanced, as an ALT or hairpin telomere would read'

fx ax_t2t_disc_p_bal_q synthetic/ax_t2t_disc_p_bal_q.fa \
   'chr_ax_t2t_disc_p_bal_q=R:TTAGGGx100+L:2800+M:CCCTAATTAGGGx50' '-' \
   'type=t2t;anom=discordant_p,balanced_q;gran=P*Q~;telo=2;labels=pq;gaps=0' \
   'One inverted arm and one balanced arm: both flags, in emission order'

fx ax_t2t_bal_p_disc_q synthetic/ax_t2t_bal_p_disc_q.fa \
   'chr_ax_t2t_bal_p_disc_q=M:CCCTAATTAGGGx50+L:2800+F:CCCTAAx100' '-' \
   'type=t2t;anom=discordant_q,balanced_p;gran=P~Q*;telo=2;labels=pq;gaps=0' \
   'The mirror: discordant_q sorts before balanced_p in the emitted order'

fx ax_t2t_extra_p synthetic/ax_t2t_extra_p.fa \
   'chr_ax_t2t_extra_p=F:CCCTAAx100+L:600+F:CCCTAAx100+L:2200+R:TTAGGGx100' '-' \
   'type=t2t;anom=misassembly;gran=PpQ;telo=2;labels=pq;gaps=0' \
   'Two arrays at the p end: the inner one is an extra terminal block'

fx ax_t2t_extra_q synthetic/ax_t2t_extra_q.fa \
   'chr_ax_t2t_extra_q=F:CCCTAAx100+L:2200+R:TTAGGGx100+L:600+R:TTAGGGx100' '-' \
   'type=t2t;anom=misassembly;gran=PqQ;telo=2;labels=pq;gaps=0' \
   'Two arrays at the q end: the >= tie-break elects the outermost'

fx ax_t2t_extra_p_disc_q synthetic/ax_t2t_extra_p_disc_q.fa \
   'chr_ax_t2t_extra_p_disc_q=F:CCCTAAx100+L:600+F:CCCTAAx100+L:2200+F:CCCTAAx100' '-' \
   'type=t2t;anom=discordant_q,misassembly;gran=PpQ*;telo=2;labels=pq;gaps=0' \
   'An extra p block and an inverted q arm together'

fx ax_incomplete_disc_p synthetic/ax_incomplete_disc_p.fa \
   'chr_ax_incomplete_disc_p=R:TTAGGGx100+L:2800' '-' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=0' \
   'One inverted arm and nothing at the other end'

fx ax_incomplete_bal_q synthetic/ax_incomplete_bal_q.fa \
   'chr_ax_incomplete_bal_q=L:2800+M:CCCTAATTAGGGx50' '-' \
   'type=incomplete;anom=balanced_q;gran=Q~;telo=1;labels=q;gaps=0' \
   'One balanced q arm and nothing at the other end'

fx ax_t2t_disc_p_gapped synthetic/ax_t2t_disc_p_gapped.fa \
   'chr_ax_t2t_disc_p_gapped=R:TTAGGGx100+L:1300+N:200+L:1300+R:TTAGGGx100' '-' \
   'type=t2t;anom=discordant_p;gran=P*Q;telo=2;labels=pq;gaps=1' \
   'An inverted p arm on a gapped scaffold: gappedness changes neither column'

fx ax_t2t_bal_pq_gapped synthetic/ax_t2t_bal_pq_gapped.fa \
   'chr_ax_t2t_bal_pq_gapped=M:CCCTAATTAGGGx50+L:1300+N:200+L:1300+M:CCCTAATTAGGGx50' '-' \
   'type=t2t;anom=balanced_p,balanced_q;gran=P~Q~;telo=2;labels=pq;gaps=1' \
   'Both arms balanced on a gapped scaffold'
