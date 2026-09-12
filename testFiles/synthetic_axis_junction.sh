# shellcheck shell=bash

# Axis F -- junction geometry. Sourced by generate_synthetic.sh.

# Axis F: opposite matches inside a p array (REG-014).
fx jn_p_then_q_k1 synthetic/jn_p_then_q_k1.fa \
   'chr_jn_p_then_q_k1=F:CCCTAAx100+R:TTAGGGx1+F:CCCTAAx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'One opposite match inside the p array is absorbed: forward share 200/201'

fx jn_p_then_q_k2 synthetic/jn_p_then_q_k2.fa \
   'chr_jn_p_then_q_k2=F:CCCTAAx100+R:TTAGGGx2+F:CCCTAAx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Two opposite matches inside the p array are absorbed: forward share 200/202'

fx jn_p_then_q_k3 synthetic/jn_p_then_q_k3.fa \
   'chr_jn_p_then_q_k3=F:CCCTAAx100+R:TTAGGGx3+F:CCCTAAx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=?;gran=?;telo=2;labels=pq;gaps=0' \
   'Three or more opposite matches cut the seed; how the split reads is decided (D14) and lands with B3b, REG-014'

fx jn_p_then_q_k4 synthetic/jn_p_then_q_k4.fa \
   'chr_jn_p_then_q_k4=F:CCCTAAx100+R:TTAGGGx4+F:CCCTAAx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=?;gran=?;telo=2;labels=pq;gaps=0' \
   'Three or more opposite matches cut the seed; how the split reads is decided (D14) and lands with B3b, REG-014'

fx jn_lone_opposite synthetic/jn_lone_opposite.fa \
   'chr_jn_lone_opposite=F:CCCTAAx100+L:60+R:TTAGGGx1+L:2400' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A single reverse match past the p array is below --min-block-counts and vanishes'

# Abutting opposite arrays on short contig (REG-017).
fx jn_short_q_then_p synthetic/jn_short_q_then_p.fa \
   'chr_jn_short_q_then_p=L:1500+R:TTAGGGx60+F:CCCTAAx60+L:1500' '-i' \
   'type=t2t;anom=discordant_p,discordant_q;gran=P*Q*;telo=2;labels=pq;gaps=0;its=0' \
   'Two abutting opposite arrays, each within tolerance of its end: the docs make both terminal arms; the binary drops the second, REG-017'

fx jn_short_p_then_q synthetic/jn_short_p_then_q.fa \
   'chr_jn_short_p_then_q=L:1500+F:CCCTAAx60+R:TTAGGGx60+L:1500' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=0' \
   'Same geometry, forward first: p arm then a concordant q arm per the docs; the binary drops the second, REG-017'

fx jn_short_p_then_p_gap synthetic/jn_short_p_then_p_gap.fa \
   'chr_jn_short_p_then_p_gap=L:1500+F:CCCTAAx60+N:100+F:CCCTAAx60+L:1500' '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1;its=0' \
   'A 100 bp N run bridges (<= -d) into one 820 bp gapped p arm; no ITS is left over'

fx jn_short_p_then_p_far synthetic/jn_short_p_then_p_far.fa \
   'chr_jn_short_p_then_p_far=L:1500+F:CCCTAAx60+L:600+F:CCCTAAx60+L:1500' '-i' \
   'type=t2t;anom=discordant_q;gran=PQ*;telo=2;labels=pq;gaps=0;its=0' \
   'A 600 bp run past -d leaves two arms, both within tolerance of their own end: p, and a discordant q'

# Interior junctions past terminal tolerance (REG-004).
fx jn_interior_q_then_p synthetic/jn_interior_q_then_p.fa \
   'chr_jn_interior_q_then_p=L:2500+R:TTAGGGx60+F:CCCTAAx60+L:2500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=?' \
   'Interior TTAGGG array then CCCTAA array, no gap: the fusion signature'

fx jn_interior_p_then_q synthetic/jn_interior_p_then_q.fa \
   'chr_jn_interior_p_then_q=L:2500+F:CCCTAAx60+R:TTAGGGx60+L:2500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=?' \
   'Interior CCCTAA array then TTAGGG array, no gap: tail to tail'

fx jn_interior_p_then_p_gap synthetic/jn_interior_p_then_p_gap.fa \
   'chr_jn_interior_p_then_p_gap=L:2500+F:CCCTAAx60+N:100+F:CCCTAAx60+L:2500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=1' \
   'Interior array split by a bridgeable N run stays one interstitial block'

fx jn_interior_p_then_p_far synthetic/jn_interior_p_then_p_far.fa \
   'chr_jn_interior_p_then_p_far=L:2500+F:CCCTAAx60+L:600+F:CCCTAAx60+L:2500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'Two interior same-orientation arrays further apart than -d are two blocks'

fx jn_fragmented_p synthetic/jn_fragmented_p.fa \
   'chr_jn_fragmented_p=F:CCCTAAx317+L:598+F:CCCTAAx2917+L:3000+R:TTAGGGx100' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'A 598 bp run past -d separates the p arm from a second forward piece that reads as one ITS'

fx jn_fragmented_p_plain synthetic/jn_fragmented_p.fa \
   'chr_jn_fragmented_p=F:CCCTAAx317+L:598+F:CCCTAAx2917+L:3000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Same file without -i: the fragment is simply not scanned, its is not checked'

fx jn_its_96bp synthetic/jn_its_96bp.fa \
   'chr_jn_its_96bp=L:2500+R:TTAGGGx16+L:2500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=0' \
   '96 bp is below --min-its-length 100: no ITS is reported (REG-003/D16 changes this later)'

fx jn_gap_inside_telomere synthetic/jn_gap_inside_telomere.fa \
   'chr_jn_gap_inside_telomere=F:CCCTAAx50+N:100+F:CCCTAAx50+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=1' \
   'A 100 bp N run inside the p array is bridged (<= -d) into one gapped block'

# Internal contig-end telomere (REG-015).
fx jn_contig_end_telomere synthetic/jn_contig_end_telomere.fa \
   'chr_jn_contig_end_telomere=F:CCCTAAx100+L:2000+N:100+L:1500+R:TTAGGGx100+N:100+L:2000+R:TTAGGGx100' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=2;its=1' \
   'An internal contig-end telomere reads as one interstitial block under current docs, REG-015'

fx jn_contig_end_telomere_n synthetic/jn_contig_end_telomere.fa \
   'chr_jn_contig_end_telomere=F:CCCTAAx100+L:2000+N:100+L:1500+R:TTAGGGx100+N:100+L:2000+R:TTAGGGx100' '-n -i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=2;its=1' \
   'Same file with -n: docs/parameters.md:103 says -n changes nothing yet, REG-015'
