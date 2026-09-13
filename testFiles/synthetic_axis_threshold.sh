# shellcheck shell=bash

# Axis C: boundary tests at threshold limits.

# ---- -l / --min-block-length, default 300 ----
fx th_len_below synthetic/th_len_below.fa \
   'chr_th_len_below=F:CCCTAAx49+L:2706' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   '294 bp array is below -l 300 and is dropped'

fx th_len_exact synthetic/th_len_exact.fa \
   'chr_th_len_exact=F:CCCTAAx50+L:6700' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   '300 bp array is exactly -l and is kept: the comparison rejects on strict <'

fx th_len_above synthetic/th_len_above.fa \
   'chr_th_len_above=F:CCCTAAx51+L:6694' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   '306 bp array is above -l and is kept'

# ---- --terminal-tolerance, default 3000, in called bases ----
fx th_tol_exact synthetic/th_tol_exact.fa \
   'chr_th_tol_exact=L:3000+F:CCCTAAx100+L:6000' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Array starting exactly at --terminal-tolerance is still terminal; far filler keeps it out of the q end'

fx th_tol_beyond synthetic/th_tol_beyond.fa \
   'chr_th_tol_beyond=L:3006+F:CCCTAAx100+L:6000' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'Array starting 6 bp past --terminal-tolerance at the near end, and far from the other'

# docs/parameters.md:91: tolerance counts called bases; leading Ns don't count against it.
fx th_tol_called_bases synthetic/th_tol_called_bases.fa \
   'chr_th_tol_called_bases=N:1000+L:1000+F:CCCTAAx100+L:2800' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1' \
   'A leading 1000 bp gap does not push a telomere out of the terminal zone'

# Tolerance 3000 (up from 2000) now reaches an array starting 2500 bp in.
fx th_tol_2500 synthetic/th_tol_2500.fa \
   'chr_th_tol_2500=L:2500+F:CCCTAAx100+L:2400' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'An array starting 2500 bp in is within --terminal-tolerance 3000 and is called'

# ---- -d / --max-block-distance, default 500 (bridges called sequence within a piece) ----
fx th_bridge_plain_exact synthetic/th_bridge_plain_exact.fa \
   'chr_th_bridge_plain_exact=F:CCCTAAx50+L:500+F:CCCTAAx50+L:6100' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A 500 bp plain-sequence run is exactly -d: the density trim bridges it into one piece'

fx th_bridge_plain_beyond synthetic/th_bridge_plain_beyond.fa \
   'chr_th_bridge_plain_beyond=F:CCCTAAx50+L:506+F:CCCTAAx50+L:6094' '-' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;telolen=600' \
   'A 506 bp run exceeds -d for the density trim, but the second array starts within --link-distance and chains as a second piece'

# N runs are hard contig boundaries: pieces never chain across one, at any length (R1/R2).
fx th_bridge_gap synthetic/th_bridge_gap.fa \
   'chr_th_bridge_gap=F:CCCTAAx50+N:500+F:CCCTAAx50+L:3600' '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1;its=1' \
   'A 500 bp N run splits the record into two contigs; the p arm is only the first piece and the second array is an interstitial row'

# ---- --min-its-length is retired: interstitial rows have no length floor (R5) ----
fx th_its_below synthetic/th_its_below.fa \
   'chr_th_its_below=F:CCCTAAx100+L:2400+R:TTAGGGx16+L:2400+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'A 96 bp interstitial array is reported: interstitial rows have no length floor any more'

fx th_its_exact synthetic/th_its_exact.fa \
   'chr_th_its_exact=F:CCCTAAx100+L:2400+R:TTAGGGx17+L:2400+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'A 102 bp interstitial array is reported too'

# minCanonicalCount threshold (REG-003); --min-its-length no longer exists.
fx th_its_canon_below synthetic/th_its_canon_below.fa \
   'chr_th_its_canon_below=F:CCCTAAx100+L:2400+R:TTAGGGx3+L:2400+R:TTAGGGx100' \
   '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=0' \
   'Three canonical matches are below the hardcoded minCanonicalCount of 4'

fx th_its_canon_exact synthetic/th_its_canon_exact.fa \
   'chr_th_its_canon_exact=F:CCCTAAx100+L:2400+R:TTAGGGx4+L:2400+R:TTAGGGx100' \
   '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'Four canonical matches are exactly minCanonicalCount and are reported'

fx th_its_canon_above synthetic/th_its_canon_above.fa \
   'chr_th_its_canon_above=F:CCCTAAx100+L:2400+V:TTAGGAx2+R:TTAGGGx1+V:TTAGGAx2+R:TTAGGGx1+V:TTAGGAx2+R:TTAGGGx1+V:TTAGGAx2+R:TTAGGGx1+V:TTAGGAx2+R:TTAGGGx1+V:TTAGGAx2+L:2400+R:TTAGGGx100' \
   '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'Five non-adjacent canonical repeats among variants, 102 bp: clears minCanonicalCount 4'

# ---- --label-threshold, default 0.667: terminal pieces are strand-pure (R1), so this is
# only observable on interstitial rows; the label itself has no report column, so these
# fixtures assert its=1 and note the label the rule predicts. ----
fx th_strand_2fwd_1rev synthetic/th_strand_2fwd_1rev.fa \
   'chr_th_strand_2fwd_1rev=L:3500+M:CCCTAACCCTAATTAGGGx40+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   'Two forward per one reverse is 666.7 per mille, just short of the 667 boundary: label b'

fx th_strand_1fwd_2rev synthetic/th_strand_1fwd_2rev.fa \
   'chr_th_strand_1fwd_2rev=L:3500+M:CCCTAATTAGGGTTAGGGx40+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   'One forward per two reverse is 33.3 per mille below 1000-667=333: label b, not q'

# ---- strand sweep at 0, 333, 667, 1000 per mille: exact fractions ----
fx th_strand_0 synthetic/th_strand_0.fa \
   'chr_th_strand_0=L:3500+M:TTAGGGx100+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   '100 reverse matches, 0 forward: label q, the low extreme'

fx th_strand_333 synthetic/th_strand_333.fa \
   'chr_th_strand_333=L:3500+M:TTAGGGTTAGGGCCCTAAx333+M:TTAGGGx1+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   '333 forward matches of 1000 total: exactly at the lower boundary (1000-667), label b'

fx th_strand_666 synthetic/th_strand_666.fa \
   'chr_th_strand_666=L:3500+M:CCCTAACCCTAATTAGGGx166+M:CCCTAATTAGGGx1+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   '333 forward matches of 500 total: 666 per mille, just short of the round(0.667*1000)=667 boundary, label b'

fx th_strand_1000 synthetic/th_strand_1000.fa \
   'chr_th_strand_1000=L:3500+M:CCCTAAx100+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   '100 forward matches, 0 reverse: label p, the high extreme'

# ---- -y / --min-block-density, default 0.5, on the interstitial (all-match) gate too ----
fx th_density_exact synthetic/th_density_exact.fa \
   'chr_th_density_exact=L:3500+M:CCCTAAACGATCx50+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   '6 bp match against 6 bp filler is exactly -y 0.5 for the interstitial gate: kept on the tie'

fx th_density_below synthetic/th_density_below.fa \
   'chr_th_density_below=L:3500+M:CCCTAAACGATCGACTx50+L:3500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=0' \
   '6 bp match against 10 bp filler is 37.5%, below -y 0.5: dropped'

# ---- --min-block-counts, default 2 ----
fx th_counts_below synthetic/th_counts_below.fa \
   'chr_th_counts_below=F:CCCTAAx1+L:2994' '-l 6' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'A single match is below --min-block-counts 2'

fx th_counts_exact synthetic/th_counts_exact.fa \
   'chr_th_counts_exact=F:CCCTAAx2+L:6988' '-l 6' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'Two matches are exactly --min-block-counts and are kept'
