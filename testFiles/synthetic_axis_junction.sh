# shellcheck shell=bash

# Axis F -- junction geometry. Sourced by generate_synthetic.sh.

# Axis F: opposite matches inside a p array. A lone opposite match (however many, as long
# as it never reaches -l on its own) is never "real" (R1) and is absorbed as bridgeable gap.
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
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Three opposite matches (18 bp) are still far short of a real reverse array (< -l) and are absorbed the same way (D14, REG-014 retired)'

fx jn_p_then_q_k4 synthetic/jn_p_then_q_k4.fa \
   'chr_jn_p_then_q_k4=F:CCCTAAx100+R:TTAGGGx4+F:CCCTAAx100+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Four opposite matches (24 bp) are still far short of a real reverse array and are absorbed the same way'

fx jn_lone_opposite synthetic/jn_lone_opposite.fa \
   'chr_jn_lone_opposite=F:CCCTAAx100+L:60+R:TTAGGGx1+L:2400' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A single reverse match past the p array is below --min-block-counts, is not real, and vanishes'

# Abutting opposite arrays on a short contig: each clamps the other's piece at the exact
# junction (R1), so both become real, independent terminal rows (REG-017 retired).
fx jn_short_q_then_p synthetic/jn_short_q_then_p.fa \
   'chr_jn_short_q_then_p=L:1500+R:TTAGGGx60+F:CCCTAAx60+L:1500' '-i' \
   'type=t2t;anom=discordant_p,discordant_q;gran=P*Q*;telo=2;labels=pq;gaps=0;its=0' \
   'Two abutting opposite arrays, each clamped at the junction and each within tolerance of its own end: both are real terminal rows'

fx jn_short_p_then_q synthetic/jn_short_p_then_q.fa \
   'chr_jn_short_p_then_q=L:1500+F:CCCTAAx60+R:TTAGGGx60+L:1500' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=0' \
   'Same geometry, forward first: a concordant p arm then a concordant q arm, both real terminal rows'

# N runs are hard contig boundaries (R1/R2): re-spaced so the second contig's array falls
# beyond --terminal-tolerance and stays a plain interstitial row.
fx jn_short_p_then_p_gap synthetic/jn_short_p_then_p_gap.fa \
   'chr_jn_short_p_then_p_gap=L:1500+F:CCCTAAx60+N:100+F:CCCTAAx60+L:3200' '-i' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1;its=1' \
   'A 100 bp N run splits the array into two contigs; pieces never chain across it, so the p arm is only the first piece and the second contig'"'"'s array is an interstitial row'

# Both arrays are forward and 600 bp apart; the q-side chain finds the identical two
# pieces the p-side chain does, so R3 gives the whole fragmented structure to the p end.
fx jn_short_p_then_p_far synthetic/jn_short_p_then_p_far.fa \
   'chr_jn_short_p_then_p_far=L:1500+F:CCCTAAx60+L:600+F:CCCTAAx60+L:1500' '-i' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;its=0' \
   'Two forward arrays 600 bp apart, both within tolerance of their own end: strand decides (R3), one fragmented p arm'

# Interior junctions, re-spaced to at least 4000 bp from both record ends so neither array
# is ever reachable by a terminal chain, keeping these pure interstitial junction tests.
fx jn_interior_q_then_p synthetic/jn_interior_q_then_p.fa \
   'chr_jn_interior_q_then_p=L:4000+R:TTAGGGx60+F:CCCTAAx60+L:4000' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'Interior TTAGGG array then CCCTAA array, abutting: the fusion signature, two plain interstitial rows classed fusion'

fx jn_interior_p_then_q synthetic/jn_interior_p_then_q.fa \
   'chr_jn_interior_p_then_q=L:4000+F:CCCTAAx60+R:TTAGGGx60+L:4000' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=2' \
   'Interior CCCTAA array then TTAGGG array, abutting: tail to tail, two plain interstitial rows'

fx jn_interior_p_then_p_gap synthetic/jn_interior_p_then_p_gap.fa \
   'chr_jn_interior_p_then_p_gap=L:4000+F:CCCTAAx60+N:100+F:CCCTAAx60+L:4000' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1;its=2' \
   'Interior array split by an N run: two separate single-contig interstitial rows, since spans never bridge a gap'

# Left at 2500 bp, not re-spaced: tolerance 3000 now reaches the first array, and the
# second chains within --link-distance, so this one flips into a fragmented p arm.
fx jn_interior_p_then_p_far synthetic/jn_interior_p_then_p_far.fa \
   'chr_jn_interior_p_then_p_far=L:2500+F:CCCTAAx60+L:600+F:CCCTAAx60+L:2500' '-i' \
   'type=incomplete;anom=fragmented_p;gran=P;telo=1;labels=p;gaps=0;its=0;telolen=720' \
   'Two forward arrays 2500 bp and 3460 bp in, 600 bp apart: both reachable at tolerance 3000, one fragmented p arm'

fx jn_fragmented_p synthetic/jn_fragmented_p.fa \
   'chr_jn_fragmented_p=F:CCCTAAx317+L:598+F:CCCTAAx2917+L:3000+R:TTAGGGx100' '-i' \
   'type=t2t;anom=fragmented_p;gran=PQ;telo=2;labels=pq;gaps=0;its=0;telolen=19404' \
   'A 598 bp run exceeds -d for density-trim but chains within --link-distance: one fragmented p arm, teloLen is the sum of the two pieces, no ITS left over'

fx jn_fragmented_p_plain synthetic/jn_fragmented_p.fa \
   'chr_jn_fragmented_p=F:CCCTAAx317+L:598+F:CCCTAAx2917+L:3000+R:TTAGGGx100' '-' \
   'type=t2t;anom=fragmented_p;gran=PQ;telo=2;labels=pq;gaps=0;telolen=19404' \
   'Same file without -i: fast mode covers this short record whole, identically to full scan; its is not checked'

fx jn_its_96bp synthetic/jn_its_96bp.fa \
   'chr_jn_its_96bp=L:2500+R:TTAGGGx16+L:2500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0;its=1' \
   '96 bp is within --terminal-tolerance but far too short for -l to anchor a piece; reported as a plain interstitial row now that there is no length floor for one'

fx jn_gap_inside_telomere synthetic/jn_gap_inside_telomere.fa \
   'chr_jn_gap_inside_telomere=F:CCCTAAx50+N:100+F:CCCTAAx50+L:2000+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=1;telolen=300' \
   'A 100 bp N run splits the p array into two contigs; pieces never chain across a gap, so the p arm is only the first piece (teloLen 300, not the old gap-bridged 700)'

# Internal contig-end telomere: an ordinary interstitial row without -n, a contig row in
# the terminal BED (excluded from the interstitial BED) with -n.
fx jn_contig_end_telomere synthetic/jn_contig_end_telomere.fa \
   'chr_jn_contig_end_telomere=F:CCCTAAx100+L:2000+N:100+L:1500+R:TTAGGGx100+N:100+L:2000+R:TTAGGGx100' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=2;its=1' \
   'An internal contig-end telomere reads as an ordinary interstitial row without -n'

fx jn_contig_end_telomere_n synthetic/jn_contig_end_telomere.fa \
   'chr_jn_contig_end_telomere=F:CCCTAAx100+L:2000+N:100+L:1500+R:TTAGGGx100+N:100+L:2000+R:TTAGGGx100' '-n -i' \
   'type=t2t;anom=.;gran=PqQ;telo=2;labels=pq;gaps=2;its=0' \
   'Same file with -n: the internal contig-end telomere moves to the terminal BED as a lowercase contig row and leaves the interstitial BED (its 1 -> 0)'

# A three-contig scaffold: true p and q arms on the outer contigs, an internal contig-end
# telomere on the middle one.
fx jn_contig_middle synthetic/jn_contig_middle.fa \
   'chr_jn_contig_middle=F:CCCTAAx100+L:2000+N:100+L:1500+R:TTAGGGx100+N:100+L:2000+R:TTAGGGx100' '-i' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=2;its=1' \
   'Without -n the middle contig'"'"'s telomere is a plain interstitial row; the scaffold arms are unaffected'

fx jn_contig_middle_n synthetic/jn_contig_middle.fa \
   'chr_jn_contig_middle=F:CCCTAAx100+L:2000+N:100+L:1500+R:TTAGGGx100+N:100+L:2000+R:TTAGGGx100' '-n -i' \
   'type=t2t;anom=.;gran=PqQ;telo=2;labels=pq;gaps=2;its=0' \
   'With -n the middle contig'"'"'s telomere becomes a lowercase contig row and leaves the interstitial BED'

# An array flanked by N on both sides forms a contig of its own; -n builds chains from both
# of that contig's own ends, which R3 collapses to one contig row (same strand, same span).
fx jn_flanked_both synthetic/jn_flanked_both.fa \
   'chr_jn_flanked_both=L:1500+N:100+F:CCCTAAx60+N:100+L:1500' '-i' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=2;its=1' \
   'Without -n the flanked array is a plain interstitial row; neither outer contig has a telomere'

fx jn_flanked_both_n synthetic/jn_flanked_both.fa \
   'chr_jn_flanked_both=L:1500+N:100+F:CCCTAAx60+N:100+L:1500' '-n -i' \
   'type=none;anom=.;gran=p;telo=0;labels=none;gaps=2;its=0' \
   'With -n both of the contig'"'"'s own internal chains find the same array; R3 keeps one lowercase contig row'

# A 300 bp inverted tip anchors the p end (R2: outermost, whatever its strand); the far
# larger forward telomere 1500 bp behind it is outside --link-distance and is never even
# tested as a chain-ending array, so it is left as a plain interstitial row.
fx jn_tip_beats_giant synthetic/jn_tip_beats_giant.fa \
   'chr_jn_tip_beats_giant=R:TTAGGGx50+L:1500+F:CCCTAAx500+L:3200' '-i' \
   'type=incomplete;anom=discordant_p;gran=P*;telo=1;labels=p;gaps=0;its=1;telolen=300' \
   'A 300 bp reverse tip is the outermost qualifying piece and becomes the (discordant) p arm; the 3 kb forward telomere behind it is too far to end the chain and is reported as an interstitial row'
