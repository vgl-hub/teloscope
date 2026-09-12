# shellcheck shell=bash

# Axis D/E: degenerate records, masking, motifs, and assembly shapes.

# ---- degenerate records ----
fx st_single_base synthetic/st_single_base.fa \
   'chr_st_single_base=S:A' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'A one-base record must classify, not crash'

fx st_sub_window synthetic/st_sub_window.fa \
   'chr_st_sub_window=F:CCCTAAx10' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=0' \
   'A 60 bp record is shorter than one window and below -l'

fx st_all_gap synthetic/st_all_gap.fa \
   'chr_st_all_gap=N:3000' '-' \
   'type=none;anom=.;gran=;telo=0;labels=none;gaps=1' \
   'A record that is entirely assembly gap'

fx st_all_telomere synthetic/st_all_telomere.fa \
   'chr_st_all_telomere=F:CCCTAAx500' '-' \
   'type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0' \
   'A record that is nothing but forward repeat is one array, not two arms'

fx st_gap_flanked synthetic/st_gap_flanked.fa \
   'chr_st_gap_flanked=N:200+F:CCCTAAx100+L:2400+R:TTAGGGx100+N:200' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=2' \
   'Telomeres behind leading and trailing gaps are still terminal'

fx st_x_mask synthetic/st_x_mask.fa \
   'chr_st_x_mask=F:CCCTAAx100+L:1200+X:200+L:1200+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=1' \
   'X is an assembly gap character too, and counts in the gaps BED'

# Soft-masked sequence is ordinary sequence: repeat masking must not hide a telomere.
fx st_soft_masked synthetic/st_soft_masked.fa \
   'chr_st_soft_masked=M:ccctaax100+L:2800+R:TTAGGGx100' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'A lowercase soft-masked p arm is still a telomere'

# Multi-record coverage across completeness buckets.
fx st_six_buckets synthetic/st_six_buckets.fa \
   'rec_t2t=F:CCCTAAx100+L:2800+R:TTAGGGx100;rec_t2t_gapped=F:CCCTAAx100+L:1300+N:200+L:1300+R:TTAGGGx100;rec_incomplete=F:CCCTAAx100+L:2800;rec_incomplete_gapped=F:CCCTAAx100+L:1300+N:200+L:1300;rec_none=L:3000;rec_none_gapped=L:1400+N:200+L:1400' '-' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0|type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=1|type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=0|type=incomplete;anom=.;gran=P;telo=1;labels=p;gaps=1|type=none;anom=.;gran=;telo=0;labels=none;gaps=0|type=none;anom=.;gran=;telo=0;labels=none;gaps=1' \
   'One record in each of the six completeness buckets'

# ---- motif and pattern handling ----
fx st_plant_pair synthetic/st_plant_pair.fa \
   'chr_st_plant_pair=F:CCCTAAAx86+L:2800+R:TTTAGGGx86' '-c CCCTAAA' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Plant 7-mer canonical on the standard geometry'

fx st_edit_zero synthetic/st_edit_zero.fa \
   'chr_st_edit_zero=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-x 0' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Exact canonical arrays are found with no edit distance allowed'

fx st_edit_two synthetic/st_edit_two.fa \
   'chr_st_edit_two=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-x 2' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'Raising the edit distance must not lose an exact canonical arm'

fx st_explicit_patterns synthetic/st_explicit_patterns.fa \
   'chr_st_explicit_patterns=F:CCCTAAx100+L:2800+R:TTAGGGx100' '-c TTAGGG -p TTAGGG,CCCTAA' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'An explicit -p set covering the canonical pair behaves like the derived one'

# ---- shapes that occur in real assemblies ----
fx st_variant_halo synthetic/st_variant_halo.fa \
   'chr_st_variant_halo=V:CTCTAAx20+F:CCCTAAx80+L:2800+R:TTAGGGx100' '-x 1' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0' \
   'A degenerate outer edge on an otherwise canonical arm, as subtelomeres really look'

# Terminal zone limit (-t) toggles terminal vs interstitial.
fx st_inner_array_its synthetic/st_inner_array.fa \
   'chr_st_inner_array=F:CCCTAAx100+L:600+F:CCCTAAx50+L:2400+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'At -t 1000 the array at offset 1200 is outside the terminal zone and is an ITS'

fx st_inner_array_terminal synthetic/st_inner_array.fa \
   'chr_st_inner_array=F:CCCTAAx100+L:600+F:CCCTAAx50+L:2400+R:TTAGGGx100' '-i' \
   'type=t2t;anom=misassembly;gran=PpQ;telo=2;labels=pq;gaps=0;its=0' \
   'The same array at the default tolerance of 2000 is an extra terminal block'

fx st_its_abuts_terminal synthetic/st_its_abuts_terminal.fa \
   'chr_st_its_abuts_terminal=F:CCCTAAx100+L:600+F:CCCTAAx50+L:2400+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'The inner array falls just outside the -t 1000 zone and abuts the p-arm terminal block as an ITS'

fx st_its_true_interior synthetic/st_its_true_interior.fa \
   'chr_st_its_true_interior=F:CCCTAAx100+L:3000+F:CCCTAAx100+L:3000+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=1' \
   'The same array moved to the true interior is an ITS'

fx st_fusion_headtohead synthetic/st_fusion_headtohead.fa \
   'chr_st_fusion_headtohead=F:CCCTAAx100+L:1200+R:TTAGGGx100+F:CCCTAAx100+L:1200+R:TTAGGGx100' '-i -t 1000' \
   'type=t2t;anom=.;gran=PQ;telo=2;labels=pq;gaps=0;its=?' \
   'An interior head-to-head junction, the signature of a chromosome fusion: REG-004'
