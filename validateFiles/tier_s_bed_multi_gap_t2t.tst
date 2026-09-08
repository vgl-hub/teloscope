-f testFiles/multi_gap_t2t.fa -i -n -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file multi_gap_t2t.fa_terminal_telomeres.bed testFiles/expected/tier_s/multi_gap_t2t.fa_terminal_telomeres.bed
expect_file multi_gap_t2t.fa_interstitial_telomeres.bed testFiles/expected/tier_s/multi_gap_t2t.fa_interstitial_telomeres.bed
