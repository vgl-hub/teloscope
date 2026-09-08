-f testFiles/its_headtohead.fa -i -n -t 50 -x 1 -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file its_headtohead.fa_terminal_telomeres.bed testFiles/expected/tier_s/its_headtohead.fa_terminal_telomeres.bed
expect_file its_headtohead.fa_interstitial_telomeres.bed testFiles/expected/tier_s/its_headtohead.fa_interstitial_telomeres.bed
