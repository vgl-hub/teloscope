-f testFiles/its.fa -i -n -t 1000 -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file its.fa_terminal_telomeres.bed testFiles/expected/tier_s/its.fa_terminal_telomeres.bed
expect_file its.fa_interstitial_telomeres.bed testFiles/expected/tier_s/its.fa_interstitial_telomeres.bed
