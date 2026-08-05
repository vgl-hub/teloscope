-f testFiles/discordant.fa -i -n -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file discordant.fa_terminal_telomeres.bed testFiles/expected/tier_s/discordant.fa_terminal_telomeres.bed
expect_file discordant.fa_interstitial_telomeres.bed testFiles/expected/tier_s/discordant.fa_interstitial_telomeres.bed
