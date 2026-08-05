-f testFiles/misassembly.fa -i -n -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file misassembly.fa_terminal_telomeres.bed testFiles/expected/tier_s/misassembly.fa_terminal_telomeres.bed
expect_file misassembly.fa_interstitial_telomeres.bed testFiles/expected/tier_s/misassembly.fa_interstitial_telomeres.bed
