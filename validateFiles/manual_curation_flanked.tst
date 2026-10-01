-f testFiles/synthetic/jn_flanked_both.fa -n -i -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file jn_flanked_both.fa_terminal_telomeres.bed testFiles/expected/tier_s/jn_flanked_both.fa_terminal_telomeres.bed
