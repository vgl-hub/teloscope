-f testFiles/t2t.fa -a -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file t2t.fa_terminal_telomeres.fa testFiles/expected/tier_s/t2t.fa_terminal_telomeres.fa
expect_file t2t.fa_terminal_telomeres.bed testFiles/expected/tier_s/t2t.fa_terminal_telomeres.bed
