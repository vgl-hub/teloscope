-f testFiles/vgp_turtle.fa.gz -i -n -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file vgp_turtle.fa.gz_terminal_telomeres.bed testFiles/expected/tier_s/vgp_turtle.fa.gz_terminal_telomeres.bed
