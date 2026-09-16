-f testFiles/synthetic/jn_contig_middle.fa -n -i -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file jn_contig_middle.fa_terminal_telomeres.bed testFiles/expected/tier_s/jn_contig_middle.fa_terminal_telomeres.bed
