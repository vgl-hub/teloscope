-f testFiles/t2t.fa -w 100 -s 96 -r -g -e -m -i -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_file t2t.fa_terminal_telomeres.bed testFiles/expected/overlap/t2t.fa_terminal_telomeres.bed
expect_file t2t.fa_window_repeat_density.bedgraph testFiles/expected/overlap/t2t.fa_window_repeat_density.bedgraph
expect_file t2t.fa_canonical_matches.bed testFiles/expected/overlap/t2t.fa_canonical_matches.bed
expect_file t2t.fa_window_canonical_ratio.bedgraph testFiles/expected/overlap/t2t.fa_window_canonical_ratio.bedgraph
expect_file t2t.fa_window_strand_ratio.bedgraph testFiles/expected/overlap/t2t.fa_window_strand_ratio.bedgraph
expect_file t2t.fa_window_gc.bedgraph testFiles/expected/overlap/t2t.fa_window_gc.bedgraph
expect_file t2t.fa_window_entropy.bedgraph testFiles/expected/overlap/t2t.fa_window_entropy.bedgraph
expect_file t2t.fa_noncanonical_matches.bed testFiles/expected/overlap/t2t.fa_noncanonical_matches.bed
