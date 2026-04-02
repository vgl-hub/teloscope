-f testFiles/gfa_path_orient_pairs_small.gfa -j 1 -x 0 -l 60 -o %OUTDIR%
expect_exit 0
expect_stdout ignore
expect_output_name gfa_path_orient_pairs_small.gfa.telo.annotated.gfa
expect_gfa_header 1.2
gfa_expect testFiles/expected/gfa/gfa_path_orient_pairs_small.tsv
gfa_preserve_input subset
