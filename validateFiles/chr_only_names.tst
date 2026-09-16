-f testFiles/chr_only_names.fa --chr-only -o %OUTDIR%
expect_exit 0
expect_stderr_substr selected 3 of 5 paths
expect_stderr_substr prefix 'chr', 0 separator(s)
expect_stdout testFiles/expected/chr_only_names.txt
