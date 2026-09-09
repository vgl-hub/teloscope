-f testFiles/t2t.fa -p NNNNNNN -x 0 -o %OUTDIR%
expect_exit 1
expect_stdout ignore
expect_stderr_substr expands beyond 4096 combinations
