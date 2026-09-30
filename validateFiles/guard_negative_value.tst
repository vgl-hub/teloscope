-f testFiles/t2t.fa -w -5 -o %OUTDIR%
expect_exit 1
expect_stdout ignore
expect_stderr_substr -w/--window must be > 0
