-f testFiles/t2t.fa -w 5 -s 5 -r -o %OUTDIR%
expect_exit 1
expect_stdout ignore
expect_stderr_substr is smaller than pattern
