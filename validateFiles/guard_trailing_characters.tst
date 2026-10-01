-f testFiles/t2t.fa -t 5e4 -o %OUTDIR%
expect_exit 1
expect_stdout ignore
expect_stderr_substr Invalid value '5e4' for -t/--terminal-limit
