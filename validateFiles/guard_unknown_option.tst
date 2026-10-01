-f testFiles/t2t.fa --plot-reprot -o %OUTDIR%
expect_exit 1
expect_stdout ignore
expect_stderr_substr Unknown or ambiguous option --plot-reprot
