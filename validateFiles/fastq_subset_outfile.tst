--fastq-subset -x 0 -l 18 -y 0.8 -k 10 -d 10 testFiles/fastq_subset.fq -o %OUTDIR%
expect_exit 0
expect_stdout testFiles/expected/fastq_subset_empty.fq
expect_stderr_substr Wrote telomeric reads to
