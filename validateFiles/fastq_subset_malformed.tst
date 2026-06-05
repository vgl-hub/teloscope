--fastq-subset -x 0 -l 6 testFiles/fastq_malformed.fq
expect_exit 1
expect_stdout ignore
expect_stderr_substr sequence and quality length differ
