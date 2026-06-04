--fastq-subset -x 0 -l 18 -y 0.8 -k 10 -d 10 < testFiles/fastq_subset.fq
expect_exit 0
expect_stdout testFiles/expected/fastq_subset.fq
expect_stderr_substr FASTQ subset: kept 2 of 4 reads.
