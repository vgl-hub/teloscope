--fastq-subset -x 0 -y 0.8 -k 10 -d 10 testFiles/fastq_subset_threshold.fq
expect_exit 0
expect_stdout testFiles/expected/fastq_subset_threshold_default.fq
expect_stderr_substr FASTQ subset: kept 2 of 3 reads.
