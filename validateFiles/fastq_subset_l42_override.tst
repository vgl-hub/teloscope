--fastq-subset -x 0 -y 0.8 -k 10 -d 10 -l 42 testFiles/fastq_subset_threshold.fq
expect_exit 0
expect_stdout testFiles/expected/fastq_subset_threshold_l42.fq
expect_stderr_substr FASTQ subset: kept 3 of 3 reads.
