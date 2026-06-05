--fastq-subset -x 0 -l 18 -y 0.8 -k 10 -d 10 -j 8 testFiles/fastq_subset_large.fq
expect_exit 0
expect_stdout testFiles/expected/fastq_subset_large.fq
expect_stderr_substr FASTQ subset: kept 240 of 320 reads.
