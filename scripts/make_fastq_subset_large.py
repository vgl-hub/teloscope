#!/usr/bin/env python3
"""Generate a deterministic large FASTQ for the --fastq-subset determinism test.

Reads alternate telomeric (pass) and non-telomeric (fail) so the emitted subset is a
non-trivial subsequence spanning multiple batches and thread chunks. Unique headers
make any reordering detectable by the order-sensitive validator. Output is written to
stdout; redirect it to testFiles/fastq_subset_large.fq.

Intended for: teloscope --fastq-subset -x 0 -l 18 -y 0.8 -k 10 -d 10
"""
import sys

N = 320  # > recordsPerBatch (256) so at least two batches run at any thread count

def main():
    lines = []
    for i in range(N):
        is_fail = (i % 4 == 3)              # 1 in 4 reads is non-telomeric
        if is_fail:
            length = 24 + (i % 4) * 6       # 24..42 bp, no TTAGGG/CCCTAA 6-mer
            seq = ("ACGT" * 12)[:length]
            tag = "fail"
        else:
            motif = "TTAGGG" if (i % 2 == 0) else "CCCTAA"  # exercise both strands
            seq = motif * (4 + (i % 4))     # 24..42 bp telomeric tract
            tag = "pass"
        lines.append(f"@read_{i:06d} {tag}")
        lines.append(seq)
        lines.append("+")
        lines.append("I" * len(seq))
    sys.stdout.write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
