#!/usr/bin/env python3

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))
SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")


def run(args, **kwargs):
    try:
        return subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=30,
            **kwargs,
        )
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"command timed out: {' '.join(map(str, args))}") from error


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def samtools_view(path, header=False):
    args = [SAMTOOLS, "view", "--no-PG"]
    if header:
        args.append("-H")
    args.append(str(path))
    result = run(args)
    require(result.returncode == 0, result.stderr.decode())
    return result.stdout.decode().splitlines()


def main():
    require(TELOSCOPE.exists(), f"Teloscope binary not found: {TELOSCOPE}")
    require(shutil.which(SAMTOOLS) is not None, f"samtools not found: {SAMTOOLS}")

    qualities = "I" * 18
    sam_lines = [
        "@HD\tVN:1.6\tSO:unknown",
        "@SQ\tSN:chr1\tLN:1000000",
        "@CO\tgenerated independently by samtools",
        f"mapped_p\t0\tchr1\t11\t60\t18M\t*\t0\t0\tCCCTAACCCTAACCCTAA\t{qualities}\tNM:i:0\tZZ:Z:kept",
        f"ordinary\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGTAC\t{qualities}\tZZ:Z:dropped",
        f"reverse_q\t16\tchr1\t21\t40\t18M\t*\t0\t0\tTTAGGGTTAGGGTTAGGG\t{qualities}\tAS:i:18",
        f"secondary\t260\t*\t0\t0\t*\t*\t0\t0\tTTAGGGTTAGGGTTAGGG\t{qualities}\tXA:Z:secondary",
        "missing\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
        f"supplementary\t2048\tchr1\t31\t30\t18M\t*\t0\t0\tCCCTAACCCTAACCCTAA\t{qualities}\tSA:Z:chr1,31,+,18M,30,0;",
    ]
    expected_names = ["mapped_p", "reverse_q", "secondary", "supplementary"]

    with tempfile.TemporaryDirectory(prefix="teloscope_samtools_") as temp:
        tmp = pathlib.Path(temp)
        sam_path = tmp / "oracle.sam"
        input_bam = tmp / "oracle.bam"
        output_dir = tmp / "output"
        output_bam = output_dir / "oracle_telomeric.bam"
        sam_path.write_text("\n".join(sam_lines) + "\n", encoding="ascii")

        converted = run([SAMTOOLS, "view", "--no-PG", "-b", "-o", str(input_bam), str(sam_path)])
        require(converted.returncode == 0, converted.stderr.decode())

        result = run([
            str(TELOSCOPE),
            "--bam-subset",
            "-x", "0",
            "-l", "18",
            "-y", "0.8",
            "-k", "10",
            "-d", "10",
            "-j", "2",
            "-o", str(output_dir),
            str(input_bam),
        ])
        require(result.returncode == 0, result.stderr.decode())
        require(output_bam.exists(), "Teloscope did not create the BAM output")
        require(result.stdout == b"", "file output wrote BAM bytes to stdout")

        checked = run([SAMTOOLS, "quickcheck", "-v", str(output_bam)])
        require(checked.returncode == 0, checked.stdout.decode() + checked.stderr.decode())
        require(samtools_view(input_bam, header=True) == samtools_view(output_bam, header=True),
                "samtools observed a changed BAM header")

        input_records = {line.split("\t", 1)[0]: line for line in samtools_view(input_bam)}
        output_records = samtools_view(output_bam)
        actual_names = [line.split("\t", 1)[0] for line in output_records]
        require(actual_names == expected_names, "samtools observed the wrong output records")
        require(output_records == [input_records[name] for name in expected_names],
                "samtools observed changed alignment records")
        require(not list(output_dir.glob("*.bai")), "Teloscope created a BAM index")
        require(not list(output_dir.glob("*.csi")), "Teloscope created a CSI index")

    print("PASS BAM samtools interoperability")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL BAM samtools interoperability: {error}", file=sys.stderr)
        raise
