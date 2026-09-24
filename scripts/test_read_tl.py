#!/usr/bin/env python3

import gzip
import math
import os
import pathlib
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import test_bam_subset as bamlib  # reuses the BAM/BGZF synthesiser


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))

FWD = "CCCTAA"  # p: forward canonical, expected at a read start
REV = "TTAGGG"  # q: reverse canonical, expected at a read end
MAX_BLOCK_DIST = 1000  # -d default


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def run(args, stdin=None):
    result = subprocess.run(
        [str(TELOSCOPE), *[str(arg) for arg in args]],
        input=stdin,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=60,
    )
    return result


def read_bed(path):
    if not path.exists():
        return []
    text = path.read_text()
    rows = [line.split("\t") for line in text.splitlines() if line]
    return rows


def revcomp(sequence):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement[base] for base in reversed(sequence))


def run_fastq_reads(tmp, sequences, extra_args=(), prefix="reads"):
    # writes a FASTQ, runs it (positional path, no flag), and returns rows/report/stdout/kept
    path = tmp / f"{prefix}.fq"
    with open(path, "w") as f:
        for name, sequence in sequences.items():
            f.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    out = tmp / f"{prefix}_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out), *extra_args])
    require(result.returncode == 0, f"reads-mode FASTQ run failed: {result.stderr.decode()}")
    rows = read_bed(out / f"{path.name}_terminal_telomeres.bed")
    report = (out / f"{path.name}_report.tsv").read_text()
    kept_path = out / f"{path.name}_telomeric.fastq"
    kept_names = [line[1:] for line in kept_path.read_text().splitlines()[::4]] if kept_path.exists() else []
    return rows, report, result.stdout.decode(), kept_names


def run_bam_reads(tmp, sequences, extra_args=(), prefix="reads", flags=None):
    names = list(sequences)
    flags = flags or {name: 0x4 for name in names}
    records = [bamlib.bam_record(name, sequences[name], flag=flags[name]) for name in names]
    path = tmp / f"{prefix}.bam"
    path.write_bytes(bamlib.bgzf(bamlib.bam_payload(records)))
    out = tmp / f"{prefix}_bam_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out), *extra_args])
    require(result.returncode == 0, f"reads-mode BAM run failed: {result.stderr.decode()}")
    rows = read_bed(out / f"{path.name}_terminal_telomeres.bed")
    bam_path = out / f"{path.stem}_telomeric.bam"
    _, kept_records = bamlib.split_bam(bamlib.unpack_bgzf(bam_path.read_bytes()))
    kept_names = [bamlib.record_name(record) for record in kept_records]
    return rows, kept_names


def percentile_type7(sorted_values, p):
    # R type 7 / numpy default: linear interpolation between the two closest ranks
    h = (len(sorted_values) - 1) * p
    lo = math.floor(h)
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (h - lo) * (sorted_values[hi] - sorted_values[lo])


def row_by_name(rows, name):
    return [row for row in rows if row[0] == name]


def make_p_read(tl_len, tip_junk=0, tail_len=2000):
    # p: canonical forward motif at the read start, tip_junk bp of unrelated bases before it
    units = round(tl_len / len(FWD))
    return "A" * tip_junk + FWD * units + "G" * tail_len


def make_q_read(tl_len, tip_junk=0, head_len=2000):
    # q: canonical reverse motif at the read end, tip_junk bp of unrelated bases after it
    units = round(tl_len / len(REV))
    return "A" * head_len + REV * units + "A" * tip_junk


def test_p_and_q_lengths(tmp):
    sequences = {}
    expected_len = {}
    for tl_len in (300, 1002, 3000, 10002):
        pname = f"p_{tl_len}"
        qname = f"q_{tl_len}"
        sequences[pname] = make_p_read(tl_len)
        sequences[qname] = make_q_read(tl_len)
        expected_len[pname] = tl_len
        expected_len[qname] = tl_len

    rows, _, _, _ = run_fastq_reads(tmp, sequences, prefix="lengths")
    require(len(rows) == len(sequences), f"expected one row per read, got {len(rows)}")
    for row in rows:
        name = row[0]
        teloLen = int(row[3])
        require(teloLen == expected_len[name], f"{name}: expected TL {expected_len[name]}, got {teloLen}")
        label = row[4]
        closest = row[5]
        require(label == closest, f"{name}: expected concordant, got label={label} closestEnd={closest}")
        expectedSide = 'p' if name.startswith('p_') else 'q'
        require(closest == expectedSide, f"{name}: expected closestEnd={expectedSide}, got {closest}")


def test_read_tl_percentiles(tmp):
    # complete-read length sets for n = 1..5; all lengths are multiples of len(FWD), so make_p_read's rounding is a no-op and the BED teloLen equals the requested length
    cases = {
        1: [600],
        2: [300, 900],
        3: [300, 900, 2100],
        4: [300, 600, 1200, 2400],
        5: [300, 600, 900, 1200, 1500],  # hand-checked below
    }
    for n, lengths in cases.items():
        sequences = {f"p_{i}": make_p_read(tl_len) for i, tl_len in enumerate(lengths)}
        rows, report_text, stdout_text, _ = run_fastq_reads(tmp, sequences, prefix=f"pct_n{n}")
        require(len(rows) == n, f"n={n}: expected {n} rows, got {len(rows)}")
        teloLens = sorted(int(row[3]) for row in rows)
        require(teloLens == sorted(lengths), f"n={n}: BED teloLen mismatch: {teloLens} vs {lengths}")
        require("Reads kept:" in report_text, "report is missing the 'Reads kept:' line")
        require(f"Reads kept:\t{n}" in report_text, f"n={n}: wrong Reads kept value")

        values = sorted(lengths)
        expected = {
            "Mean length": sum(values) / n,
            "Median length": percentile_type7(values, 0.5),
            "25th percentile length": percentile_type7(values, 0.25),
            "75th percentile length": percentile_type7(values, 0.75),
            "90th percentile length": percentile_type7(values, 0.90),
        }
        for label, value in expected.items():
            line = f"{label}:\t{value:.2f}"
            require(line in report_text, f"n={n}: report missing {line!r}")
            require(line in stdout_text, f"n={n}: stdout missing {line!r}")
        require(f"Min length:\t{values[0]}" in report_text, f"n={n}: wrong Min length")
        require(f"Max length:\t{values[-1]}" in report_text, f"n={n}: wrong Max length")

    # hand-checked case: lengths 300/600/900/1200/1500, n=5, h = 4p
    hand_report = run_fastq_reads(tmp, {f"h_{i}": make_p_read(l)
                                         for i, l in enumerate([300, 600, 900, 1200, 1500])},
                                   prefix="pct_hand")[1]
    require("Mean length:\t900.00" in hand_report, "hand-checked mean wrong")
    require("Median length:\t900.00" in hand_report, "hand-checked median wrong")
    require("25th percentile length:\t600.00" in hand_report, "hand-checked p25 wrong")
    require("75th percentile length:\t1200.00" in hand_report, "hand-checked p75 wrong")
    require("90th percentile length:\t1380.00" in hand_report, "hand-checked p90 wrong")
    require("Min length:\t300" in hand_report, "hand-checked min wrong")
    require("Max length:\t1500" in hand_report, "hand-checked max wrong")

    # mixed p/q case exercising a non-trivial interpolation at both p25 and p75
    mixed_lengths = [300, 600, 1002, 1800, 3000, 4200]
    mixed_sequences = {}
    for i, l in enumerate(mixed_lengths):
        mixed_sequences[f"mp_{i}"] = make_p_read(l) if i % 2 == 0 else make_q_read(l)
    mixed_report = run_fastq_reads(tmp, mixed_sequences, prefix="pct_mixed")[1]
    values = sorted(mixed_lengths)
    for label, p in (("25th percentile length", 0.25), ("75th percentile length", 0.75)):
        line = f"{label}:\t{percentile_type7(values, p):.2f}"
        require(line in mixed_report, f"mixed p/q case missing {line!r}")


def test_discordant(tmp):
    # reverse motif planted at the read start: teloLabel 'q' != closestEnd 'p'
    sequences = {"discordant_p": REV * 100 + "A" * 2000}
    rows, report, _, _ = run_fastq_reads(tmp, sequences, prefix="discordant")
    require(len(rows) == 1, "expected exactly one row")
    row = rows[0]
    require(row[4] == 'q' and row[5] == 'p', f"expected discordant q/p row, got {row[4]}/{row[5]}")
    require("Discordant:\t1" in report, "report did not count the discordant row")


def test_reaching_read_end(tmp):
    # the telomere runs to the very end of the read: flank < -d, so it is not "complete"
    sequences = {"short_p": make_p_read(600, tail_len=200)}
    rows, report, _, _ = run_fastq_reads(tmp, sequences, prefix="reaching_end")
    require(len(rows) == 1, "expected exactly one row")
    require("Complete:\t0" in report, "short flank was wrongly counted complete")
    require("Reaching read end:\t1" in report, "short flank was not counted as reaching the read end")


def test_interstitial_only_no_row(tmp):
    # telomeric block sitting in the read's interior: no terminal row is written
    sequences = {"interstitial": "A" * 500 + FWD * 60 + "A" * 500}
    rows, _, _, kept_names = run_fastq_reads(tmp, sequences, prefix="interstitial")
    require(len(rows) == 0, f"interstitial-only read produced a terminal row: {rows}")
    # the OR rule: no row, but the (much larger) keep-scan zone still finds the block
    require(kept_names == ["interstitial"], "interstitial-only read should still be kept")


def test_or_rule_kept_is_superset_of_measured(tmp):
    # every measured read is kept; the reverse is not true (kept-with-zero-rows, shown above)
    sequences = {
        "p_short": make_p_read(300),
        "q_long": make_q_read(3000),
        "discordant": REV * 80 + "A" * 2000,
        "reaching_end": make_p_read(600, tail_len=100),
        "interstitial": "A" * 400 + FWD * 60 + "A" * 400,
        "no_telomere": "A" * 2000,
    }
    rows, _, _, kept_names = run_fastq_reads(tmp, sequences, prefix="or_rule")
    measured_names = {row[0] for row in rows}
    require(measured_names, "expected at least one measured row in this mix")
    require(measured_names.issubset(set(kept_names)),
            f"measured names {measured_names} are not all kept ({kept_names})")
    require("no_telomere" not in kept_names, "a read with no telomeric content at all should not be kept")


def test_tip_junk_tolerance(tmp):
    # default read tip allowance is 300bp; junk past that hides the telomere entirely
    sequences = {
        "junk_0": make_p_read(600, tip_junk=0),
        "junk_150": make_p_read(600, tip_junk=150),
        "junk_400": make_p_read(600, tip_junk=400),
        "q_junk_150": make_q_read(600, tip_junk=150),
    }
    rows, _, _, _ = run_fastq_reads(tmp, sequences, prefix="tip_junk")
    require(row_by_name(rows, "junk_0"), "junk_0 should have a row")
    p150 = row_by_name(rows, "junk_150")
    require(p150 and p150[0][1] == "150" and p150[0][3] == "600", f"junk_150 should start at 150 with 600 bp: {p150}")
    q150 = row_by_name(rows, "q_junk_150")
    require(q150 and q150[0][2] == str(2000 + 600) and q150[0][3] == "600", f"q_junk_150 should end 150 bp from the tip with 600 bp: {q150}")
    require(not row_by_name(rows, "junk_400"), "junk_400 should have no row past the default tolerance")


def test_lowercase(tmp):
    sequences = {"lower_p": make_p_read(600).lower()}
    rows, _, _, _ = run_fastq_reads(tmp, sequences, prefix="lowercase")
    require(len(rows) == 1 and int(rows[0][3]) == 600, "lowercase read was not scanned like uppercase")


def test_crlf(tmp):
    sequence = make_p_read(600)
    path = tmp / "crlf.fq"
    with open(path, "wb") as f:
        f.write(f"@crlf_p\r\n{sequence}\r\n+\r\n{'I' * len(sequence)}\r\n".encode())
    out = tmp / "crlf_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 0, result.stderr.decode())
    rows = read_bed(out / f"{path.name}_terminal_telomeres.bed")
    require(len(rows) == 1, "CRLF FASTQ produced the wrong row count")
    require(rows[0][0] == "crlf_p", f"CRLF header parsing kept stray characters: {rows[0][0]!r}")
    require(int(rows[0][3]) == 600, "CRLF read length was miscounted")


def test_gz(tmp):
    sequence = make_p_read(600)
    path = tmp / "reads.fq.gz"
    with gzip.open(path, "wt") as f:
        f.write(f"@gz_p\n{sequence}\n+\n{'I' * len(sequence)}\n")
    out = tmp / "gz_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 0, result.stderr.decode())
    rows = read_bed(out / f"{path.name}_terminal_telomeres.bed")
    require(len(rows) == 1 and int(rows[0][3]) == 600, "gzipped FASTQ was not scanned correctly")


def test_fastq_bam_parity(tmp):
    sequences = {
        "p_short": make_p_read(300),
        "q_long": make_q_read(3000),
        "discordant": REV * 80 + "A" * 2000,
        "reaching_end": make_p_read(600, tail_len=100),
        "interstitial": "A" * 400 + FWD * 60 + "A" * 400,
    }
    fastq_rows, _, _, _ = run_fastq_reads(tmp, sequences, prefix="parity")
    bam_rows, _ = run_bam_reads(tmp, sequences, prefix="parity")
    require(fastq_rows == bam_rows,
            f"FASTQ and unaligned-BAM read TL rows differ:\nFASTQ={fastq_rows}\nBAM={bam_rows}")


def test_bam_reverse_strand_restores_sequencing_orientation(tmp):
    # a 0x10 record stores SEQ reverse-complemented; the row must match the un-reversed read
    sequencing_orientation = make_p_read(600)
    fastq_rows, _, _, _ = run_fastq_reads(tmp, {"rev_p": sequencing_orientation}, prefix="orient_fastq")

    stored = revcomp(sequencing_orientation)
    bam_rows, _ = run_bam_reads(tmp, {"rev_p": stored}, prefix="orient_bam", flags={"rev_p": 0x10})
    require(fastq_rows == bam_rows,
            f"reverse-strand BAM row should match the FASTQ row: {bam_rows} vs {fastq_rows}")


def test_bam_secondary_supplementary_kept_not_measured(tmp):
    sequences = {"primary": make_p_read(600)}
    records = [
        bamlib.bam_record("primary", sequences["primary"], flag=0x4),
        bamlib.bam_record("secondary", sequences["primary"], flag=0x104),
        bamlib.bam_record("supplementary", sequences["primary"], flag=0x804),
    ]
    path = tmp / "flags.bam"
    path.write_bytes(bamlib.bgzf(bamlib.bam_payload(records)))
    out = tmp / "flags_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 0, result.stderr.decode())
    rows = read_bed(out / f"{path.name}_terminal_telomeres.bed")
    require([row[0] for row in rows] == ["primary"],
            f"secondary/supplementary records should not be measured: {rows}")
    _, kept_records = bamlib.split_bam(bamlib.unpack_bgzf((out / "flags_telomeric.bam").read_bytes()))
    kept_names = [bamlib.record_name(record) for record in kept_records]
    require(kept_names == ["primary", "secondary", "supplementary"],
            f"secondary/supplementary records should still be kept: {kept_names}")


def test_bam_hard_clip_skips_measurement_but_can_be_kept(tmp):
    sequence = make_p_read(600)
    record = bamlib.bam_record(
        "clipped_p", sequence, flag=0x4, ref_id=0, pos=0,
        cigar=((10, 5), (len(sequence) - 10, 0)),  # 5 = hard clip, at the leading end
    )
    path = tmp / "hardclip.bam"
    path.write_bytes(bamlib.bgzf(bamlib.bam_payload([record])))
    out = tmp / "hardclip_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 0, result.stderr.decode())
    require(b"1 hard-clipped record not measured" in result.stderr, "hard-clip count was not reported on stderr")
    rows = read_bed(out / f"{path.name}_terminal_telomeres.bed")
    require(rows == [], f"a hard-clipped record should not be measured: {rows}")
    _, kept_records = bamlib.split_bam(bamlib.unpack_bgzf((out / "hardclip_telomeric.bam").read_bytes()))
    require([bamlib.record_name(r) for r in kept_records] == ["clipped_p"],
            "a hard-clipped record should still be kept when the keep scan passes")


def test_bam_reads_measured_and_kept_counts(tmp):
    # measured only counts the plain primary record; kept also passes secondary/supplementary/hard-clipped ones
    sequence = make_p_read(600)
    records = [
        bamlib.bam_record("primary", sequence, flag=0x4),
        bamlib.bam_record("secondary", sequence, flag=0x104),
        bamlib.bam_record("supplementary", sequence, flag=0x804),
        bamlib.bam_record("clipped", sequence, flag=0x4, ref_id=0, pos=0,
                           cigar=((10, 5), (len(sequence) - 10, 0))),
    ]
    path = tmp / "counts.bam"
    path.write_bytes(bamlib.bgzf(bamlib.bam_payload(records)))
    out = tmp / "counts_bam_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 0, result.stderr.decode())
    report = (out / f"{path.name}_report.tsv").read_text()
    require("Reads measured:\t1" in report, f"expected 1 measured read, got: {report}")
    require("Reads kept:\t4" in report, f"expected Reads kept to exceed Reads measured: {report}")


def test_thread_count_identical(tmp):
    random.seed(11)
    sequences = {}
    for index in range(120):
        tl_len = random.choice([300, 600, 1000, 3000])
        kind = index % 3
        if kind == 0:
            sequences[f"p_{index}"] = make_p_read(tl_len)
        elif kind == 1:
            sequences[f"q_{index}"] = make_q_read(tl_len)
        else:
            sequences[f"none_{index}"] = "A" * 400 + FWD * 40 + "A" * 400

    rows1, _, _, _ = run_fastq_reads(tmp, sequences, extra_args=["-j", "1"], prefix="threads1")
    rows8, _, _, _ = run_fastq_reads(tmp, sequences, extra_args=["-j", "8"], prefix="threads8")
    require(rows1 == rows8, "-j 1 and -j 8 produced different read TL rows")
    bam1, _ = run_bam_reads(tmp, sequences, extra_args=["-j", "1"], prefix="threads1")
    bam8, _ = run_bam_reads(tmp, sequences, extra_args=["-j", "8"], prefix="threads8")
    require(bam1 == bam8 == rows1, "BAM rows changed with -j or differ from FASTQ")


def test_tiled_matches_whole_read_scan(tmp):
    # the key fast-mode invariant: growing the tile in scanSegment must never change
    # the result versus scanning the whole read in one shot (-t default vs -t 1000000)
    random.seed(29)

    def junk(n):
        return "".join(random.choice("ACGT") for _ in range(n))

    sequences = {}
    for index in range(40):
        motif = FWD if index % 2 == 0 else REV
        parts = [motif * random.randint(60, 400)]
        for _ in range(random.randint(1, 4)):
            gap = random.randint(5, 950)  # within -d, so the chain should keep bridging
            parts += [junk(gap), motif * random.randint(20, 400)]
        parts.append(junk(random.randint(0, 1500)))  # pushes some reads past the 2kb default tile
        sequences[f"long_{index}"] = "".join(parts if motif == FWD else parts[::-1])

    rows_default, _, _, _ = run_fastq_reads(tmp, sequences, prefix="tiled")
    rows_whole, _, _, _ = run_fastq_reads(tmp, sequences, extra_args=["-t", "5000000"], prefix="whole")
    require(rows_default == rows_whole,
            "tiled fast-mode scan (default -t) differs from a whole-read scan (-t 5000000)")


def test_empty_input_error(tmp):
    path = tmp / "empty.fq"
    path.write_bytes(b"")
    out = tmp / "empty_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 1, f"expected exit 1 for empty input, got {result.returncode}")
    require(f"{path.name}' is empty.".encode() in result.stderr,  # the path is printed canonical (macOS /private/var)
            f"expected the empty-input diagnostic, got: {result.stderr!r}")


def test_empty_stdin_error(tmp):
    result = run([], stdin=b"")
    require(result.returncode == 1, f"expected exit 1 for empty stdin, got {result.returncode}")
    require(b"Error: input on stdin is empty." in result.stderr,
            f"expected the empty-stdin diagnostic, got: {result.stderr!r}")


def test_reads_measured_counts_every_fastq_read(tmp):
    # "Reads measured" counts every FASTQ read, whether or not it has telomeric content
    sequences = {
        "p_read": make_p_read(600),
        "discordant": REV * 80 + "A" * 2000,
        "no_telomere": "A" * 2000,
    }
    _, report, _, _ = run_fastq_reads(tmp, sequences, prefix="measured_fastq")
    require(f"Reads measured:\t{len(sequences)}" in report,
            f"expected Reads measured to count every FASTQ read: {report}")


def test_no_complete_telomeres_statistics_omitted(tmp):
    # all-discordant and all-reaching-end reads: no complete concordant telomere, so no length stats
    sequences = {
        "discordant": REV * 80 + "A" * 2000,
        "reaching_end": make_p_read(600, tail_len=100),
    }
    rows, report, stdout_text, _ = run_fastq_reads(tmp, sequences, prefix="no_stats")
    require(len(rows) == 2, f"expected two measured rows, got {rows}")
    require("No complete concordant telomeres for statistics." in report, "missing the empty-statistics line")
    require("No complete concordant telomeres for statistics." in stdout_text,
            "empty-statistics line missing from stdout")
    for label in ("Mean length", "Median length", "25th percentile length", "75th percentile length",
                  "90th percentile length", "Min length", "Max length"):
        require(f"{label}:" not in report, f"{label} row should be omitted with no complete telomeres")


def test_output_path_blocked_by_existing_directory(tmp):
    # a directory already sitting at the BED path blocks that ofstream from opening
    in_dir = tmp / "blocked_in"
    in_dir.mkdir()
    sequence = make_p_read(600)
    path = in_dir / "blocked.fq"
    path.write_text(f"@p_600\n{sequence}\n+\n{'I' * len(sequence)}\n")
    out = tmp / "blocked_out"
    out.mkdir()
    (out / f"{path.name}_terminal_telomeres.bed").mkdir()
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 1, f"expected exit 1, got {result.returncode}")
    require(f"Error: cannot write telomeric records to '{out}'.".encode() in result.stderr,
            f"expected the cannot-write diagnostic, got: {result.stderr!r}")
    require(not list(out.iterdir()), f"blocked BED path left outputs behind: {list(out.iterdir())}")


def test_scanned_kept_summary_line_removed(tmp):
    sequence = make_p_read(600)
    path = tmp / "summary_line.fq"
    path.write_text(f"@p_600\n{sequence}\n+\n{'I' * len(sequence)}\n")
    out = tmp / "summary_line_out"
    out.mkdir(exist_ok=True)
    result = run([str(path), "-o", str(out)])
    require(result.returncode == 0, result.stderr.decode())
    require(b"Reads:" not in result.stderr, "expected the 'Reads: N measured, M kept...' line to be gone")


def main():
    require(TELOSCOPE.exists(), f"Teloscope binary not found: {TELOSCOPE}")
    with tempfile.TemporaryDirectory(prefix="teloscope_read_tl_") as temp:
        tmp = pathlib.Path(temp)
        test_p_and_q_lengths(tmp)
        test_read_tl_percentiles(tmp)
        test_discordant(tmp)
        test_reaching_read_end(tmp)
        test_interstitial_only_no_row(tmp)
        test_or_rule_kept_is_superset_of_measured(tmp)
        test_tip_junk_tolerance(tmp)
        test_lowercase(tmp)
        test_crlf(tmp)
        test_gz(tmp)
        test_fastq_bam_parity(tmp)
        test_bam_reverse_strand_restores_sequencing_orientation(tmp)
        test_bam_secondary_supplementary_kept_not_measured(tmp)
        test_bam_hard_clip_skips_measurement_but_can_be_kept(tmp)
        test_bam_reads_measured_and_kept_counts(tmp)
        test_thread_count_identical(tmp)
        test_tiled_matches_whole_read_scan(tmp)
        test_empty_input_error(tmp)
        test_empty_stdin_error(tmp)
        test_reads_measured_counts_every_fastq_read(tmp)
        test_no_complete_telomeres_statistics_omitted(tmp)
        test_output_path_blocked_by_existing_directory(tmp)
        test_scanned_kept_summary_line_removed(tmp)
    print("PASS read TL integration")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL read TL integration: {error}", file=sys.stderr)
        raise
