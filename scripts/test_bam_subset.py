#!/usr/bin/env python3

import binascii
import os
import pathlib
import random
import struct
import subprocess
import sys
import tempfile
import zlib


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))
EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
BASE_CODES = {base: i for i, base in enumerate("=ACMGRSVTWYHKDBN")}


def bgzf_block(data, extra=b"", filename=b"", comment=b"", header_crc=False):
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = compressor.compress(data) + compressor.flush()
    xlen = 6 + len(extra)
    flags = 4 | (8 if filename else 0) | (16 if comment else 0) | (2 if header_crc else 0)
    header = bytearray(b"\x1f\x8b\x08" + bytes([flags]) + b"\x00\x00\x00\x00\x00\xff")
    header += struct.pack("<H", xlen)
    header += b"BC\x02\x00\x00\x00" + extra
    if filename:
        header += filename + b"\0"
    if comment:
        header += comment + b"\0"
    if header_crc:
        header += b"\0\0"
    total = len(header) + len(payload) + 8
    if total > 65536:
        raise ValueError("BGZF block is too large")
    struct.pack_into("<H", header, 16, total - 1)
    if header_crc:
        struct.pack_into("<H", header, len(header) - 2, binascii.crc32(header[:-2]) & 0xFFFF)
    footer = struct.pack("<II", binascii.crc32(data) & 0xFFFFFFFF, len(data))
    return bytes(header) + payload + footer


def bgzf(data, chunk_size=60000, extra=b"", eof=True, **header_fields):
    blocks = []
    for start in range(0, len(data), chunk_size):
        first_fields = header_fields if start == 0 else {}
        blocks.append(bgzf_block(data[start:start + chunk_size], extra if start == 0 else b"", **first_fields))
    if eof:
        blocks.append(EOF_BLOCK)
    return b"".join(blocks)


def unpack_bgzf(data):
    out = bytearray()
    pos = 0
    while pos < len(data):
        if data[pos:pos + 3] != b"\x1f\x8b\x08":
            raise AssertionError("invalid BGZF output")
        xlen = struct.unpack_from("<H", data, pos + 10)[0]
        extra_end = pos + 12 + xlen
        extra = data[pos + 12:extra_end]
        block_size = None
        sub = 0
        while sub < len(extra):
            length = struct.unpack_from("<H", extra, sub + 2)[0]
            if extra[sub:sub + 2] == b"BC":
                block_size = struct.unpack_from("<H", extra, sub + 4)[0] + 1
            sub += 4 + length
        if block_size is None:
            raise AssertionError("missing BGZF BC field")
        block = data[pos:pos + block_size]
        out += zlib.decompress(block[extra_end - pos:-8], -15)
        pos += block_size
    return bytes(out)


def pack_sequence(sequence):
    encoded = bytearray((len(sequence) + 1) // 2)
    for index, base in enumerate(sequence.upper()):
        code = BASE_CODES[base]
        if index % 2:
            encoded[index // 2] |= code
        else:
            encoded[index // 2] = code << 4
    return bytes(encoded)


def bam_record(name, sequence, flag=0, ref_id=-1, pos=-1, cigar=(), tags=b"", qualities=None):
    read_name = name.encode() + b"\0"
    cigar_data = b"".join(struct.pack("<I", (length << 4) | op) for length, op in cigar)
    if qualities is None:
        quality_data = b"\xff" * len(sequence)
    else:
        quality_data = bytes(qualities)
    bin_mq_nl = (0 << 16) | (60 << 8) | len(read_name)
    flag_nc = (flag << 16) | len(cigar)
    core = struct.pack(
        "<iiIIiiii",
        ref_id,
        pos,
        bin_mq_nl,
        flag_nc,
        len(sequence),
        -1,
        -1,
        0,
    )
    body = core + read_name + cigar_data + pack_sequence(sequence) + quality_data + tags
    return struct.pack("<i", len(body)) + body


def bam_payload(records, header_text=b"@HD\tVN:1.6\tSO:unknown\n", references=(("chr1", 1000000),)):
    reference_data = bytearray(struct.pack("<i", len(references)))
    for name, length in references:
        encoded_name = name.encode() + b"\0"
        reference_data += struct.pack("<i", len(encoded_name)) + encoded_name + struct.pack("<i", length)
    return b"BAM\1" + struct.pack("<i", len(header_text)) + header_text + reference_data + b"".join(records)


def split_bam(payload):
    if payload[:4] != b"BAM\1":
        raise AssertionError("invalid BAM payload")
    text_length = struct.unpack_from("<i", payload, 4)[0]
    pos = 8 + text_length
    reference_count = struct.unpack_from("<i", payload, pos)[0]
    pos += 4
    for _ in range(reference_count):
        name_length = struct.unpack_from("<i", payload, pos)[0]
        pos += 4 + name_length + 4
    header = payload[:pos]
    records = []
    while pos < len(payload):
        block_size = struct.unpack_from("<i", payload, pos)[0]
        end = pos + 4 + block_size
        records.append(payload[pos:end])
        pos = end
    return header, records


def run(args, stdin=None, cwd=None):
    try:
        return subprocess.run(
            [str(TELOSCOPE), *[str(arg) for arg in args]],
            input=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=cwd,
            check=False,
            timeout=30,
        )
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"Teloscope timed out: {' '.join(str(a) for a in args)}") from error


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def record_name(record):
    name_length = record[12]
    return record[36:36 + name_length - 1].decode()


# Reads mode always writes files under -o; the report also prints to stdout.

def run_reads(out_dir, options, data, path=None):
    """Run teloscope over reads-mode input, either from a file at `path` or from stdin.
    Returns (result, input_name) where input_name is the base used for output filenames."""
    if path is None:
        result = run([*options, "-o", str(out_dir)], stdin=data)
        return result, "stdin"
    path.write_bytes(data)
    result = run([*options, str(path), "-o", str(out_dir)])
    return result, path.name


def bed_rows(out_dir, name):
    bed_path = out_dir / f"{name}_terminal_telomeres.bed"
    if not bed_path.exists():
        return []
    return [line.split("\t") for line in bed_path.read_text().splitlines() if line]


def report_text(out_dir, name):
    return (out_dir / f"{name}_report.tsv").read_text()


def fastq_subset_bytes(out_dir, name):
    return (out_dir / f"{name}_telomeric.fastq").read_bytes()


def same_records(actual, golden):
    # Windows reads text files with CRLF translated, as the retired validator's goldens did
    if os.name == "nt":
        actual, golden = actual.replace(b"\r\n", b"\n"), golden.replace(b"\r\n", b"\n")
    return actual == golden


def bam_subset_path(out_dir, name):
    return out_dir / f"{pathlib.Path(name).stem}_telomeric.bam"


def bam_subset_records(out_dir, name):
    path = bam_subset_path(out_dir, name)
    require(path.exists(), f"missing BAM output: {path}")
    return split_bam(unpack_bgzf(path.read_bytes()))


def fastq_payload(sequences):
    data = bytearray()
    for name, sequence in sequences.items():
        data += f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()
    return bytes(data)


# the keep-scan floor is fixed at 42bp regardless of -l; -y/-k/-d/-x still gate it
def explicit_args():
    return ["-x", "0", "-y", "0.8", "-k", "10", "-d", "10"]


def assert_fastq_bam_parity(tmp, label, sequences, options):
    bam_dir = tmp / f"{label}_bam_out"
    bam_records_input = [bam_record(name, sequence, flag=0x4) for name, sequence in sequences.items()]
    bam_result, bam_name = run_reads(bam_dir, options, bgzf(bam_payload(bam_records_input)))
    require(bam_result.returncode == 0, bam_result.stderr.decode())
    _, bam_out_records = bam_subset_records(bam_dir, bam_name)
    bam_names = [record_name(record) for record in bam_out_records]

    fastq_dir = tmp / f"{label}_fastq_out"
    fastq_result, fastq_name = run_reads(fastq_dir, options, fastq_payload(sequences))
    require(fastq_result.returncode == 0, fastq_result.stderr.decode())
    fastq_bytes = fastq_subset_bytes(fastq_dir, fastq_name)
    fastq_names = [line[1:].decode() for line in fastq_bytes.splitlines()[::4]]
    require(bam_names == fastq_names, f"FASTQ/BAM scoring parity failed ({label}): {bam_names} vs {fastq_names}")
    return bam_names


def test_basic_and_record_preservation(tmp):
    # 8 repeat units (48bp) clear the fixed 42bp keep floor regardless of -l
    passing = [
        bam_record("mapped_p", "CCCTAA" * 8, ref_id=0, pos=10, cigar=((48, 0),), tags=b"NM\x69\x00\x00\x00\x00"),
        bam_record("reverse_q", "TTAGGG" * 8, flag=0x10, ref_id=0, pos=20, cigar=((48, 0),)),
        bam_record("secondary", "TTAGGG" * 8, flag=0x100),
        bam_record("supplementary", "CCCTAA" * 8, flag=0x800),
        bam_record("odd_ambiguous", "N" + "TTAGGG" * 8, flag=0x4),
        bam_record("paired_first", "TTAGGG" * 8, flag=0x41),
        bam_record("paired_second", "CCCTAA" * 8, flag=0x81),
        bam_record("all_codes", "=ACMGRSVTWYHKDBN" + "TTAGGG" * 8, flag=0x4),
    ]
    failing = [
        bam_record("ordinary", "ACGTACGTACGTACGTAC", flag=0x4),
        bam_record("missing", "", flag=0x4),
    ]
    input_payload = bam_payload(
        [passing[0], failing[0], passing[1], passing[2], failing[1], *passing[3:]],
        references=(("chr1", 1000000), ("chrM", 16569)),
    )
    input_bam = bgzf(
        input_payload,
        chunk_size=37,
        extra=b"XY\x03\x00abc",
        filename=b"reads.bam",
        comment=b"fixture",
        header_crc=True,
    )
    out_dir = tmp / "basic_out"
    result, name = run_reads(out_dir, explicit_args(), input_bam)
    require(result.returncode == 0, result.stderr.decode())
    header, records = bam_subset_records(out_dir, name)
    expected_header, _ = split_bam(input_payload)
    require(header == expected_header, "BAM header changed")
    require(records == passing, "passing BAM records were not preserved exactly")
    require(b"kept 8 of 10 records" in result.stderr, "summary count mismatch")
    require(b"skipped 1 record without SEQ" in result.stderr, "missing SEQ count mismatch")


def test_measured_vs_kept_counts(tmp):
    # measured (BED-row) rows skip secondary/supplementary/hard-clipped, but all can be kept
    sequence = "TTAGGG" * 100 + "A" * 2000
    records = [
        bam_record("primary", sequence, flag=0x4),
        bam_record("secondary", sequence, flag=0x104),
        bam_record("supplementary", sequence, flag=0x804),
        bam_record("clipped", sequence, flag=0x4, ref_id=0, pos=0,
                   cigar=((10, 5), (len(sequence) - 10, 0))),  # op 5 = hard clip, leading end
    ]
    out_dir = tmp / "measured_vs_kept_out"
    result, name = run_reads(out_dir, [], bgzf(bam_payload(records)))
    require(result.returncode == 0, result.stderr.decode())
    require(b"1 secondary/supplementary record not measured" not in result.stderr and
            b"2 secondary/supplementary records not measured" in result.stderr,
            "secondary/supplementary not-measured count is wrong")
    require(b"1 hard-clipped record not measured" in result.stderr, "hard-clip not-measured count is wrong")
    require(b"kept 4 of 4 records" in result.stderr, "all four records should still be kept")
    rows = bed_rows(out_dir, name)
    require([row[0] for row in rows] == ["primary"], f"only the primary record should be measured: {rows}")


def test_file_output_and_threads(tmp):
    records = []
    expected = []
    for index in range(700):
        if index % 4:
            record = bam_record(f"pass_{index}", "TTAGGG" * (7 + index % 7), flag=0x4)
            expected.append(record)
        else:
            record = bam_record(f"fail_{index}", "ACGT" * 20, flag=0x4)
        records.append(record)
    input_path = tmp / "many.bam"
    input_path.write_bytes(bgzf(bam_payload(records), chunk_size=113))

    single_dir, multi_dir = tmp / "many_j1", tmp / "many_j8"
    single = run([str(input_path), "-j", "1", "-o", str(single_dir)])
    multi = run([str(input_path), "-j", "8", "-o", str(multi_dir)])
    require(single.returncode == 0, single.stderr.decode())
    require(multi.returncode == 0, multi.stderr.decode())
    single_bytes = bam_subset_path(single_dir, input_path.name).read_bytes()
    multi_bytes = bam_subset_path(multi_dir, input_path.name).read_bytes()
    require(single_bytes == multi_bytes, "thread count changed BAM output")
    _, actual = split_bam(unpack_bgzf(single_bytes))
    require(actual == expected, "large-batch filtering mismatch")
    require(b"BAM\1" not in single.stdout, "stdout carried binary BAM bytes instead of the text report")


def test_multiblock_determinism(tmp):
    generator = random.Random(7)
    records = []
    expected = []
    for index in range(640):
        sequence = "".join(generator.choices("ACGT", k=16000))
        if index % 5 == 0:
            sequence = sequence[:-300] + "TTAGGG" * 50
            records.append(bam_record(f"pass_{index}", sequence, flag=0x4))
            expected.append(records[-1])
        else:
            records.append(bam_record(f"fail_{index}", sequence, flag=0x4))
    input_path = tmp / "multiblock.bam"
    input_path.write_bytes(bgzf(bam_payload(records)))

    single_dir, multi_dir = tmp / "multi_j1", tmp / "multi_j8"
    single = run([str(input_path), "-j", "1", "-o", str(single_dir)])
    multi = run([str(input_path), "-j", "8", "-o", str(multi_dir)])
    require(single.returncode == 0, single.stderr.decode())
    require(multi.returncode == 0, multi.stderr.decode())
    single_bytes = bam_subset_path(single_dir, input_path.name).read_bytes()
    multi_bytes = bam_subset_path(multi_dir, input_path.name).read_bytes()
    require(single_bytes == multi_bytes, "thread count changed multi-block BAM output")
    _, actual = split_bam(unpack_bgzf(single_bytes))
    require(actual == expected, "multi-block filtering mismatch")


def test_header_only_and_missing_eof(tmp):
    payload = bam_payload([])
    out_dir = tmp / "header_only_out"
    result, name = run_reads(out_dir, [], bgzf(payload, eof=False))
    require(result.returncode == 0, result.stderr.decode())
    require(b"missing the BGZF EOF marker" in result.stderr, "missing EOF warning absent")
    header, records = bam_subset_records(out_dir, name)
    require(header == payload and records == [], "header-only BAM changed")


def test_header_variants_and_empty_results(tmp):
    passing = [
        bam_record("first", "TTAGGG" * 8, flag=0x4),
        bam_record("second", "CCCTAA" * 8, flag=0x4),
    ]
    payload = bam_payload(
        passing,
        header_text=b"",
        references=(("a", 1), ("long_reference_name", 2**31 - 1)),
    )
    out_dir = tmp / "variant_out"
    result, name = run_reads(out_dir, [], bgzf(payload))
    header, records = bam_subset_records(out_dir, name)
    expected_header, _ = split_bam(payload)
    require(header == expected_header and records == passing, "BAM header variant changed")

    no_reference_payload = bam_payload(passing, header_text=b"@CO\tempty references\n", references=())
    no_reference_dir = tmp / "no_reference_out"
    no_reference, name2 = run_reads(no_reference_dir, [], bgzf(no_reference_payload))
    header, records = bam_subset_records(no_reference_dir, name2)
    expected_header, _ = split_bam(no_reference_payload)
    require(header == expected_header and records == passing, "reference-free BAM changed")

    failing = [bam_record(f"fail_{index}", "ACGT" * (index + 2), flag=0x4) for index in range(5)]
    fail_payload = bam_payload(failing)
    all_fail_dir = tmp / "all_fail_out"
    all_fail, name3 = run_reads(all_fail_dir, [], bgzf(fail_payload))
    header, records = bam_subset_records(all_fail_dir, name3)
    expected_header, _ = split_bam(fail_payload)
    require(header == expected_header and records == [], "all-fail BAM was not header-only")
    require(b"kept 0 of 5 records" in all_fail.stderr, "all-fail count mismatch")


def test_embedded_eof(tmp):
    records = [
        bam_record("first", "TTAGGG" * 8, flag=0x4),
        bam_record("second", "CCCTAA" * 8, flag=0x4),
    ]
    payload = bam_payload(records)
    split = len(payload) // 2
    input_bam = bgzf(payload[:split]) + bgzf(payload[split:])
    out_dir = tmp / "embedded_eof_out"
    result, name = run_reads(out_dir, [], input_bam)
    require(result.returncode == 0, result.stderr.decode())
    _, actual = bam_subset_records(out_dir, name)
    require(actual == records, "embedded BGZF EOF interrupted input")


def test_cross_block_large_record(tmp):
    generator = random.Random(7)
    sequence = "".join(generator.choice("ACGT") for _ in range(70000)) + "TTAGGG" * 10
    qualities = [generator.randrange(94) for _ in sequence]
    record = bam_record("large", sequence, flag=0x4, qualities=qualities)
    payload = bam_payload([record])
    out_dir = tmp / "cross_block_out"
    result, name = run_reads(out_dir, [], bgzf(payload, chunk_size=60000))
    require(result.returncode == 0, result.stderr.decode())
    _, records = bam_subset_records(out_dir, name)
    require(records == [record], "record spanning BGZF blocks changed")


def test_byte_bounded_batch(tmp):
    aux_count = 33 * 1024 * 1024
    huge_aux = b"ZZBC" + struct.pack("<i", aux_count) + bytes(aux_count)
    huge = bam_record("huge_aux", "TTAGGG" * 8, flag=0x4, tags=huge_aux)
    tail = bam_record("tail", "CCCTAA" * 8, flag=0x4)
    out_dir = tmp / "byte_bounded_out"
    result, name = run_reads(out_dir, ["-j", "2"], bgzf(bam_payload([huge, tail])))
    require(result.returncode == 0, result.stderr.decode())
    _, records = bam_subset_records(out_dir, name)
    require(records == [huge, tail], "byte-bounded batch changed records")


def test_default_threshold_parity(tmp):
    # the fixed keep floor is 42bp: 6 repeats (36bp) fail, 7 (42bp) and up pass
    sequences = {
        "short": "TTAGGG" * 6,
        "default_pass": "TTAGGG" * 7,
        "long": "CCCTAA" * 15,
        "fail": "ACGT" * 20,
    }
    names = assert_fastq_bam_parity(tmp, "default_threshold", sequences, [])
    require(names == ["default_pass", "long"], "default threshold boundary failed")

    sequence = "TTAGGG" * 10
    crlf = f"@crlf\r\n{sequence}\r\n+\r\n{'I' * len(sequence)}\r\n".encode()
    crlf_dir = tmp / "crlf_stdin_out"
    result, name = run_reads(crlf_dir, [], crlf)
    require(result.returncode == 0, result.stderr.decode())
    require(fastq_subset_bytes(crlf_dir, name).startswith(b"@crlf\r\n"), "CRLF FASTQ filtering failed")


def test_l_changes_measurement_not_subset(tmp):
    # -l no longer gates the subset (fixed at 42bp); it still gates the measured BED rows
    fixture = ROOT / "testFiles" / "fastq_subset_threshold.fq"
    low_dir, high_dir = tmp / "l_low_out", tmp / "l_high_out"
    low = run([str(fixture), "-y", "0.8", "-k", "10", "-d", "10", "-l", "20", "-o", str(low_dir)])
    high = run([str(fixture), "-y", "0.8", "-k", "10", "-d", "10", "-l", "5000", "-o", str(high_dir)])
    require(low.returncode == 0, low.stderr.decode())
    require(high.returncode == 0, high.stderr.decode())
    require(
        fastq_subset_bytes(low_dir, fixture.name) == fastq_subset_bytes(high_dir, fixture.name),
        "-l changed which reads were kept",
    )
    require(len(bed_rows(low_dir, fixture.name)) == 3, "-l 20 should measure all three reads")
    require(len(bed_rows(high_dir, fixture.name)) == 0, "-l 5000 should measure none of the three reads")
    golden = (ROOT / "testFiles" / "expected" / "fastq_subset_threshold_default.fq").read_bytes()
    require(same_records(fastq_subset_bytes(low_dir, fixture.name), golden), "kept reads changed from the checked-in golden")


def test_checked_in_fastq_fixtures(tmp):
    # the file-output equivalent of the retired fastq_subset_*.tst goldens
    expected_dir = ROOT / "testFiles" / "expected"
    default_golden = (expected_dir / "fastq_subset.fq").read_bytes()

    for name, path in (
        ("blanklines", ROOT / "testFiles" / "fastq_subset_blanklines.fq"),
        ("file", ROOT / "testFiles" / "fastq_subset.fq"),
        ("gz", ROOT / "testFiles" / "fastq_subset.fq.gz"),
    ):
        out_dir = tmp / f"fixture_{name}_out"
        result = run([str(path), "-o", str(out_dir)])
        require(result.returncode == 0, f"{name}: {result.stderr.decode()}")
        require(same_records(fastq_subset_bytes(out_dir, path.name), default_golden), f"{name}: kept reads changed")

    # stdin, piped from the same fixture
    stdin_dir = tmp / "fixture_stdin_out"
    stdin_result, stdin_name = run_reads(stdin_dir, [], (ROOT / "testFiles" / "fastq_subset.fq").read_bytes())
    require(stdin_result.returncode == 0, stdin_result.stderr.decode())
    require(same_records(fastq_subset_bytes(stdin_dir, stdin_name), default_golden), "stdin: kept reads changed")

    # CRLF: its own golden, since the kept records carry \r\n
    crlf_path = ROOT / "testFiles" / "fastq_subset_crlf.fq"
    crlf_golden = (expected_dir / "fastq_subset_crlf.fq").read_bytes()
    crlf_out = tmp / "fixture_crlf_out"
    crlf_result = run([str(crlf_path), "-o", str(crlf_out)])
    require(crlf_result.returncode == 0, crlf_result.stderr.decode())
    require(same_records(fastq_subset_bytes(crlf_out, crlf_path.name), crlf_golden), "crlf: kept reads changed")

    # realistic: a wider mix of terminal/internal/softmasked/N-flanked/scattered reads
    realistic_path = ROOT / "testFiles" / "fastq_subset_realistic.fq"
    realistic_golden = (expected_dir / "fastq_subset_realistic.fq").read_bytes()
    realistic_out = tmp / "fixture_realistic_out"
    realistic_result = run([str(realistic_path), "-y", "0.8", "-k", "10", "-d", "10", "-o", str(realistic_out)])
    require(realistic_result.returncode == 0, realistic_result.stderr.decode())
    require(same_records(fastq_subset_bytes(realistic_out, realistic_path.name), realistic_golden),
            "realistic: kept reads changed")


def test_fastq_malformed_fixture_aborts(tmp):
    path = ROOT / "testFiles" / "fastq_malformed.fq"
    out_dir = tmp / "malformed_out"
    result = run([str(path), "-o", str(out_dir)])
    require(result.returncode == 1, "malformed FASTQ should exit 1")
    require(b"sequence and quality length differ" in result.stderr, "malformed FASTQ diagnostic changed")
    require(not out_dir.exists() or not list(out_dir.iterdir()), "malformed FASTQ left partial output")


def test_fastq_file_output_and_threads(tmp):
    # the FASTQ equivalent of test_file_output_and_threads: batching boundaries and -j parity
    lines = []
    expected_names = []
    for index in range(700):
        if index % 4:
            sequence = "TTAGGG" * (7 + index % 7)
            expected_names.append(f"pass_{index}")
        else:
            sequence = "ACGT" * 20
        lines.append(f"@{'pass' if index % 4 else 'fail'}_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    input_path = tmp / "many.fq"
    input_path.write_text("".join(lines))

    single_dir, multi_dir = tmp / "many_j1", tmp / "many_j8"
    single = run([str(input_path), "-j", "1", "-o", str(single_dir)])
    multi = run([str(input_path), "-j", "8", "-o", str(multi_dir)])
    require(single.returncode == 0, single.stderr.decode())
    require(multi.returncode == 0, multi.stderr.decode())
    single_bytes = fastq_subset_bytes(single_dir, input_path.name)
    multi_bytes = fastq_subset_bytes(multi_dir, input_path.name)
    require(single_bytes == multi_bytes, "thread count changed FASTQ output")
    kept_names = [line[1:] for line in single_bytes.decode().splitlines()[::4]]
    require(kept_names == expected_names, "large-batch FASTQ filtering mismatch")


def test_exact_math_boundaries(tmp):
    # the keep floor sits at exactly 7 canonical repeat units (42bp); a threshold match is kept
    lengths = {"six_repeats": "TTAGGG" * 6, "seven_repeats": "TTAGGG" * 7}
    names = assert_fastq_bam_parity(tmp, "unit_boundary", lengths, ["-x", "0", "-y", "1", "-k", "10", "-d", "10"])
    require(names == ["seven_repeats"], "42bp/7-unit boundary failed")

    # density boundary, scaled up so the whole block clears the 42bp floor either way
    density = {"two_thirds": "TTAGGG" * 6 + "AAAAAA" * 6 + "TTAGGG" * 6}
    pass_options = ["-x", "0", "-y", "0.666", "-k", "40", "-d", "40"]
    fail_options = ["-x", "0", "-y", "0.667", "-k", "40", "-d", "40"]
    require(assert_fastq_bam_parity(tmp, "density_pass", density, pass_options) == ["two_thirds"],
            "density lower boundary failed")
    require(assert_fastq_bam_parity(tmp, "density_fail", density, fail_options) == [],
            "density upper boundary failed")

    plant = {
        "plant_pass": "TTTAGGG" * 6,
        "vertebrate_fail": "TTAGGG" * 7,
    }
    plant_options = ["-c", "CCCTAAA", "-x", "0", "-y", "1"]
    require(assert_fastq_bam_parity(tmp, "plant", plant, plant_options) == ["plant_pass"],
            "custom canonical failed")


def test_randomized_fastq_bam_parity(tmp):
    generator = random.Random(23)
    sequences = {}
    for index in range(240):
        length = generator.randrange(18, 250)
        sequence = "".join(generator.choice("ACGTN") for _ in range(length))
        if index % 3 == 0:
            insert = generator.randrange(len(sequence) + 1)
            repeat = generator.choice(("TTAGGG", "CCCTAA")) * generator.randrange(2, 18)
            sequence = sequence[:insert] + repeat + sequence[insert:]
        if index % 11 == 0:
            sequence += "TCAGGG" * 8 + "TTAGGG"
        sequences[f"random_{index:03d}"] = sequence

    option_sets = [
        ["-x", "0", "-y", "0.8", "-k", "10", "-d", "10"],
        ["-x", "1", "-y", "0.5", "-k", "50", "-d", "50"],
        ["-x", "0", "-y", "1", "-k", "10", "-d", "10"],
    ]
    for index, options in enumerate(option_sets):
        assert_fastq_bam_parity(tmp, f"randomized_{index}", sequences, options)


def test_cli_guards_and_cleanup(tmp):
    record = bam_record("pass", "TTAGGG" * 8, flag=0x4)
    input_bam = bgzf(bam_payload([record]))
    cmd_dir = tmp / "cmd_out"
    command = run(["--cmd", "-o", str(cmd_dir)], input_bam)
    require(command.returncode == 0, command.stderr.decode())
    _, records = bam_subset_records(cmd_dir, "stdin")
    require(records == [record], "--cmd run corrupted the BAM output file")
    require(b"teloscope" in command.stdout, "--cmd did not echo the command line to stdout")

    for flag in ("--bam-subset", "--fastq-subset"):
        removed = run([flag], input_bam)
        require(removed.returncode != 0, f"{flag} unexpectedly succeeded")
        require(
            b"--fastq-subset and --bam-subset were removed" in removed.stderr,
            f"{flag} lacked the removal diagnostic",
        )

    input_path = tmp / "broken.bam"
    input_path.write_bytes(input_bam[:20])
    out_dir = tmp / "broken_out"
    failed = run([str(input_path), "-o", str(out_dir)])
    require(failed.returncode != 0, "truncated file unexpectedly succeeded")
    require(not list(out_dir.iterdir()), "partial reads-mode output was retained")

    if os.name != "nt" and os.geteuid() != 0:
        valid_path = tmp / "valid.bam"
        valid_path.write_bytes(input_bam)
        unwritable = tmp / "unwritable"
        unwritable.mkdir()
        unwritable.chmod(0o500)
        try:
            failed = run([str(valid_path), "-o", str(unwritable)])
            require(failed.returncode != 0, "unwritable output unexpectedly succeeded")
            require(b"is not writable" in failed.stderr, "output open failure was unclear")
        finally:
            unwritable.chmod(0o700)


def test_detection_boundary(tmp):
    # detection is a literal byte check: '@' -> FASTQ, 'BAM\1' -> BAM, else the bytes must still look like FASTA/GFA
    record = bam_record("pass", "TTAGGG" * 8, flag=0x4)
    payload = bam_payload([record])

    not_bgzf_path = tmp / "not_bgzf.bam"
    not_bgzf_path.write_bytes(b"not bam")
    out_dir = tmp / "not_bgzf_out"
    result = run([str(not_bgzf_path), "-o", str(out_dir)])
    require(result.returncode != 0, "raw bytes that are not gzip, FASTA or GFA should be refused, not run as an assembly")
    require(b"is not FASTA, GFA, FASTQ or BAM" in result.stderr, "non-BAM bytes lacked the detection diagnostic")
    require(not list(out_dir.iterdir()), "non-BAM bytes left output behind")

    bad_magic_path = tmp / "bad_magic.bam"
    bad_magic_path.write_bytes(bgzf(b"BAD\1" + payload[4:]))
    out_dir2 = tmp / "bad_magic_out"
    result2 = run([str(bad_magic_path), "-o", str(out_dir2)])
    require(result2.returncode != 0, "valid BGZF with the wrong magic bytes should be refused, not run as an assembly")
    require(b"is not FASTA, GFA, FASTQ or BAM" in result2.stderr, "mismatched magic bytes lacked the detection diagnostic")
    require(not list(out_dir2.iterdir()), "mismatched magic bytes left output behind")

    # detection cannot even open an unreadable file to sniff it, so it also falls through
    if os.name != "nt" and os.geteuid() != 0:
        unreadable_path = tmp / "unreadable.bam"
        unreadable_path.write_bytes(bgzf(payload))
        unreadable_path.chmod(0)
        out_dir3 = tmp / "unreadable_out"
        try:
            result3 = run([str(unreadable_path), "-o", str(out_dir3)])
            require(result3.returncode == 0, "an unreadable file should fall through to assembly mode")
            require(not list(out_dir3.glob("*_telomeric.*")), "an unreadable file was treated as reads-mode input")
        finally:
            unreadable_path.chmod(0o600)


# inputs that keep the literal "BAM\1" magic: detected as BAM, must abort cleanly, no partial output
def test_failures(tmp):
    record = bam_record("pass", "TTAGGG" * 8, flag=0x4)
    payload = bam_payload([record])
    valid = bytearray(bgzf(payload))
    header, records = split_bam(payload)
    record_offset = len(header)

    cases = {}
    gzip_compressor = zlib.compressobj(6, zlib.DEFLATED, 31)
    cases["plain_gzip"] = gzip_compressor.compress(payload) + gzip_compressor.flush()  # gzip, not BGZF-chunked
    cases["truncated"] = bytes(valid[:20])

    # a corrupt second block surfaces from the threaded pipeline, not the single-block sniff
    many_records = [bam_record(f"multi_{i}", "TTAGGG" * 8, flag=0x4) for i in range(50)]
    multiblock = bytearray(bgzf(bam_payload(many_records), chunk_size=200))
    first_block_len = struct.unpack_from("<H", multiblock, 16)[0] + 1
    second_block_len = struct.unpack_from("<H", multiblock, first_block_len + 16)[0] + 1
    second_block_crc_offset = first_block_len + second_block_len - 8
    multiblock[second_block_crc_offset] ^= 1
    cases["second_block_bad_crc"] = bytes(multiblock)

    # a valid block followed by non-gzip bytes: the reader thread's block-boundary check
    cases["trailing_junk"] = bgzf(payload, eof=False) + b"NOTAGZIPBLOCKTRAILINGJUNK1234567890"

    corrupt_crc = bytearray(valid)
    first_size = struct.unpack_from("<H", corrupt_crc, 16)[0] + 1
    corrupt_crc[first_size - 8] ^= 1
    cases["bad_crc"] = bytes(corrupt_crc)

    corrupt_size = bytearray(valid)
    corrupt_size[first_size - 4] ^= 1
    cases["bad_isize"] = bytes(corrupt_size)

    missing_bc = bytearray(valid)
    missing_bc[12:14] = b"XY"
    cases["missing_bc"] = bytes(missing_bc)

    duplicate_bc = bgzf(payload, extra=b"BC\x02\x00\x00\x00")
    cases["duplicate_bc"] = duplicate_bc

    malformed_extra = bgzf(payload, extra=b"X")
    cases["malformed_extra"] = malformed_extra

    malformed_subfield = bgzf(payload, extra=b"XY\x05\x00Z")
    cases["malformed_extra_subfield"] = malformed_subfield

    bad_block_size = bytearray(valid)
    struct.pack_into("<H", bad_block_size, 16, 1)
    cases["bad_bgzf_block_size"] = bytes(bad_block_size)

    reserved_flag = bytearray(valid)
    reserved_flag[3] |= 0x20
    cases["reserved_gzip_flag"] = bytes(reserved_flag)

    bad_header_crc = bytearray(bgzf(payload, header_crc=True))
    bad_header_crc[18] ^= 1
    cases["bad_header_crc"] = bytes(bad_header_crc)

    truncated_header_crc = bytearray(bgzf_block(b"", header_crc=True))
    struct.pack_into("<H", truncated_header_crc, 16, 26)
    cases["truncated_header_crc"] = bytes(truncated_header_crc)

    unterminated_filename = bytearray(bgzf_block(payload, filename=b"x"))
    filename_footer = len(unterminated_filename) - 8
    for index in range(18, filename_footer):
        if unterminated_filename[index] == 0:
            unterminated_filename[index] = 1
    cases["unterminated_filename"] = bytes(unterminated_filename)

    unterminated_comment = bytearray(bgzf_block(payload, comment=b"x"))
    comment_footer = len(unterminated_comment) - 8
    for index in range(18, comment_footer):
        if unterminated_comment[index] == 0:
            unterminated_comment[index] = 1
    cases["unterminated_comment"] = bytes(unterminated_comment)

    bad_deflate = bytearray(valid)
    bad_deflate[18] ^= 0x80
    cases["bad_deflate"] = bytes(bad_deflate)

    oversized_uncompressed = bytearray(valid)
    struct.pack_into("<I", oversized_uncompressed, first_size - 4, 65537)
    cases["oversized_uncompressed_block"] = bytes(oversized_uncompressed)

    cases["negative_header_length"] = bgzf(b"BAM\1" + struct.pack("<i", -1))
    cases["truncated_header_text"] = bgzf(b"BAM\1" + struct.pack("<i", 100) + b"short")
    cases["negative_reference_count"] = bgzf(b"BAM\1" + struct.pack("<i", 0) + struct.pack("<i", -1))
    cases["zero_reference_name"] = bgzf(
        b"BAM\1" + struct.pack("<i", 0) + struct.pack("<i", 1) + struct.pack("<i", 0)
    )
    cases["unterminated_reference_name"] = bgzf(
        b"BAM\1" + struct.pack("<i", 0) + struct.pack("<i", 1) +
        struct.pack("<i", 4) + b"chr1" + struct.pack("<i", 100)
    )
    cases["negative_reference_length"] = bgzf(
        b"BAM\1" + struct.pack("<i", 0) + struct.pack("<i", 1) +
        struct.pack("<i", 2) + b"x\0" + struct.pack("<i", -1)
    )

    bad_record = bytearray(payload)
    struct.pack_into("<i", bad_record, record_offset, 31)
    cases["bad_record_size"] = bgzf(bytes(bad_record))

    oversized_record = header + struct.pack("<I", 256 * 1024 * 1024 + 1)
    cases["oversized_record"] = bgzf(oversized_record)

    negative_record = header + struct.pack("<i", -1)
    cases["negative_record"] = bgzf(negative_record)

    truncated_record_size = header + b"\x20\x00"
    cases["truncated_record_size"] = bgzf(truncated_record_size)

    undersized_record_core = header + struct.pack("<i", 31) + bytes(31)
    cases["undersized_record_core"] = bgzf(undersized_record_core)

    bad_offsets = bytearray(payload)
    bad_offsets[record_offset + 12] = 250
    cases["bad_record_offsets"] = bgzf(bytes(bad_offsets))

    zero_read_name = bytearray(payload)
    zero_read_name[record_offset + 12] = 0
    cases["zero_read_name"] = bgzf(bytes(zero_read_name))

    unterminated_read_name = bytearray(payload)
    read_name_length = unterminated_read_name[record_offset + 12]
    unterminated_read_name[record_offset + 36 + read_name_length - 1] = ord("X")
    cases["unterminated_read_name"] = bgzf(bytes(unterminated_read_name))

    negative_sequence_length = bytearray(payload)
    struct.pack_into("<i", negative_sequence_length, record_offset + 20, -1)
    cases["negative_sequence_length"] = bgzf(bytes(negative_sequence_length))

    oversized_sequence = bytearray(payload)
    struct.pack_into("<i", oversized_sequence, record_offset + 20, 10000)
    cases["oversized_sequence"] = bgzf(bytes(oversized_sequence))

    oversized_cigar = bytearray(payload)
    struct.pack_into("<H", oversized_cigar, record_offset + 16, 65535)
    cases["oversized_cigar"] = bgzf(bytes(oversized_cigar))

    truncated_record = bytearray(records[0][:-5])
    struct.pack_into("<i", truncated_record, 0, len(truncated_record) - 4)
    cases["truncated_record_fields"] = bgzf(header + bytes(truncated_record))

    for name, data in cases.items():
        out_dir = tmp / f"fail_{name}"
        result = run(["-o", str(out_dir)], data)
        require(result.returncode != 0, f"{name} unexpectedly succeeded")
        require(b"Error:" in result.stderr, f"{name} lacked a clear error")
        require(not list(out_dir.iterdir()), f"{name} left partial output behind")


def test_gzipped_non_bam_stdin_aborts_cleanly(tmp):
    # gzip's leading byte routes stdin straight to BAM mode; a non-BAM payload trips the magic check inside
    out_dir = tmp / "gz_non_bam_out"
    result, _ = run_reads(out_dir, [], bgzf(b"this is plain text, not a BAM payload, once BGZF-decompressed"))
    require(result.returncode == 1, f"expected exit 1, got {result.returncode}")
    require(b"Error: invalid BAM magic." in result.stderr, f"missing BAM magic error: {result.stderr!r}")
    require(b"Compressed FASTQ or FASTA on stdin or a pipe is not supported" in result.stderr,
            f"missing compressed-stdin hint: {result.stderr!r}")
    require(not list(out_dir.iterdir()), "gzipped non-BAM stdin left partial output behind")


def test_mutation_robustness(tmp):
    case_count = int(os.environ.get("BAM_MUTATION_CASES", "48"))
    generator = random.Random(91)
    records = [
        bam_record(f"record_{index}", "TTAGGG" * (8 + index % 5), flag=0x4)
        for index in range(8)
    ]
    payload = bam_payload(records)
    header, split_records = split_bam(payload)
    record_offset = len(header)

    for index in range(case_count):
        mode = index % 7
        if mode == 0:
            # real magic + random tail, so this still routes to BAM (see test_detection_boundary)
            size = generator.randrange(0, 2048)
            garbage = bytes(generator.getrandbits(8) for _ in range(size))
            data = bgzf(b"BAM\1" + garbage)
        elif mode == 1:
            data = bgzf(payload[:generator.randrange(len(payload) + 1)])
        elif mode == 2:
            mutated = bytearray(payload)
            target = generator.randrange(record_offset + 36, len(mutated))
            mutated[target] ^= 1 << generator.randrange(8)
            data = bgzf(bytes(mutated))
        elif mode == 3:
            mutated = bytearray(payload)
            struct.pack_into("<i", mutated, record_offset, generator.randrange(-16, 129))
            data = bgzf(bytes(mutated))
        elif mode == 4:
            mutated = bytearray(payload)
            mutated[record_offset + 12] = generator.randrange(256)
            data = bgzf(bytes(mutated))
        elif mode == 5:
            mutated = bytearray(payload)
            struct.pack_into("<H", mutated, record_offset + 16, generator.randrange(65536))
            data = bgzf(bytes(mutated))
        else:
            kept = split_records[:generator.randrange(len(split_records) + 1)]
            data = bgzf(header + b"".join(kept))

        out_dir = tmp / f"mutation_{index}"
        result = run(["-j", "1", "-o", str(out_dir)], data)
        if result.returncode == 0:
            # either a real BAM decode, or content that fell through to the assembly path
            bam_path = bam_subset_path(out_dir, "stdin")
            if bam_path.exists():
                split_bam(unpack_bgzf(bam_path.read_bytes()))
        else:
            require(b"Error:" in result.stderr, f"mutation {index} (mode {mode}) lacked a clear error")
            require(not list(out_dir.iterdir()), f"mutation {index} (mode {mode}) left partial output behind")


def main():
    require(TELOSCOPE.exists(), f"Teloscope binary not found: {TELOSCOPE}")
    with tempfile.TemporaryDirectory(prefix="teloscope_bam_") as temp:
        tmp = pathlib.Path(temp)
        test_basic_and_record_preservation(tmp)
        test_measured_vs_kept_counts(tmp)
        test_file_output_and_threads(tmp)
        test_multiblock_determinism(tmp)
        test_header_only_and_missing_eof(tmp)
        test_header_variants_and_empty_results(tmp)
        test_embedded_eof(tmp)
        test_cross_block_large_record(tmp)
        test_byte_bounded_batch(tmp)
        test_default_threshold_parity(tmp)
        test_l_changes_measurement_not_subset(tmp)
        test_checked_in_fastq_fixtures(tmp)
        test_fastq_malformed_fixture_aborts(tmp)
        test_fastq_file_output_and_threads(tmp)
        test_exact_math_boundaries(tmp)
        test_randomized_fastq_bam_parity(tmp)
        test_cli_guards_and_cleanup(tmp)
        test_detection_boundary(tmp)
        test_failures(tmp)
        test_gzipped_non_bam_stdin_aborts_cleanly(tmp)
        test_mutation_robustness(tmp)
    print("PASS BAM subset integration")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL BAM subset integration: {error}", file=sys.stderr)
        raise
