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


def run(args, stdin=None):
    try:
        return subprocess.run(
            [str(TELOSCOPE), *args],
            input=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=30,
        )
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"Teloscope timed out: {' '.join(args)}") from error


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def record_name(record):
    name_length = record[12]
    return record[36:36 + name_length - 1].decode()


def output_records(result):
    require(result.returncode == 0, result.stderr.decode())
    return split_bam(unpack_bgzf(result.stdout))


def fastq_payload(sequences):
    data = bytearray()
    for name, sequence in sequences.items():
        data += f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()
    return bytes(data)


def assert_fastq_bam_parity(sequences, options):
    records = [bam_record(name, sequence, flag=0x4) for name, sequence in sequences.items()]
    bam_result = run(["--bam-subset", *options], bgzf(bam_payload(records)))
    _, bam_records = output_records(bam_result)
    fastq_result = run(["--fastq-subset", *options], fastq_payload(sequences))
    require(fastq_result.returncode == 0, fastq_result.stderr.decode())
    bam_names = [record_name(record) for record in bam_records]
    fastq_names = [line[1:].decode() for line in fastq_result.stdout.splitlines()[::4]]
    require(bam_names == fastq_names, "FASTQ/BAM scoring parity failed")
    return bam_names


def explicit_args(path="-"):
    args = ["--bam-subset", "-x", "0", "-l", "18", "-y", "0.8", "-k", "10", "-d", "10"]
    if path != "-":
        args.append(str(path))
    return args


def test_basic_and_record_preservation(tmp):
    passing = [
        bam_record("mapped_p", "CCCTAACCCTAACCCTAA", ref_id=0, pos=10, cigar=((18, 0),), tags=b"NM\x69\x00\x00\x00\x00"),
        bam_record("reverse_q", "TTAGGGTTAGGGTTAGGG", flag=0x10, ref_id=0, pos=20, cigar=((18, 0),)),
        bam_record("secondary", "TTAGGGTTAGGGTTAGGG", flag=0x100),
        bam_record("supplementary", "CCCTAACCCTAACCCTAA", flag=0x800),
        bam_record("odd_ambiguous", "NTTAGGGTTAGGGTTAGGG", flag=0x4),
        bam_record("paired_first", "TTAGGGTTAGGGTTAGGG", flag=0x41),
        bam_record("paired_second", "CCCTAACCCTAACCCTAA", flag=0x81),
        bam_record("all_codes", "=ACMGRSVTWYHKDBNTTAGGGTTAGGGTTAGGG", flag=0x4),
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
    result = run(explicit_args(), input_bam)
    require(result.returncode == 0, result.stderr.decode())
    header, records = split_bam(unpack_bgzf(result.stdout))
    expected_header, _ = split_bam(input_payload)
    require(header == expected_header, "BAM header changed")
    require(records == passing, "passing BAM records were not preserved exactly")
    require(b"kept 8 of 10 records" in result.stderr, "summary count mismatch")
    require(b"skipped 1 record without SEQ" in result.stderr, "missing SEQ count mismatch")


def test_file_output_and_threads(tmp):
    records = []
    expected = []
    for index in range(700):
        if index % 4:
            record = bam_record(f"pass_{index}", "TTAGGG" * (3 + index % 7), flag=0x4)
            expected.append(record)
        else:
            record = bam_record(f"fail_{index}", "ACGT" * 20, flag=0x4)
        records.append(record)
    input_path = tmp / "many.bam"
    input_path.write_bytes(bgzf(bam_payload(records), chunk_size=113))

    single = run([*explicit_args(input_path), "-j", "1"])
    multi = run([*explicit_args(input_path), "-j", "8"])
    require(single.returncode == 0, single.stderr.decode())
    require(multi.returncode == 0, multi.stderr.decode())
    require(single.stdout == multi.stdout, "thread count changed BAM output")
    _, actual = split_bam(unpack_bgzf(single.stdout))
    require(actual == expected, "large-batch filtering mismatch")

    out_dir = tmp / "out"
    saved = run([*explicit_args(input_path), "-o", str(out_dir)])
    require(saved.returncode == 0, saved.stderr.decode())
    output_path = out_dir / "many_telomeric.bam"
    require(output_path.exists(), "BAM -o output name mismatch")
    require(output_path.read_bytes() == single.stdout, "file and stdout BAM differ")
    require(saved.stdout == b"", "BAM -o wrote binary data to stdout")


def test_header_only_and_missing_eof(tmp):
    payload = bam_payload([])
    result = run(explicit_args(), bgzf(payload, eof=False))
    require(result.returncode == 0, result.stderr.decode())
    require(b"missing the BGZF EOF marker" in result.stderr, "missing EOF warning absent")
    header, records = split_bam(unpack_bgzf(result.stdout))
    require(header == payload and records == [], "header-only BAM changed")


def test_header_variants_and_empty_results(tmp):
    passing = [
        bam_record("first", "TTAGGG" * 3, flag=0x4),
        bam_record("second", "CCCTAA" * 4, flag=0x4),
    ]
    payload = bam_payload(
        passing,
        header_text=b"",
        references=(("a", 1), ("long_reference_name", 2**31 - 1)),
    )
    result = run(explicit_args(), bgzf(payload))
    header, records = output_records(result)
    expected_header, _ = split_bam(payload)
    require(header == expected_header and records == passing, "BAM header variant changed")

    no_reference_payload = bam_payload(passing, header_text=b"@CO\tempty references\n", references=())
    no_reference = run(explicit_args(), bgzf(no_reference_payload))
    header, records = output_records(no_reference)
    expected_header, _ = split_bam(no_reference_payload)
    require(header == expected_header and records == passing, "reference-free BAM changed")

    failing = [bam_record(f"fail_{index}", "ACGT" * (index + 2), flag=0x4) for index in range(5)]
    fail_payload = bam_payload(failing)
    all_fail = run(explicit_args(), bgzf(fail_payload))
    header, records = output_records(all_fail)
    expected_header, _ = split_bam(fail_payload)
    require(header == expected_header and records == [], "all-fail BAM was not header-only")
    require(b"kept 0 of 5 records" in all_fail.stderr, "all-fail count mismatch")


def test_embedded_eof(tmp):
    records = [
        bam_record("first", "TTAGGG" * 3, flag=0x4),
        bam_record("second", "CCCTAA" * 3, flag=0x4),
    ]
    payload = bam_payload(records)
    split = len(payload) // 2
    input_bam = bgzf(payload[:split]) + bgzf(payload[split:])
    result = run(explicit_args(), input_bam)
    require(result.returncode == 0, result.stderr.decode())
    _, actual = split_bam(unpack_bgzf(result.stdout))
    require(actual == records, "embedded BGZF EOF interrupted input")


def test_cross_block_large_record(tmp):
    generator = random.Random(7)
    sequence = "".join(generator.choice("ACGT") for _ in range(70000)) + "TTAGGG" * 10
    qualities = [generator.randrange(94) for _ in sequence]
    record = bam_record("large", sequence, flag=0x4, qualities=qualities)
    payload = bam_payload([record])
    result = run(explicit_args(), bgzf(payload, chunk_size=60000))
    require(result.returncode == 0, result.stderr.decode())
    _, records = split_bam(unpack_bgzf(result.stdout))
    require(records == [record], "record spanning BGZF blocks changed")


def test_byte_bounded_batch(tmp):
    aux_count = 33 * 1024 * 1024
    huge_aux = b"ZZBC" + struct.pack("<i", aux_count) + bytes(aux_count)
    huge = bam_record("huge_aux", "TTAGGG" * 3, flag=0x4, tags=huge_aux)
    tail = bam_record("tail", "CCCTAA" * 3, flag=0x4)
    result = run([*explicit_args(), "-j", "2"], bgzf(bam_payload([huge, tail])))
    _, records = output_records(result)
    require(records == [huge, tail], "byte-bounded batch changed records")


def test_default_threshold_parity(tmp):
    sequences = {
        "short": "TTAGGG" * 9,
        "default_pass": "TTAGGG" * 10,
        "long": "CCCTAA" * 15,
        "fail": "ACGT" * 20,
    }
    names = assert_fastq_bam_parity(sequences, [])
    require(names == ["default_pass", "long"], "default threshold boundary failed")

    sequence = "TTAGGG" * 10
    crlf = f"@crlf\r\n{sequence}\r\n+\r\n{'I' * len(sequence)}\r\n".encode()
    result = run(["--fastq-subset"], crlf)
    require(result.returncode == 0, result.stderr.decode())
    require(result.stdout.startswith(b"@crlf\r\n"), "CRLF FASTQ filtering failed")


def test_exact_math_boundaries(tmp):
    lengths = {
        "one_repeat": "TTAGGG",
        "exact_12": "TTAGGG" * 2,
        "flanked_exact": "ACGT" + "CCCTAA" * 2 + "TGCA",
        "exact_18": "TTAGGG" * 3,
    }
    length_12 = ["-x", "0", "-l", "12", "-y", "1", "-k", "10", "-d", "10"]
    length_18 = ["-x", "0", "-l", "18", "-y", "1", "-k", "10", "-d", "10"]
    names = assert_fastq_bam_parity(lengths, length_12)
    require(names == ["exact_12", "flanked_exact", "exact_18"], "12 bp length boundary failed")
    names = assert_fastq_bam_parity(lengths, length_18)
    require(names == ["exact_18"], "18 bp length boundary failed")

    density = {"two_thirds": "TTAGGGAAAAAATTAGGG"}
    pass_options = ["-x", "0", "-l", "18", "-y", "0.666", "-k", "20", "-d", "10"]
    fail_options = ["-x", "0", "-l", "18", "-y", "0.667", "-k", "20", "-d", "10"]
    require(assert_fastq_bam_parity(density, pass_options) == ["two_thirds"], "density lower boundary failed")
    require(assert_fastq_bam_parity(density, fail_options) == [], "density upper boundary failed")

    plant = {
        "plant_pass": "TTTAGGG" * 3,
        "vertebrate_fail": "TTAGGG" * 4,
    }
    plant_options = ["-c", "CCCTAAA", "-x", "0", "-l", "21", "-y", "1"]
    require(assert_fastq_bam_parity(plant, plant_options) == ["plant_pass"], "custom canonical failed")


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
        ["-x", "0", "-l", "18", "-y", "0.8", "-k", "10", "-d", "10"],
        ["-x", "1", "-l", "42", "-y", "0.5", "-k", "50", "-d", "50"],
        ["-x", "0", "-l", "60", "-y", "1", "-k", "10", "-d", "10"],
    ]
    for options in option_sets:
        assert_fastq_bam_parity(sequences, options)


def test_cli_guards_and_cleanup(tmp):
    record = bam_record("pass", "TTAGGG" * 3, flag=0x4)
    input_bam = bgzf(bam_payload([record]))
    command = run([*explicit_args(), "--cmd"], input_bam)
    require(command.returncode == 0, command.stderr.decode())
    _, records = split_bam(unpack_bgzf(command.stdout))
    require(records == [record], "--cmd corrupted BAM stdout")
    require(b"--bam-subset" in command.stderr, "--cmd did not use stderr")

    conflict = run(["--bam-subset", "--fastq-subset"], input_bam)
    require(conflict.returncode != 0, "conflicting subset modes succeeded")

    input_path = tmp / "broken.bam"
    input_path.write_bytes(input_bam[:20])
    out_dir = tmp / "broken_out"
    failed = run([*explicit_args(input_path), "-o", str(out_dir)])
    require(failed.returncode != 0, "truncated file unexpectedly succeeded")
    require(not (out_dir / "broken_telomeric.bam").exists(), "partial BAM output was retained")

    if os.name != "nt" and os.geteuid() != 0:
        unreadable = tmp / "unreadable.bam"
        unreadable.write_bytes(input_bam)
        unreadable.chmod(0)
        try:
            failed = run(explicit_args(unreadable))
            require(failed.returncode != 0, "unreadable BAM unexpectedly succeeded")
            require(b"cannot open BAM input" in failed.stderr, "input open failure was unclear")
        finally:
            unreadable.chmod(0o600)

        unwritable = tmp / "unwritable"
        unwritable.mkdir()
        unwritable.chmod(0o500)
        try:
            failed = run([*explicit_args(input_path), "-o", str(unwritable)])
            require(failed.returncode != 0, "unwritable output unexpectedly succeeded")
            require(b"cannot write telomeric records" in failed.stderr, "output open failure was unclear")
        finally:
            unwritable.chmod(0o700)


def test_failures(tmp):
    record = bam_record("pass", "TTAGGG" * 3, flag=0x4)
    payload = bam_payload([record])
    valid = bytearray(bgzf(payload))
    header, records = split_bam(payload)
    record_offset = len(header)

    cases = {}
    gzip_compressor = zlib.compressobj(6, zlib.DEFLATED, 31)
    cases["not_bgzf"] = b"not bam"
    cases["plain_gzip"] = gzip_compressor.compress(payload) + gzip_compressor.flush()
    cases["bad_magic"] = bgzf(b"BAD\1" + payload[4:])
    cases["truncated"] = bytes(valid[:20])
    cases["empty"] = b""

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
        result = run(explicit_args(), data)
        require(result.returncode != 0, f"{name} unexpectedly succeeded")
        require(b"BAM subset failed" in result.stderr, f"{name} lacked a clear error")


def test_mutation_robustness(tmp):
    case_count = int(os.environ.get("BAM_MUTATION_CASES", "48"))
    generator = random.Random(91)
    records = [
        bam_record(f"record_{index}", "TTAGGG" * (3 + index % 5), flag=0x4)
        for index in range(8)
    ]
    payload = bam_payload(records)
    header, split_records = split_bam(payload)
    record_offset = len(header)

    for index in range(case_count):
        mode = index % 7
        if mode == 0:
            size = generator.randrange(0, 2048)
            data = bytes(generator.getrandbits(8) for _ in range(size))
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

        result = run(["--bam-subset", "-j", "1"], data)
        if result.returncode == 0:
            split_bam(unpack_bgzf(result.stdout))
        else:
            require(b"BAM subset failed" in result.stderr,
                    f"mutation {index} lacked a clear error")


def main():
    require(TELOSCOPE.exists(), f"Teloscope binary not found: {TELOSCOPE}")
    with tempfile.TemporaryDirectory(prefix="teloscope_bam_") as temp:
        tmp = pathlib.Path(temp)
        test_basic_and_record_preservation(tmp)
        test_file_output_and_threads(tmp)
        test_header_only_and_missing_eof(tmp)
        test_header_variants_and_empty_results(tmp)
        test_embedded_eof(tmp)
        test_cross_block_large_record(tmp)
        test_byte_bounded_batch(tmp)
        test_default_threshold_parity(tmp)
        test_exact_math_boundaries(tmp)
        test_randomized_fastq_bam_parity(tmp)
        test_cli_guards_and_cleanup(tmp)
        test_failures(tmp)
        test_mutation_robustness(tmp)
    print("PASS BAM subset integration")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL BAM subset integration: {error}", file=sys.stderr)
        raise
