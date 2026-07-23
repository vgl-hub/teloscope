#!/usr/bin/env python3

import gzip
import os
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))
MULTI_FASTA = ROOT / "testFiles" / "multi.fa"
PATH_GFA = ROOT / "testFiles" / "gfa_path_orient_pairs_small.gfa"
PATHLESS_GFA = ROOT / "testFiles" / "gfa_pathless_small.gfa"
SHARED_GFA = ROOT / "testFiles" / "gfa_single_seg_paths_small.gfa"
LARGE_GZIP_FASTA = ROOT / "testFiles" / "bTaeGut7_chr33_mat.fa.gz"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def decode(data):
    return data.decode("utf-8", errors="replace")


def run(args, stdin=None):
    try:
        return subprocess.run(
            [str(TELOSCOPE), *[str(arg) for arg in args]],
            input=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=60,
        )
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"Teloscope timed out: {' '.join(str(arg) for arg in args)}") from error


def require_success(result, context):
    require(
        result.returncode == 0,
        f"{context} failed with exit {result.returncode}\nstdout:\n{decode(result.stdout)}"
        f"\nstderr:\n{decode(result.stderr)}",
    )


def require_failure(result, context, expected=None):
    require(result.returncode != 0, f"{context} unexpectedly succeeded")
    if expected is not None:
        combined = decode(result.stdout) + decode(result.stderr)
        require(expected in combined, f"{context} lacked diagnostic {expected!r}:\n{combined}")


def parse_report_names(text):
    names = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[0].isdigit():
            names.append(fields[1])
    return names


def parse_report_table(text):
    header = None
    rows = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[:2] == ["pos", "header"]:
            header = fields
        elif len(fields) >= 2 and fields[0].isdigit():
            rows.append(fields)
    require(header is not None, "tabular report header is absent")
    for fields in rows:
        require(
            len(fields) == len(header),
            f"tabular report row changed schema ({len(fields)} fields, expected {len(header)}): {fields}",
        )
    return header, rows


def summary_value(text, label):
    prefix = f"{label}:\t"
    for line in text.splitlines():
        if line.startswith(prefix):
            return line[len(prefix):]
    raise AssertionError(f"summary field {label!r} was absent")


def only_report(out_dir):
    reports = sorted(out_dir.glob("*_report.tsv"))
    require(len(reports) == 1, f"expected one report in {out_dir}, found {reports}")
    return reports[0]


def assert_fasta_selection(result, out_dir, expected_names, input_count):
    require_success(result, "filtered FASTA run")
    expected_names = list(expected_names)
    stdout = decode(result.stdout)
    stderr = decode(result.stderr)
    report = only_report(out_dir).read_text(encoding="utf-8")

    parse_report_table(stdout)
    parse_report_table(report)
    require(parse_report_names(stdout) == expected_names, "stdout contains the wrong selected record order")
    require(parse_report_names(report) == expected_names, "report contains the wrong selected record order")
    require(
        f"Sequence filter: selected {len(expected_names)} of {input_count} paths." in stderr,
        "stderr selection count is missing or wrong",
    )
    require(summary_value(stdout, "Total paths") == str(len(expected_names)), "Total paths is wrong")
    require(summary_value(stdout, "Filter input paths") == str(input_count), "input count is wrong")
    require(
        summary_value(stdout, "Filter selected paths") == str(len(expected_names)),
        "selected count is wrong",
    )


def assert_excluded_absent(out_dir, excluded_names):
    outputs = sorted(
        path for path in out_dir.iterdir()
        if path.suffix in {".bed", ".bedgraph", ".tsv"}
    )
    require(outputs, f"no text outputs found in {out_dir}")
    for path in outputs:
        text = path.read_text(encoding="utf-8")
        for name in excluded_names:
            require(name not in text, f"excluded ID {name!r} leaked into {path.name}")


def run_fasta(fasta, out_dir, options):
    return run(["-f", fasta, "-o", out_dir, "-j", "1", *options])


def read_fasta(path):
    records = []
    header = None
    sequence = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header = line[1:]
            sequence = []
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path, records):
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for header, sequence in records:
            stream.write(f">{header}\n{sequence}\n")


def test_exact_id_and_all_text_outputs(tmp):
    out_dir = tmp / "exact_out"
    selector = tmp / "include.ids"
    selector.write_text("contig_t2t\n", encoding="utf-8")
    result = run_fasta(
        MULTI_FASTA,
        out_dir,
        ["--include-bed", selector, "-r", "-g", "-e", "-m", "-i"],
    )

    assert_fasta_selection(result, out_dir, ["contig_t2t"], 3)
    assert_excluded_absent(out_dir, ["contig_none", "contig_incomplete"])
    density = out_dir / "multi.fa_window_repeat_density.bedgraph"
    require(density.exists(), "repeat-density output was not generated")
    density_ids = {
        fields[0]
        for line in density.read_text(encoding="utf-8").splitlines()
        if len(fields := line.split("\t")) >= 4
    }
    require(density_ids == {"contig_t2t"}, "BEDgraph contains an unselected record")


def test_bed3_comments_directives_crlf_and_duplicates(tmp):
    out_dir = tmp / "bed3_out"
    selector = tmp / "include.bed"
    selector.write_bytes(
        b"# generated selector\r\n"
        b"\r\n"
        b"track name=selected\r\n"
        b"browser position contig_t2t:1-10\r\n"
        b"contig_t2t\t0\t0\r\n"
        b"contig_t2t\t0\t100\textra\r\n"
    )
    result = run_fasta(MULTI_FASTA, out_dir, ["--include-bed", selector])
    assert_fasta_selection(result, out_dir, ["contig_t2t"], 3)

    bom_selector = tmp / "bom.ids"
    bom_selector.write_bytes(b"\xef\xbb\xbfcontig_t2t\r\n")
    bom_out = tmp / "bom_out"
    bom_result = run_fasta(MULTI_FASTA, bom_out, ["--include-bed", bom_selector])
    assert_fasta_selection(bom_result, bom_out, ["contig_t2t"], 3)


def test_include_union_repetition_and_exclusion_precedence(tmp):
    include_t2t = tmp / "include_t2t.ids"
    include_incomplete = tmp / "include_incomplete.bed"
    exclude_t2t = tmp / "exclude_t2t.ids"
    include_t2t.write_text("contig_t2t\ncontig_t2t\n", encoding="utf-8")
    include_incomplete.write_text("contig_incomplete\t0\t1\n", encoding="utf-8")
    exclude_t2t.write_text("contig_t2t\n", encoding="utf-8")

    union_out = tmp / "union_out"
    union_result = run_fasta(
        MULTI_FASTA,
        union_out,
        [
            "--include-bed", include_t2t,
            "--include-bed", include_incomplete,
            "--include-bed", include_t2t,
            "--include-prefix", " contig_inc , contig_t2 ",
            "--include-prefix", "contig_inc",
            "--exclude-prefix", " contig_none ",
        ],
    )
    assert_fasta_selection(union_result, union_out, ["contig_t2t", "contig_incomplete"], 3)

    subtract_out = tmp / "subtract_out"
    subtract_result = run_fasta(
        MULTI_FASTA,
        subtract_out,
        [
            "--include-bed", include_t2t,
            "--include-bed", include_incomplete,
            "--exclude-bed", exclude_t2t,
        ],
    )
    assert_fasta_selection(subtract_result, subtract_out, ["contig_incomplete"], 3)

    exclude_only_out = tmp / "exclude_only_out"
    exclude_only_result = run_fasta(
        MULTI_FASTA,
        exclude_only_out,
        ["--exclude-prefix", "contig_none"],
    )
    assert_fasta_selection(
        exclude_only_result,
        exclude_only_out,
        ["contig_t2t", "contig_incomplete"],
        3,
    )


def test_prefixes_are_trimmed_literal_and_case_sensitive(tmp):
    source_records = read_fasta(MULTI_FASTA)
    sequences = [sequence for _, sequence in source_records]
    fasta = tmp / "prefixes.fa"
    write_fasta(
        fasta,
        [
            ("CM.alpha chromosome one", sequences[0]),
            ("CMXalpha chromosome two", sequences[1]),
            ("cm.alpha lower case", sequences[2]),
            ("star*one literal wildcard", sequences[0]),
        ],
    )

    dot_out = tmp / "literal_dot_out"
    dot_result = run_fasta(fasta, dot_out, ["--include-prefix", "CM."])
    assert_fasta_selection(dot_result, dot_out, ["CM.alpha"], 4)

    case_out = tmp / "case_out"
    case_result = run_fasta(fasta, case_out, ["--include-prefix", "cm."])
    assert_fasta_selection(case_result, case_out, ["cm.alpha"], 4)

    grouped_out = tmp / "grouped_out"
    grouped_result = run_fasta(
        fasta,
        grouped_out,
        ["--include-prefix", " star* , CM. ", "--include-prefix", "CM."],
    )
    assert_fasta_selection(grouped_result, grouped_out, ["CM.alpha", "star*one"], 4)


def test_fasta_primary_ids_with_crlf_and_tab_descriptions(tmp):
    source_records = read_fasta(MULTI_FASTA)
    sequences = [sequence.encode() for _, sequence in source_records]

    crlf_fasta = tmp / "crlf.fa"
    crlf_fasta.write_bytes(
        b">CMCRLF\r\n" + sequences[0] + b"\r\n"
        b">NWCRLF\r\n" + sequences[1] + b"\r\n"
    )
    crlf_selector = _write(tmp / "crlf.ids", "CMCRLF\n")
    crlf_out = tmp / "crlf_out"
    crlf_result = run_fasta(crlf_fasta, crlf_out, ["--include-bed", crlf_selector])
    assert_fasta_selection(crlf_result, crlf_out, ["CMCRLF"], 2)
    expected_length = str(len(sequences[0]))
    require(
        summary_value(decode(crlf_result.stdout), "Scaffold N50") == expected_length,
        "CRLF changed the selected scaffold length",
    )
    require(
        summary_value(decode(crlf_result.stdout), "Contig N50") == expected_length,
        "CRLF changed the selected contig length",
    )

    terminal_records = [
        (b"CM_terminal_N", b"ACGT" * 400 + b"N"),
        (b"CM_terminal_X", b"TGCA" * 400 + b"X"),
    ]
    lf_terminal = tmp / "terminal_lf.fa"
    crlf_terminal = tmp / "terminal_crlf.fa"
    lf_terminal.write_bytes(
        b"".join(b">" + name + b"\n" + sequence + b"\n" for name, sequence in terminal_records)
    )
    crlf_terminal.write_bytes(
        b"".join(b">" + name + b"\r\n" + sequence + b"\r\n" for name, sequence in terminal_records)
    )
    lf_out = tmp / "terminal_lf_out"
    crlf_terminal_out = tmp / "terminal_crlf_out"
    lf_result = run_fasta(lf_terminal, lf_out, ["--include-prefix", "CM_terminal_"])
    crlf_terminal_result = run_fasta(
        crlf_terminal,
        crlf_terminal_out,
        ["--include-prefix", "CM_terminal_"],
    )
    terminal_names = ["CM_terminal_N", "CM_terminal_X"]
    assert_fasta_selection(lf_result, lf_out, terminal_names, 2)
    assert_fasta_selection(crlf_terminal_result, crlf_terminal_out, terminal_names, 2)
    for label in ("Total gaps", "Scaffold N50", "Contig N50"):
        require(
            summary_value(decode(crlf_terminal_result.stdout), label)
            == summary_value(decode(lf_result.stdout), label),
            f"CRLF changed {label} for a terminal N/X record",
        )
    require(
        parse_report_table(only_report(crlf_terminal_out).read_text(encoding="utf-8"))
        == parse_report_table(only_report(lf_out).read_text(encoding="utf-8")),
        "CRLF changed report classification or component counts for terminal N/X records",
    )
    lf_gaps = next(lf_out.glob("*_gaps.bed")).read_text(encoding="utf-8")
    crlf_gaps = next(crlf_terminal_out.glob("*_gaps.bed")).read_text(encoding="utf-8")
    require(crlf_gaps == lf_gaps, "CRLF changed terminal N/X gap coordinates")

    tab_fasta = tmp / "tab_header.fa"
    tab_fasta.write_bytes(
        b">CMTAB\tdescription with tabs\n" + sequences[0] + b"\n"
        b">NWTAB\tdescription with tabs\n" + sequences[1] + b"\n"
    )
    tab_selector = _write(tmp / "tab.ids", "CMTAB\n")
    tab_out = tmp / "tab_out"
    tab_result = run_fasta(tab_fasta, tab_out, ["--include-bed", tab_selector])
    assert_fasta_selection(tab_result, tab_out, ["CMTAB"], 2)
    tab_report = only_report(tab_out).read_text(encoding="utf-8")
    header, rows = parse_report_table(tab_report)
    require(len(rows) == 1, "tab-description FASTA should produce one selected report row")
    require(rows[0][header.index("header")] == "CMTAB", "FASTA description leaked into the report ID")

    duplicate_fasta = tmp / "duplicate_primary_ids.fa"
    duplicate_fasta.write_bytes(
        b">dup first description\n" + sequences[0] + b"\n"
        b">dup\tsecond description\n" + sequences[1] + b"\n"
    )
    duplicate_selector = _write(tmp / "duplicate.ids", "dup\n")
    duplicate_result = run_fasta(
        duplicate_fasta,
        tmp / "duplicate_out",
        ["--include-bed", duplicate_selector],
    )
    require_failure(duplicate_result, "duplicate normalized FASTA primary IDs")
    duplicate_diagnostic = (decode(duplicate_result.stdout) + decode(duplicate_result.stderr)).lower()
    require("dup" in duplicate_diagnostic, "duplicate-ID diagnostic did not name the colliding ID")
    require(
        "duplicate" in duplicate_diagnostic or "already exists" in duplicate_diagnostic,
        "duplicate-ID failure lacked a clear diagnostic",
    )


def test_wrapped_gap_lf_crlf_parity(tmp):
    lf_fasta = tmp / "wrapped_lf.fa"
    crlf_fasta = tmp / "wrapped_crlf.fa"
    lf_fasta.write_bytes(b">CMgap\nCCCTAANNN\nNNNTTAGGG\n")
    crlf_fasta.write_bytes(b">CMgap\r\nCCCTAANNN\r\nNNNTTAGGG\r\n")

    lf_out = tmp / "wrapped_lf_out"
    crlf_out = tmp / "wrapped_crlf_out"
    lf_result = run_fasta(lf_fasta, lf_out, ["--include-prefix", "CMgap"])
    crlf_result = run_fasta(crlf_fasta, crlf_out, ["--include-prefix", "CMgap"])
    assert_fasta_selection(lf_result, lf_out, ["CMgap"], 1)
    assert_fasta_selection(crlf_result, crlf_out, ["CMgap"], 1)

    for label in ("Total gaps", "Scaffold N50", "Contig N50"):
        require(
            summary_value(decode(crlf_result.stdout), label)
            == summary_value(decode(lf_result.stdout), label),
            f"wrapped CRLF input changed {label}",
        )
    require(summary_value(decode(crlf_result.stdout), "Total gaps") == "1", "wrapped N run split into multiple gaps")
    require(summary_value(decode(crlf_result.stdout), "Scaffold N50") == "18", "wrapped CRLF scaffold length is wrong")
    require(summary_value(decode(crlf_result.stdout), "Contig N50") == "6", "wrapped CRLF contig N50 is wrong")

    lf_report = parse_report_table(only_report(lf_out).read_text(encoding="utf-8"))
    crlf_report = parse_report_table(only_report(crlf_out).read_text(encoding="utf-8"))
    require(crlf_report == lf_report, "wrapped CRLF input changed the report row")
    lf_gaps = next(lf_out.glob("*_gaps.bed")).read_text(encoding="utf-8")
    crlf_gaps = next(crlf_out.glob("*_gaps.bed")).read_text(encoding="utf-8")
    require(lf_gaps == "CMgap\t6\t12\n", f"wrapped LF gap coordinates are wrong: {lf_gaps!r}")
    require(crlf_gaps == lf_gaps, "wrapped CRLF input changed gap coordinates")


def test_plain_gzip_and_stdin_parity(tmp):
    payload = MULTI_FASTA.read_bytes()
    compressed = tmp / "multi.fa.gz"
    with gzip.open(compressed, "wb") as stream:
        stream.write(payload)

    plain_out = tmp / "plain_out"
    plain = run_fasta(MULTI_FASTA, plain_out, ["--include-prefix", "contig_t2t"])
    assert_fasta_selection(plain, plain_out, ["contig_t2t"], 3)

    gzip_out = tmp / "gzip_out"
    zipped = run_fasta(compressed, gzip_out, ["--include-prefix", "contig_t2t"])
    assert_fasta_selection(zipped, gzip_out, ["contig_t2t"], 3)

    stdin_out = tmp / "stdin_out"
    piped = run(
        ["-o", stdin_out, "-j", "1", "--include-prefix", "contig_t2t"],
        stdin=payload,
    )
    assert_fasta_selection(piped, stdin_out, ["contig_t2t"], 3)
    require(only_report(stdin_out).name == "stdin_report.tsv", "stdin output basename changed")

    for result in (zipped, piped):
        require(
            summary_value(decode(result.stdout), "Scaffold N50")
            == summary_value(decode(plain.stdout), "Scaffold N50"),
            "plain/gzip/stdin filtering changed selected assembly statistics",
        )


def test_large_gzip_fasta_is_not_truncated_or_stalled(tmp):
    out_dir = tmp / "large_gzip_out"
    result = run_fasta(
        LARGE_GZIP_FASTA,
        out_dir,
        ["--include-prefix", "chr33_mat", "-u"],
    )
    assert_fasta_selection(result, out_dir, ["chr33_mat"], 1)
    require(
        summary_value(decode(result.stdout), "Scaffold N50") == "4246341",
        "large gzip FASTA was truncated",
    )


def test_selector_validation_failures(tmp):
    malformed = {
        "empty": (b"", "contains no sequence IDs"),
        "directives_only": (
            b"# no IDs\ntrack name=none\nbrowser position contig_t2t:1-2\n",
            "contains no sequence IDs",
        ),
        "two_columns": (b"contig_t2t\t0\n", "must contain either one ID column or at least three BED columns"),
        "negative_start": (b"contig_t2t\t-1\t2\n", "invalid BED start/end coordinates"),
        "non_numeric": (b"contig_t2t\tzero\t2\n", "invalid BED start/end coordinates"),
        "overflow": (b"contig_t2t\t18446744073709551616\t18446744073709551616\n", "invalid BED start/end coordinates"),
        "reversed": (b"contig_t2t\t3\t2\n", "invalid BED start/end coordinates"),
    }
    for name, (contents, diagnostic) in malformed.items():
        selector = tmp / f"{name}.bed"
        selector.write_bytes(contents)
        result = run_fasta(MULTI_FASTA, tmp / f"{name}_out", ["--include-bed", selector])
        require_failure(result, name, diagnostic)

    missing = tmp / "does_not_exist.bed"
    result = run_fasta(MULTI_FASTA, tmp / "missing_out", ["--include-bed", missing])
    require_failure(result, "missing selector")
    require(result.stdout == b"", "missing selector diagnostic polluted stdout")
    require("--include-bed" in decode(result.stderr), "missing selector diagnostic omitted the option name")
    require("does not exist" in decode(result.stderr), "missing selector lacked a stderr diagnostic")

    unmatched_id = tmp / "unmatched.ids"
    unmatched_id.write_text("absent_record\n", encoding="utf-8")
    result = run_fasta(MULTI_FASTA, tmp / "unmatched_id_out", ["--include-bed", unmatched_id])
    require_failure(result, "unmatched exact ID", "matched no input paths")

    result = run_fasta(MULTI_FASTA, tmp / "unmatched_exclude_out", ["--exclude-bed", unmatched_id])
    require_failure(result, "unmatched exact exclusion ID", "matched no input paths")

    result = run_fasta(
        MULTI_FASTA,
        tmp / "unmatched_prefix_out",
        ["--include-prefix", "absent_prefix"],
    )
    require_failure(result, "unmatched prefix", "matched no input paths")

    result = run_fasta(MULTI_FASTA, tmp / "empty_selection_out", ["--exclude-prefix", "contig_"])
    require_failure(result, "all records excluded", "excluded all input paths")

    overlap = _write(tmp / "overlap.ids", "contig_t2t\n")
    result = run_fasta(
        MULTI_FASTA,
        tmp / "overlap_out",
        ["--include-bed", overlap, "--exclude-bed", overlap],
    )
    require_failure(result, "include/exclude overlap", "excluded all input paths")

    for index, value in enumerate(("", ",contig", "contig,,none", "contig,  ,none", "contig,")):
        result = run_fasta(
            MULTI_FASTA,
            tmp / f"empty_prefix_{index}_out",
            ["--include-prefix", value],
        )
        require_failure(result, f"empty prefix token {value!r}", "contains an empty prefix")


def gfa_lines(path):
    return [line for line in path.read_text(encoding="utf-8").splitlines() if line]


def gfa_telomere_nodes(lines):
    return {
        fields[1]
        for line in lines
        if len(fields := line.split("\t")) >= 2
        and fields[0] == "S"
        and fields[1].startswith("telomere_")
    }


def assert_original_gfa_records_preserved(input_path, output_lines, exact_types):
    original = gfa_lines(input_path)
    output_set = set(output_lines)
    for line in original:
        if line.split("\t", 1)[0] in exact_types:
            require(line in output_set, f"original GFA record was not preserved: {line}")


def assert_colors_match_nodes(out_dir, nodes):
    colors_files = sorted(out_dir.glob("*.telo.annotated.colors.csv"))
    require(len(colors_files) == 1, "expected one GFA colors file")
    rows = colors_files[0].read_text(encoding="utf-8").splitlines()
    require(rows and rows[0] == "node\tcolor", "GFA colors header is wrong")
    colored_nodes = {row.split("\t", 1)[0] for row in rows[1:] if row}
    require(colored_nodes == nodes, "colors file and annotated telomere nodes disagree")


def test_gfa_path_filtering_and_preservation(tmp):
    out_dir = tmp / "gfa_path_out"
    result = run_fasta(
        PATH_GFA,
        out_dir,
        ["--include-bed", _write(tmp / "paths.ids", "path_pp\n"), "-x", "0", "-l", "60"],
    )
    require_success(result, "path-aware GFA filter")
    require("selected 1 of 2 paths" in decode(result.stderr), "path-aware GFA count is wrong")

    outputs = sorted(out_dir.glob("*.telo.annotated.gfa"))
    require(len(outputs) == 1, "expected one annotated GFA")
    lines = gfa_lines(outputs[0])
    nodes = gfa_telomere_nodes(lines)
    require(len(nodes) == 2, f"selected path should create two telomere nodes, found {nodes}")
    require(all("seg_fp" in node or "seg_lp" in node for node in nodes), "excluded path was annotated")
    require(not any("seg_fn" in node or "seg_ln" in node for node in nodes), "excluded path was annotated")
    assert_original_gfa_records_preserved(PATH_GFA, lines, {"H", "S", "L"})
    require(
        "P\tpath_nn\tseg_fn-,seg_ln-\t*" in lines,
        "excluded path topology was changed",
    )
    require({line.split("\t")[1] for line in lines if line.startswith("P\t")} == {"path_pp", "path_nn"},
            "an original GFA path disappeared")
    assert_colors_match_nodes(out_dir, nodes)


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def test_gfa_pathless_filtering_and_preservation(tmp):
    uppercase_gfa = tmp / "upper.GFA"
    uppercase_gfa.write_bytes(PATHLESS_GFA.read_bytes())
    out_dir = tmp / "gfa_pathless_out"
    result = run_fasta(
        uppercase_gfa,
        out_dir,
        ["--include-prefix", "seg_t2t", "-x", "0", "-l", "60"],
    )
    require_success(result, "pathless GFA filter")
    require("selected 1 of 4 segments" in decode(result.stderr), "pathless GFA count is wrong")

    outputs = sorted(out_dir.glob("*.telo.annotated.gfa"))
    require(len(outputs) == 1, "expected one pathless annotated GFA")
    lines = gfa_lines(outputs[0])
    nodes = gfa_telomere_nodes(lines)
    require(len(nodes) == 2, f"selected pathless segment should create two telomere nodes, found {nodes}")
    require(all("seg_t2t" in node for node in nodes), "an excluded pathless segment was annotated")
    for excluded in ("seg_p", "seg_q", "seg_none"):
        require(not any(f"telomere_{excluded}" in line for line in lines), f"{excluded} was annotated")
    assert_original_gfa_records_preserved(PATHLESS_GFA, lines, {"H", "S", "L"})
    assert_colors_match_nodes(out_dir, nodes)


def test_gfa_shared_terminal_selection_is_orientation_specific(tmp):
    expected_by_path = {
        "path_plus": {
            "telomere_seg_shared+_start",
            "telomere_seg_shared+_end",
        },
        "path_minus": {
            "telomere_seg_shared-_start",
            "telomere_seg_shared-_end",
        },
    }
    observed_by_path = {}
    for selected_path, expected_nodes in expected_by_path.items():
        out_dir = tmp / f"shared_out_{selected_path}"
        result = run_fasta(
            SHARED_GFA,
            out_dir,
            ["--include-prefix", selected_path, "-x", "0", "-l", "60", "-j", "8"],
        )
        require_success(result, "shared-terminal GFA filter")
        require("selected 1 of 2 paths" in decode(result.stderr), "shared-terminal GFA count is wrong")
        outputs = sorted(out_dir.glob("*.telo.annotated.gfa"))
        require(len(outputs) == 1, "expected one shared-terminal annotated GFA")
        lines = gfa_lines(outputs[0])
        nodes = gfa_telomere_nodes(lines)
        require(nodes == expected_nodes, f"wrong shared-terminal annotations: {nodes}")
        excluded_nodes = set().union(*expected_by_path.values()) - expected_nodes
        require(not any(node in line for node in excluded_nodes for line in lines), "excluded orientation was annotated")
        path_names = {line.split("\t")[1] for line in lines if line.startswith("P\t")}
        require(path_names == {"path_plus", "path_minus"}, "shared-segment path was lost")
        telomere_links = {
            line for line in lines
            if line.startswith("L\t") and any(node in line.split("\t") for node in nodes)
        }
        require(len(telomere_links) == 2, f"expected two shared-terminal links, found {telomere_links}")
        assert_colors_match_nodes(out_dir, nodes)
        observed_by_path[selected_path] = nodes
    require(
        observed_by_path["path_plus"].isdisjoint(observed_by_path["path_minus"]),
        "selected and excluded shared-segment orientations produced the same annotations",
    )


def test_unsupported_gfa_records_are_rejected(tmp):
    cases = {
        "gfa2": (
            "fixture.gfa2",
            "H\tVN:Z:2.0\nS\tsegA\t12\tTTAGGGTTAGGG\nO\tpathA\tsegA+\n",
            "do not support GFA2",
        ),
        "gfa2_header": (
            "header_disguised.gfa",
            "H\tVN:Z:2.0\nS\tsegA\t12\tTTAGGGTTAGGG\nO\tpathA\tsegA+\n",
            "do not support GFA2 at line",
        ),
        "gfa2_o_record": (
            "o_record.gfa",
            "H\tVN:Z:1.0\nS\tsegA\tTTAGGGTTAGGG\nO\tpathA\tsegA+\n",
            "do not support GFA2 record type 'O'",
        ),
        "walk": (
            "walk.gfa",
            "H\tVN:Z:1.1\nS\tsegA\tTTAGGGTTAGGG\nW\tsample\t0\tchr1\t0\t12\t>segA\n",
            "do not support GFA1 W walks",
        ),
        "containment": (
            "containment.gfa",
            "H\tVN:Z:1.0\nS\tsegA\tTTAGGGTTAGGG\nS\tsegB\tTTAGGG\nC\tsegA\t+\tsegB\t+\t0\t6M\n",
            "do not support GFA1 C containment",
        ),
        "unknown": (
            "unknown.gfa",
            "H\tVN:Z:1.0\nS\tsegA\tTTAGGGTTAGGG\nZ\tunsupported\n",
            "do not support GFA record type 'Z'",
        ),
    }
    for name, (filename, contents, diagnostic) in cases.items():
        fixture = _write(tmp / filename, contents)
        out_dir = tmp / f"{name}_out"
        result = run_fasta(fixture, out_dir, ["--include-prefix", "seg"])
        require_failure(result, f"unsupported {name} GFA", diagnostic)
        require(not list(out_dir.glob("*.telo.annotated.gfa")), f"unsupported {name} GFA was written")
        require(not list(out_dir.glob("*.telo.annotated.colors.csv")), f"unsupported {name} colors were written")


def test_fasta_in_gfa_named_directory_is_not_misclassified(tmp):
    misleading_dir = tmp / "inputs.gfa.archive"
    misleading_dir.mkdir()
    fasta = misleading_dir / "assembly.fa"
    fasta.write_bytes(MULTI_FASTA.read_bytes())
    out_dir = tmp / "misleading_directory_out"
    result = run_fasta(fasta, out_dir, ["--include-prefix", "contig_t2t"])
    assert_fasta_selection(result, out_dir, ["contig_t2t"], 3)
    require(not list(out_dir.glob("*.telo.annotated.gfa")), "FASTA path was misclassified as GFA")


def test_cli_surface_and_read_subset_guards(tmp):
    help_result = run(["--help"])
    require_success(help_result, "--help")
    help_text = decode(help_result.stdout)
    for option in ("--include-bed", "--exclude-bed", "--include-prefix", "--exclude-prefix"):
        require(option in help_text, f"{option} is absent from --help")

    version_result = run(["--version"])
    require_success(version_result, "--version")
    require("Teloscope v0.1.5" in decode(version_result.stdout), "binary version changed before release")

    for mode in ("--fastq-subset", "--bam-subset"):
        result = run([mode, "--include-prefix", "contig"], stdin=b"")
        require_failure(result, f"{mode} filter guard", "cannot be used in read subset mode")

    assembly_fastq = tmp / "assembly_mode.fastq"
    assembly_fastq.write_text(
        "@read1\nTTAGGGTTAGGG\n+\nIIIIIIIIIIII\n",
        encoding="utf-8",
    )
    fastq_out = tmp / "assembly_fastq_out"
    assembly_fastq_result = run_fasta(
        assembly_fastq,
        fastq_out,
        ["--include-prefix", "read"],
    )
    require_failure(
        assembly_fastq_result,
        "filtered FASTQ in assembly mode",
        "require FASTA input or a recognized GFA file",
    )
    require(not list(fastq_out.glob("*_report.tsv")), "filtered assembly-mode FASTQ produced a report")

    for option in ("--include-bed", "--exclude-bed", "--include-prefix", "--exclude-prefix"):
        missing_argument = run(["-f", MULTI_FASTA, option])
        require_failure(
            missing_argument,
            f"missing {option} argument",
            f"Option {option} is missing a required argument",
        )


TESTS = (
    test_exact_id_and_all_text_outputs,
    test_bed3_comments_directives_crlf_and_duplicates,
    test_include_union_repetition_and_exclusion_precedence,
    test_prefixes_are_trimmed_literal_and_case_sensitive,
    test_fasta_primary_ids_with_crlf_and_tab_descriptions,
    test_wrapped_gap_lf_crlf_parity,
    test_plain_gzip_and_stdin_parity,
    test_large_gzip_fasta_is_not_truncated_or_stalled,
    test_selector_validation_failures,
    test_gfa_path_filtering_and_preservation,
    test_gfa_pathless_filtering_and_preservation,
    test_gfa_shared_terminal_selection_is_orientation_specific,
    test_unsupported_gfa_records_are_rejected,
    test_fasta_in_gfa_named_directory_is_not_misclassified,
    test_cli_surface_and_read_subset_guards,
)


def main():
    require(TELOSCOPE.exists(), f"Teloscope binary not found: {TELOSCOPE}")
    require(MULTI_FASTA.exists(), f"FASTA fixture not found: {MULTI_FASTA}")
    require(SHARED_GFA.exists(), f"shared GFA fixture not found: {SHARED_GFA}")
    require(LARGE_GZIP_FASTA.exists(), f"large gzip FASTA fixture not found: {LARGE_GZIP_FASTA}")
    with tempfile.TemporaryDirectory(prefix="teloscope_sequence_filters_") as temp:
        temp_root = pathlib.Path(temp)
        for test in TESTS:
            case_dir = temp_root / test.__name__
            case_dir.mkdir()
            try:
                test(case_dir)
            except Exception as error:
                raise AssertionError(f"{test.__name__}: {error}") from error
    print(f"PASS sequence filter integration ({len(TESTS)} groups)")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL sequence filter integration: {error}", file=sys.stderr)
        raise
