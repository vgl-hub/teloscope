#!/usr/bin/env python3

import argparse
import gzip
import json
import os
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCES = (
    ROOT / "src" / "bam.cpp",
    ROOT / "src" / "bgzf.cpp",
    ROOT / "src" / "read-filter.cpp",
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-dir", required=True, type=pathlib.Path)
    return parser.parse_args()


def canonical(path):
    return pathlib.Path(path).resolve()


def main():
    args = parse_args()
    object_dir = args.object_dir.resolve()
    show_branches = os.environ.get("COVERAGE_BRANCH_DETAILS") == "1"
    for source in SOURCES:
        if not source.exists():
            raise RuntimeError(f"coverage source not found: {source}")

    with tempfile.TemporaryDirectory(prefix="teloscope_gcov_") as temp:
        result = subprocess.run(
            ["gcov", "-j", "-b", "-o", str(object_dir), *map(str, SOURCES)],
            cwd=temp,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)

        reports = {}
        for report_path in pathlib.Path(temp).glob("*.gcov.json.gz"):
            with gzip.open(report_path, "rt", encoding="utf-8") as report_file:
                report = json.load(report_file)
            for source_report in report.get("files", []):
                reports[canonical(source_report["file"])] = source_report

    failed = False
    for source in SOURCES:
        source_report = reports.get(canonical(source))
        if source_report is None:
            print(f"FAIL {source.relative_to(ROOT)}: gcov produced no report", file=sys.stderr)
            failed = True
            continue

        executable = []
        branch_total = 0
        branch_covered = 0
        for line in source_report.get("lines", []):
            count = int(line.get("count", 0))
            if count > 0 or line.get("unexecuted_block", False):
                executable.append(line)
            branches = line.get("branches", [])
            branch_total += len(branches)
            branch_covered += sum(int(branch.get("count", 0)) > 0 for branch in branches)

        if not executable:
            print(f"FAIL {source.relative_to(ROOT)}: no executable lines", file=sys.stderr)
            failed = True
            continue

        missed = [
            int(line["line_number"])
            for line in executable
            if int(line.get("count", 0)) == 0
        ]
        covered = len(executable) - len(missed)
        branch_text = (
            f", raw gcov branches {branch_covered}/{branch_total} informational"
            if show_branches and branch_total
            else ""
        )
        print(f"{source.relative_to(ROOT)}: lines {covered}/{len(executable)}{branch_text}")
        if missed:
            print("  missed executable lines: " + ", ".join(map(str, missed)), file=sys.stderr)
            failed = True

    if failed:
        raise RuntimeError("BAM hardening coverage is below 100% executable lines")
    print("PASS BAM hardening coverage: 100% executable lines")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL BAM hardening coverage: {error}", file=sys.stderr)
        sys.exit(1)
