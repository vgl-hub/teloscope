#!/usr/bin/env python3
"""Check every fixture against the intent declared for it in the synthetic manifest.
A disagreement means the fixture, the expectation, or the code is wrong; fix it
deliberately rather than editing the manifest to match.

Usage:
    python3 scripts/test_synthetic_intent.py [--only GLOB]
    TELOSCOPE=/path/to/other/teloscope python3 scripts/test_synthetic_intent.py
"""

import argparse
import fnmatch
import os
import pathlib
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))
MANIFEST = ROOT / "testFiles" / "synthetic" / "manifest.tsv"

# Manifest column -> report column; expect_its is checked only when the report has an its column.
CHECKED = [
    ("expect_type", "type"),
    ("expect_anomaly", "anomaly"),
    ("expect_granular", "granular"),
    ("expect_telomeres", "telomeres"),
    ("expect_labels", "labels"),
    ("expect_gaps", "gaps"),
    ("expect_its", "its"),
]

BED_FIELDS = [
    "chrom", "start", "end", "teloLen", "teloLabel", "closestEnd", "fwdCan", "revCan",
    "fwdNonCan", "revNonCan", "chrSize", "teloType",
]
INT_FIELDS = {"start", "end", "teloLen", "fwdCan", "revCan", "fwdNonCan", "revNonCan", "chrSize"}


def read_block_bed(path):
    """Terminal or interstitial BED -> list of dicts, in file order."""
    rows = []
    if not path.exists():
        return rows
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip() or s.lstrip().startswith(("#", "track", "browser")):
                continue
            f = s.split("\t")
            row = dict(zip(BED_FIELDS, f))
            for k in INT_FIELDS:
                row[k] = int(row[k])
            rows.append(row)
    return rows


def load_manifest(path):
    rows, header = [], None
    for raw in path.read_text().splitlines():
        if raw.startswith("#") or not raw.strip():
            continue
        fields = raw.split("\t")
        if header is None:
            header = fields
            continue
        rows.append(dict(zip(header, fields)))
    if header is None:
        raise SystemExit(f"{path}: no header row")
    return rows


def parse_report(text):
    """Read the per-path table by column name; the column set differs between scan modes."""
    columns, out = None, {}
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[:2] == ["pos", "header"]:
            columns = fields
            continue
        if columns and fields and fields[0].isdigit():
            row = dict(zip(columns, fields))
            out[row["header"]] = row
    return out


def run_fixture(fixture_path, flags, outdir):
    args = [str(TELOSCOPE), "-f", str(ROOT / "testFiles" / fixture_path), "-o", str(outdir)]
    if flags and flags != "-":
        args += shlex.split(flags)
    return subprocess.run(args, capture_output=True, timeout=120, check=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="*", help="glob over fixture ids")
    args = ap.parse_args()

    if not TELOSCOPE.exists():
        raise SystemExit(f"teloscope binary not found at {TELOSCOPE}; run `make head`")
    if not MANIFEST.exists():
        raise SystemExit(f"{MANIFEST} missing; run testFiles/generate_synthetic.sh")

    all_rows = load_manifest(MANIFEST)
    if not all_rows:
        raise SystemExit(f"{MANIFEST} has no rows; run testFiles/generate_synthetic.sh")
    rows = [r for r in all_rows if fnmatch.fnmatch(r["id"], args.only)]
    if not rows:
        raise SystemExit(f"no fixtures matched --only {args.only!r}; refusing to report PASS")
    # One run per (id, flags); a multi-record fixture checks every record from that run.
    by_run = {}
    for row in rows:
        by_run.setdefault((row["id"], row["path"], row["flags"]), []).append(row)

    failures, checks = [], 0

    for (fid, path, flags), group in sorted(by_run.items()):
        outdir = tempfile.mkdtemp(prefix="telo_intent_")
        try:
            proc = run_fixture(path, flags, outdir)
            if proc.returncode != 0:
                failures.append(
                    f"{fid} [{flags}]: teloscope exited {proc.returncode}\n"
                    f"    {proc.stderr.decode('utf-8', 'replace').strip()[:400]}")
                continue
            report = parse_report(proc.stdout.decode("utf-8", "replace"))
            terminal_bed = {}
            bed_path = pathlib.Path(outdir) / f"{pathlib.Path(path).name}_terminal_telomeres.bed"
            for b in read_block_bed(bed_path):
                terminal_bed.setdefault(b["chrom"], []).append(b)

            for row in group:
                scaffold = row["scaffold"]
                if scaffold not in report:
                    failures.append(
                        f"{fid} [{flags}]: scaffold {scaffold!r} absent from the report; "
                        f"saw {sorted(report)}")
                    continue
                actual = report[scaffold]
                for key, column in CHECKED:
                    want = row[key]
                    if want == "-":
                        continue
                    if want == "?":
                        failures.append(
                            f"{fid} [{flags}] {scaffold}: {column} is '?', an unstated "
                            f"expectation")
                        continue
                    if column not in actual:
                        # `its` is absent under ultra-fast; that is expected, not a failure.
                        if column == "its":
                            continue
                        failures.append(
                            f"{fid} [{flags}] {scaffold}: report has no column {column!r}")
                        continue
                    checks += 1
                    got = actual[column]
                    if got != want:
                        failures.append(
                            f"{fid} [{flags}] {scaffold}: {column} expected {want!r}, "
                            f"got {got!r}\n    intent: {row['intent']}")

                want_telolen = row.get("expect_telolen", "-")
                if want_telolen == "-":
                    continue
                if want_telolen == "?":
                    failures.append(
                        f"{fid} [{flags}] {scaffold}: teloLen is '?', an unstated expectation")
                    continue
                arm_rows = [b for b in terminal_bed.get(scaffold, [])
                            if b["teloType"] == "scaffold"]
                checks += 1
                if not arm_rows:
                    failures.append(
                        f"{fid} [{flags}] {scaffold}: expect_telolen needs a scaffold arm "
                        f"in the terminal BED, found none")
                    continue
                want_int = int(want_telolen)
                got_telolens = [b["teloLen"] for b in arm_rows]
                # A record can carry two arms; the intent names one teloLen, not a slot.
                got_telolen = want_int if want_int in got_telolens else got_telolens[0]
                if want_int not in got_telolens:
                    failures.append(
                        f"{fid} [{flags}] {scaffold}: teloLen expected {want_int!r}, "
                        f"got {got_telolen!r}\n    intent: {row['intent']}")
        finally:
            shutil.rmtree(outdir, ignore_errors=True)

    print(f"binary:  {TELOSCOPE}")
    print(f"fixtures: {len(by_run)} runs, {len(rows)} scaffold expectations, {checks} field checks")
    if failures:
        print(f"\n{len(failures)} disagreement(s):\n")
        for f in failures:
            print("  " + f)
        print("\nA disagreement means the fixture, the expectation, or the code is wrong; "
              "fix it deliberately, do not edit the manifest to match.")
        return 1
    print("PASS: every stated intent holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
