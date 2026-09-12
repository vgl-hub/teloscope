#!/usr/bin/env python3
"""Check every fixture against the intent declared for it in the synthetic manifest.
A disagreement is a finding (fixture, expectation, or code is wrong); record it in
docs/conflicts.md rather than editing the manifest to match.

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
WAIVERS = ROOT / "validateFiles" / "intent_waivers.tsv"
BLOCKED = ROOT / "validateFiles" / "intent_blocked.tsv"
REGISTER = ROOT / "docs" / "conflicts.md"

# Manifest column -> report column; expect_its is checked only in full-scan mode.
CHECKED = [
    ("expect_type", "type"),
    ("expect_anomaly", "anomaly"),
    ("expect_granular", "granular"),
    ("expect_telomeres", "telomeres"),
    ("expect_labels", "labels"),
    ("expect_gaps", "gaps"),
    ("expect_its", "its"),
]


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


def load_waivers(path):
    """(id, scaffold, column) -> register id for expectations the binary does not meet yet."""
    out = {}
    if not path.exists():
        return out
    for raw in path.read_text().splitlines():
        if raw.startswith("#") or not raw.strip():
            continue
        f = raw.split("\t")
        if len(f) < 5:
            raise SystemExit(f"{path}: malformed waiver row: {raw!r}")
        out[(f[0], f[1], f[2])] = f[3]
    return out


def open_register_ids(path):
    if not path.exists():
        return None
    import re
    return set(re.findall(r"REG-\d{3}", path.read_text()))


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

    waivers = load_waivers(WAIVERS)
    blocked_reg = load_waivers(BLOCKED)
    register = open_register_ids(REGISTER)
    if register is not None:
        unknown = {v for v in list(waivers.values()) + list(blocked_reg.values())
                   if v not in register}
        if unknown:
            raise SystemExit(
                f"{WAIVERS}: waivers cite register ids absent from {REGISTER}: "
                f"{sorted(unknown)}")

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

    failures, blocked, checks = [], [], 0
    waived, stale = [], []

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
                        rkey = (fid, scaffold, column)
                        reg = blocked_reg.get(rkey)
                        if reg is None:
                            failures.append(
                                f"{fid} [{flags}] {scaffold}: {column} is '?' but no row in "
                                f"{BLOCKED.name} says which conflict it waits on. An absent "
                                f"expectation defaults to '?', so this is as likely to be a "
                                f"forgotten field as a deliberate one.")
                        else:
                            blocked.append(f"{fid}/{scaffold}.{column} [{reg}]")
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
                    ok = got == want
                    reg = waivers.get((fid, scaffold, column))
                    if ok and reg:
                        stale.append(f"{fid}/{scaffold}.{column} (waived under {reg})")
                    elif not ok and reg:
                        waived.append(f"{fid}/{scaffold}.{column}: {column} should be "
                                      f"{want!r}, is {got!r} [{reg}]")
                    elif not ok:
                        failures.append(
                            f"{fid} [{flags}] {scaffold}: {column} expected {want!r}, "
                            f"got {got!r}\n    intent: {row['intent']}")
        finally:
            shutil.rmtree(outdir, ignore_errors=True)

    print(f"binary:  {TELOSCOPE}")
    print(f"fixtures: {len(by_run)} runs, {len(rows)} scaffold expectations, {checks} field checks")
    if blocked:
        print(f"blocked on docs/conflicts.md: {len(blocked)} fields ({', '.join(sorted(set(blocked))[:6])}...)")
    if waived:
        print(f"\nwaived, known open in docs/conflicts.md: {len(waived)}")
        for w in waived:
            print("  " + w)
    if stale:
        print(f"\n{len(stale)} STALE waiver(s): these now pass and must be removed from")
        print(f"{WAIVERS} together with the fix that made them pass:")
        for w in stale:
            print("  " + w)
        return 1
    if failures:
        print(f"\n{len(failures)} disagreement(s):\n")
        for f in failures:
            print("  " + f)
        print("\nA disagreement is a finding. Record it in docs/conflicts.md and fix the "
              "code or the expectation deliberately; do not edit the manifest to match.")
        return 1
    print("PASS: every stated intent holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
