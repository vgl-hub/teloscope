#!/usr/bin/env python3
"""Invariant checks for teloscope, with no recorded expected values anywhere.

Every check re-derives an output field via a different code path or compares two runs of
the tool against each other.

Usage:
    python3 scripts/check_invariants.py [--only GLOB] [--quick]
    TELOSCOPE=/path/to/other/teloscope python3 scripts/check_invariants.py
"""

import argparse
import bisect
import fnmatch
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))
TESTFILES = ROOT / "testFiles"
MANIFEST = TESTFILES / "synthetic" / "manifest.tsv"

# ---------------------------------------------------------------- readers (block BED, report)

BED_FIELDS = [
    "chrom", "start", "end", "teloLen", "teloLabel", "closestEnd", "fwdCan", "revCan",
    "fwdNonCan", "revNonCan", "chrSize", "teloType",
]
INT_FIELDS = {"start", "end", "teloLen", "fwdCan", "revCan", "fwdNonCan", "revNonCan", "chrSize"}


def read_block_bed(path):
    """Terminal or interstitial BED -> list of dicts, in file order."""
    rows = []
    for line in _significant_lines(path):
        f = line.split("\t")
        if len(f) != len(BED_FIELDS):
            raise ValueError(
                f"{path}: expected {len(BED_FIELDS)} columns, got {len(f)}: {line!r}")
        row = dict(zip(BED_FIELDS, f))
        for k in INT_FIELDS:
            row[k] = int(row[k])
        rows.append(row)
    return rows


def read_gaps_bed(path):
    out = []
    for line in _significant_lines(path):
        f = line.split("\t")
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def _significant_lines(path):
    if not path or not str(path) or not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip() or s.lstrip().startswith(("#", "track", "browser")):
                continue
            out.append(s)
    return out


def read_report(text_or_path, is_text=False):
    """Parse the per-path table and the summary; columns are read by header name, not position."""
    if is_text:
        text = text_or_path
    else:
        with open(text_or_path) as fh:
            text = fh.read()

    columns, rows, summary, params = None, {}, {}, {}
    for line in text.splitlines():
        if line.startswith("#params"):
            for m in re.finditer(r"(\w+)=(\S+)", line):
                params[m.group(1)] = m.group(2)
            continue
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[:2] == ["pos", "header"]:
            columns = fields
            continue
        if columns and fields and fields[0].isdigit():
            row = dict(zip(columns, fields))
            rows[row["header"]] = row
            continue
        if len(fields) == 2 and fields[0].endswith(":"):
            summary[fields[0][:-1]] = fields[1]
    return rows, summary, params


# ---------------------------------------------------------------- small re-derivations

# closestEnd is the end a telomere belongs to, not the nearer one; arms, anomalies and granular key on it

LABEL_SCALE = 1000
DEFAULT_LABEL_THRESHOLD = 0.667


def strand_label(forward_count, total_count, threshold=DEFAULT_LABEL_THRESHOLD):
    """Strand label from counts, symmetric thirds as in computeStrandLabel(threshold)."""
    if total_count == 0:
        return "b"
    scaled = forward_count * LABEL_SCALE
    hi = round(threshold * LABEL_SCALE)
    if scaled > total_count * hi:
        return "p"
    if scaled < total_count * (LABEL_SCALE - hi):
        return "q"
    return "b"


def arms(blocks):
    """The scaffold-terminal rows: at most one with closestEnd p, one with closestEnd q."""
    scaffold = [b for b in blocks if b["teloType"] == "scaffold"]
    p = next((b for b in scaffold if b["closestEnd"] == "p"), None)
    q = next((b for b in scaffold if b["closestEnd"] == "q"), None)
    return p, q


def scaffold_type(blocks):
    """Both arms present t2t, one incomplete, neither none."""
    p, q = arms(blocks)
    n = (p is not None) + (q is not None)
    return {2: "t2t", 1: "incomplete", 0: "none"}[n]


ANOMALY_ORDER = ["discordant_p", "discordant_q", "fragmented_p", "fragmented_q"]


def anomalies(blocks):
    """discordant: arm strand != its end; fragmented: an arm has more than one piece."""
    p, q = arms(blocks)
    flags = set()
    for b, suffix in ((p, "p"), (q, "q")):
        if b is None:
            continue
        if b["teloLabel"] != b["closestEnd"]:
            flags.add(f"discordant_{suffix}")
        if b["teloLen"] < b["end"] - b["start"]:
            flags.add(f"fragmented_{suffix}")
    ordered = [f for f in ANOMALY_ORDER if f in flags]
    return ",".join(ordered) if ordered else "."


def granular(blocks):
    """One token per terminal row, ascending by start: upper scaffold, lower contig, * discordant."""
    blocks = sorted(blocks, key=lambda b: b["start"])
    out = []
    for b in blocks:
        ch = b["closestEnd"].upper() if b["teloType"] == "scaffold" else b["closestEnd"].lower()
        if b["teloLabel"] != b["closestEnd"]:
            ch += "*"
        out.append(ch)
    return "".join(out)


def labels_column(blocks):
    """closestEnd (side) of each scaffold arm, ascending by start -- the uppercase letters
    of granular, lowercased, in order -- or 'none'."""
    scaffold = sorted((b for b in blocks if b["teloType"] == "scaffold"), key=lambda b: b["start"])
    out = "".join(b["closestEnd"] for b in scaffold)
    return out or "none"


class GapIndex:
    """Gap runs of one record, indexed once for O(log n) overlap queries."""

    def __init__(self, gaps):
        self.gaps = list(gaps)
        self.starts = [s for s, _ in self.gaps]

    def has_gap_between(self, a, b):
        """True when any gap run overlaps [a, b)."""
        i = bisect.bisect_left(self.starts, b)
        return i > 0 and self.gaps[i - 1][1] > a


class Recorder:
    def __init__(self):
        self.rows = []

    def check(self, invariant, subject, ok, detail=""):
        self.rows.append({"invariant": invariant, "subject": subject,
                          "pass": bool(ok), "detail": detail})

    def eq(self, invariant, subject, got, want, note=""):
        self.check(invariant, subject, got == want,
                   f"expected {want!r}, got {got!r}{(' ' + note) if note else ''}")


def run(args, outdir):
    argv = [str(TELOSCOPE), *[str(a) for a in args], "-o", str(outdir)]
    return subprocess.run(argv, capture_output=True, timeout=600, check=False)


def outputs(outdir, stem):
    d = pathlib.Path(outdir)
    return {
        "terminal": d / f"{stem}_terminal_telomeres.bed",
        "interstitial": d / f"{stem}_interstitial_telomeres.bed",
        "gaps": d / f"{stem}_gaps.bed",
        "report": d / f"{stem}_report.tsv",
    }


def block_spans(path):
    """(chrom, start, end) of every row in a block BED, read fresh from disk."""
    return {(b["chrom"], b["start"], b["end"]) for b in read_block_bed(path)}


# ---------------------------------------------------------------- families

def contig_bounds(gaps_sorted, pos, seq_len):
    """[start, end) of the contig containing pos, from a sorted list of gap runs."""
    cs = max((ge for gs, ge in gaps_sorted if ge <= pos), default=0)
    ce = min((gs for gs, ge in gaps_sorted if gs >= pos), default=seq_len)
    return cs, ce


def check_one_arm_per_end(rec, subject, terminal, gaps_by_chrom, manual_curation):
    """At most one scaffold arm per end; with -n, at most one terminal row per end per contig."""
    by_chrom = {}
    for b in terminal:
        by_chrom.setdefault(b["chrom"], []).append(b)
    for chrom, blocks in by_chrom.items():
        scaffold = [b for b in blocks if b["teloType"] == "scaffold"]
        for end in ("p", "q"):
            n = sum(1 for b in scaffold if b["closestEnd"] == end)
            rec.check("DER-19-one-arm-per-end", f"{subject}:{chrom}:{end}", n <= 1,
                      f"{n} scaffold rows with closestEnd {end}")
        if not manual_curation:
            continue
        seq_len = blocks[0]["chrSize"]
        gaps_sorted = sorted(gaps_by_chrom.get(chrom, []))
        by_contig = {}
        for b in blocks:
            bounds = contig_bounds(gaps_sorted, b["start"], seq_len)
            by_contig.setdefault(bounds, []).append(b)
        for bounds, cblocks in by_contig.items():
            for end in ("p", "q"):
                n = sum(1 for b in cblocks if b["closestEnd"] == end)
                rec.check("DER-19-one-arm-per-end", f"{subject}:{chrom}:{bounds}:{end}", n <= 1,
                          f"{n} terminal rows with closestEnd {end} on contig {bounds}")


def check_junction_class(rec, subject, terminal, interstitial, gaps_by_chrom, max_block_dist):
    """Recompute each interstitial row's junction class from its nearest same-contig neighbour."""
    by_chrom = {}
    for b in terminal + interstitial:
        by_chrom.setdefault(b["chrom"], []).append(b)
    its_ids = {id(b) for b in interstitial}
    for chrom, blocks in by_chrom.items():
        ordered = sorted(blocks, key=lambda b: b["start"])
        gap_index = GapIndex(sorted(gaps_by_chrom.get(chrom, [])))
        for i, b in enumerate(ordered):
            if id(b) not in its_ids:
                continue
            candidates = []
            if i > 0 and not gap_index.has_gap_between(ordered[i - 1]["end"], b["start"]):
                candidates.append((b["start"] - ordered[i - 1]["end"], ordered[i - 1], b))
            if i + 1 < len(ordered) and not gap_index.has_gap_between(b["end"], ordered[i + 1]["start"]):
                candidates.append((ordered[i + 1]["start"] - b["end"], b, ordered[i + 1]))
            candidates = [c for c in candidates if c[0] <= max_block_dist]
            if not candidates:
                expect = "single"
            else:
                _, left_b, right_b = min(candidates, key=lambda c: c[0])
                l, r = left_b["teloLabel"], right_b["teloLabel"]
                expect = ("single" if "b" in (l, r) else
                          "fusion" if (l, r) == ("q", "p") else
                          "tail_to_tail" if (l, r) == ("p", "q") else
                          "fragmentation" if l == r else "single")
            rec.eq("DER-20-junction-class", f"{subject}:{chrom}:{b['start']}", b["teloType"], expect)


def check_derivations(rec, subject, terminal, interstitial, gaps, rows, summary, params):
    by_chrom = {}
    for b in terminal:
        by_chrom.setdefault(b["chrom"], []).append(b)
    its_by_chrom = {}
    for b in interstitial:
        its_by_chrom.setdefault(b["chrom"], []).append(b)
    gaps_by_chrom = {}
    for c, s, e in gaps:
        gaps_by_chrom.setdefault(c, []).append((s, e))

    threshold = float(params.get("label_threshold", 0.667))
    for b in terminal + interstitial:
        where = f"{subject}:{b['chrom']}:{b['start']}"
        forward = b["fwdCan"] + b["fwdNonCan"]
        total = forward + b["revCan"] + b["revNonCan"]
        rec.eq("DER-01-strand-label-from-counts", where, b["teloLabel"],
               strand_label(forward, total, threshold),
               f"(fwd {forward} of {total})")
    for b in terminal:
        where = f"{subject}:{b['chrom']}:{b['start']}"
        rec.check("DER-01-terminal-label-is-strand-pure", where, b["teloLabel"] != "b",
                  "a terminal row's teloLabel must be p or q, never b (R1)")

    for chrom, row in rows.items():
        blocks = sorted(by_chrom.get(chrom, []), key=lambda b: b["start"])
        where = f"{subject}:{chrom}"
        rec.eq("DER-03-type-from-blocks", where, row.get("type"), scaffold_type(blocks))
        if "anomaly" in row:
            rec.eq("DER-04-anomaly-from-blocks", where, row["anomaly"], anomalies(blocks))
        rec.eq("DER-05-granular-from-blocks", where, row.get("granular", ""), granular(blocks))
        rec.eq("DER-06-labels-from-blocks", where, row.get("labels"), labels_column(blocks))
        # DER-03 tests the same predicate; this checks the granular string, built separately.
        rec.eq("DER-07-telomere-count-matches-granular", where,
               int(row.get("telomeres", -1)),
               sum(1 for c in row.get("granular", "") if c in "PQ"))
        rec.eq("DER-08-gap-count", where, int(row.get("gaps", -1)),
               len(gaps_by_chrom.get(chrom, [])))
        if "its" in row:
            rec.eq("DER-09-its-count", where, int(row["its"]),
                   len(its_by_chrom.get(chrom, [])))

    manual_curation = params.get("manual_curation") == "true"
    check_one_arm_per_end(rec, subject, terminal, gaps_by_chrom, manual_curation)
    max_block_dist = int(params.get("max_block_dist", 1000))
    check_junction_class(rec, subject, terminal, interstitial, gaps_by_chrom, max_block_dist)

    if summary:
        total_paths = int(summary.get("Total paths", -1))
        buckets = ["T2T", "Gapped T2T", "Incomplete", "Gapped incomplete",
                   "No telomeres", "Gapped no telomeres"]
        if all(b in summary for b in buckets):
            rec.eq("DER-10-completeness-buckets-partition", subject,
                   sum(int(summary[b]) for b in buckets), total_paths)
        counts = ["Two telomeres", "One telomere", "Zero telomeres"]
        if all(c in summary for c in counts):
            rec.eq("DER-11-telomere-counts-partition", subject,
                   sum(int(summary[c]) for c in counts), total_paths)
        if "Total gaps" in summary:
            rec.eq("DER-12-total-gaps", subject, int(summary["Total gaps"]), len(gaps))


FULL_SCAN_FLAGS = {"-r", "-g", "-e", "-m", "-i", "--out-win-repeats", "--out-gc",
                   "--out-entropy", "--out-matches", "--out-its"}


def shared_report(rows):
    keys = ["telomeres", "labels", "gaps", "type", "anomaly", "granular"]
    return {h: tuple(r.get(k) for k in keys) for h, r in rows.items()}


def check_metamorphic(rec, subject, fasta_path, base_flags, stem):
    """Relations between runs. No oracle needed, so these hold on any input."""
    def once_with(flags, extra):
        d = tempfile.mkdtemp(prefix="telo_inv_")
        p = run(["-f", str(fasta_path), *flags, *extra], d)
        if p.returncode != 0:
            return d, None
        rows, summary, _ = read_report(p.stdout.decode("utf-8", "replace"), is_text=True)
        return d, (rows, summary)

    def once(extra):
        return once_with(base_flags, extra)

    def match_set(d, name):
        p = pathlib.Path(d) / name
        if not p.exists():
            return None
        return {tuple(line.split("\t")[:3]) for line in p.read_text().splitlines()
                if line and not line.startswith(("#", "track"))}

    dirs = []
    try:
        d1, a = once([])
        dirs.append(d1)
        if a is None:
            rec.check("MET-00-run-succeeds", subject, False, "teloscope exited non-zero")
            return

        d2, b = once([])
        dirs.append(d2)
        rec.check("MET-01-deterministic", subject, b is not None and shared_report(a[0]) == shared_report(b[0]),
                  "two identical runs disagree")

        d3, c = once(["-j", "1"])
        dirs.append(d3)
        d4, e = once(["-j", "8"])
        dirs.append(d4)
        if c and e:
            rec.check("MET-02-thread-invariant", subject, shared_report(c[0]) == shared_report(e[0]),
                      "-j 1 and -j 8 disagree")

        # -n runs are compared with -n -i in MET-08, not with the generic scan-mode check
        manual_curation = any(f in ("-n", "--manual-curation") for f in base_flags)
        if manual_curation:
            ultra_run = full_run = None
        elif any(x in FULL_SCAN_FLAGS for x in base_flags):
            ultra_flags = [x for x in base_flags if x not in FULL_SCAN_FLAGS]
            d5, u = once_with(ultra_flags, [])
            dirs.append(d5)
            ultra_run, full_run, ultra_dir, full_dir = u, a, d5, d1
        else:
            d5, f = once(["-i"])
            dirs.append(d5)
            ultra_run, full_run, ultra_dir, full_dir = a, f, d1, d5
        if ultra_run and full_run:
            keys = ["telomeres", "labels", "gaps", "type", "anomaly", "granular"]
            ultra = {h: tuple(r.get(k) for k in keys) for h, r in ultra_run[0].items()}
            full = {h: tuple(r.get(k) for k in keys) for h, r in full_run[0].items()}
            bad = [h for h in ultra if ultra.get(h) != full.get(h)]
            rec.check("MET-03-scan-mode-invariant", subject, not bad,
                      f"differs on {bad[:4]}: ultra {[ultra[h] for h in bad[:2]]} "
                      f"vs full {[full.get(h) for h in bad[:2]]}")
            ub = block_spans(outputs(ultra_dir, stem)["terminal"])
            fb = block_spans(outputs(full_dir, stem)["terminal"])
            rec.check("MET-04-scan-mode-blocks", subject, ub == fb,
                      f"terminal blocks differ: only ultra {sorted(ub - fb)[:3]}, "
                      f"only full {sorted(fb - ub)[:3]}")

        # -x is monotonic at the match level, not at the block level
        if "-x" not in base_flags and "-m" not in base_flags:
            sets = {}
            for x in ("0", "1", "2"):
                dx, rx = once(["-m", "-x", x])
                dirs.append(dx)
                if rx is None:
                    break
                sets[x] = (match_set(dx, f"{stem}_canonical_matches.bed"),
                           match_set(dx, f"{stem}_noncanonical_matches.bed"))
            if len(sets) == 3:
                if any(v is None for pair in sets.values() for v in pair):
                    rec.check("MET-05-edit-distance-matches", subject, False,
                              "-m wrote no match BED")
                else:
                    rec.check("MET-05-edit-distance-matches", subject,
                              sets["0"][0] == sets["1"][0] == sets["2"][0]
                              and sets["0"][1] <= sets["1"][1] <= sets["2"][1],
                              f"canonical {[len(sets[x][0]) for x in '012']}, "
                              f"non-canonical {[len(sets[x][1]) for x in '012']}")

        if not any(f in ("-t", "--terminal-limit") for f in base_flags):
            d6, tip = once(["-t", "5000"])
            dirs.append(d6)
            if tip:
                default_terminal = block_spans(outputs(d1, stem)["terminal"])
                tip_terminal = block_spans(outputs(d6, stem)["terminal"])
                rec.check("MET-06-tip-window-invariant", subject,
                          default_terminal == tip_terminal,
                          f"-t 5000 terminal BED differs from the default: only default "
                          f"{sorted(default_terminal - tip_terminal)[:3]}, only -t 5000 "
                          f"{sorted(tip_terminal - default_terminal)[:3]}")

        no_n_flags = [f for f in base_flags if f not in ("-n", "--manual-curation")]
        no_n_full_flags = no_n_flags if any(f in FULL_SCAN_FLAGS for f in no_n_flags) \
            else no_n_flags + ["-i"]
        n_flags = no_n_flags + ["-n"]
        d7, no_n_res = once_with(no_n_full_flags, [])
        dirs.append(d7)
        d8, n_res = once_with(n_flags, [])
        dirs.append(d8)
        if no_n_res and n_res:
            n_terminal = read_block_bed(outputs(d8, stem)["terminal"])
            n_set = {(b["chrom"], b["start"], b["end"]) for b in n_terminal}
            no_n_set = block_spans(outputs(d7, stem)["terminal"])
            contig_rows = {(b["chrom"], b["start"], b["end"])
                           for b in n_terminal if b["teloType"] == "contig"}
            rec.check("MET-07-manual-curation", subject,
                      n_set == no_n_set | contig_rows,
                      f"-n terminal BED should equal the full-scan terminal BED (without -n) "
                      f"plus contig rows: missing {sorted((no_n_set | contig_rows) - n_set)[:3]}, "
                      f"extra {sorted(n_set - (no_n_set | contig_rows))[:3]}")
            its_set = block_spans(outputs(d8, stem)["interstitial"])
            overlap = contig_rows & its_set
            rec.check("MET-07-manual-curation", subject + ":its", not overlap,
                      f"contig rows {sorted(overlap)[:3]} also appear in the -n run's "
                      f"own interstitial BED")

        if not any(f in FULL_SCAN_FLAGS for f in base_flags):
            d9, nf = once(["-n"])
            dirs.append(d9)
            d10, nfull = once(["-n", "-i"])
            dirs.append(d10)
            if nf and nfull:
                fast_terminal = block_spans(outputs(d9, stem)["terminal"])
                full_terminal = block_spans(outputs(d10, stem)["terminal"])
                rec.check("MET-08-manual-curation-fast-eq-full", subject,
                          fast_terminal == full_terminal,
                          f"-n and -n -i terminal BEDs differ: only -n "
                          f"{sorted(fast_terminal - full_terminal)[:3]}, only -n -i "
                          f"{sorted(full_terminal - fast_terminal)[:3]}")
    finally:
        for d in dirs:
            shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------- driver

def manifest_runs():
    """Every (fixture, flags) the manifest declares, deduplicated."""
    if not MANIFEST.exists():
        return []
    seen, out, header = set(), [], None
    for raw in MANIFEST.read_text().splitlines():
        if raw.startswith("#") or not raw.strip():
            continue
        f = raw.split("\t")
        if header is None:
            header = f
            continue
        row = dict(zip(header, f))
        key = (row["path"], row["flags"])
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="*", help="glob over invariant ids")
    ap.add_argument("--quick", action="store_true", help="skip the real genomes")
    args = ap.parse_args()

    if not TELOSCOPE.exists():
        raise SystemExit(f"teloscope binary not found at {TELOSCOPE}; run `make head`")

    rec = Recorder()
    subjects = list(manifest_runs())
    if not subjects:
        raise SystemExit(
            f"no subjects: {MANIFEST} is missing or empty. Run testFiles/"
            "generate_synthetic.sh. A harness with nothing to check must not report PASS.")
    if not args.quick:
        for g in ["bTaeGut7_chr33_mat.fa.gz", "bTaeGut7_chr33_pat.fa.gz",
                   "vgp_probe.fa.gz", "vgp_turtle.fa.gz"]:
            if (TESTFILES / g).exists():
                subjects.append((g, "-i"))

    ran, missing = 0, []
    for path, flags in subjects:
        fasta_path = TESTFILES / path
        if not fasta_path.exists():
            missing.append(path)
            continue
        base_flags = shlex.split(flags) if flags and flags != "-" else []
        stem = fasta_path.name
        subject = f"{stem} [{flags}]"

        ran += 1
        outdir = tempfile.mkdtemp(prefix="telo_inv_")
        try:
            proc = run(["-f", str(fasta_path), *base_flags], outdir)
            if proc.returncode != 0:
                rec.check("RUN-00-exit-zero", subject, False,
                          proc.stderr.decode("utf-8", "replace").strip()[:200])
                continue
            o = outputs(outdir, stem)
            terminal = read_block_bed(o["terminal"])
            interstitial = read_block_bed(o["interstitial"])
            gaps = read_gaps_bed(o["gaps"])
            rows, summary, params = read_report(
                proc.stdout.decode("utf-8", "replace"), is_text=True)
            if o["report"].exists():
                _, _, params = read_report(o["report"])

            check_derivations(rec, subject, terminal, interstitial, gaps, rows, summary, params)
        finally:
            shutil.rmtree(outdir, ignore_errors=True)

        check_metamorphic(rec, subject, fasta_path, base_flags, stem)

    if missing:
        raise SystemExit(
            f"{len(missing)} declared fixture(s) are missing from {TESTFILES}: "
            f"{missing[:5]}. Run testFiles/generate_synthetic.sh.")
    if ran != len(subjects):
        raise SystemExit(f"only ran {ran} of {len(subjects)} declared subjects")

    rows = [r for r in rec.rows if fnmatch.fnmatch(r["invariant"], args.only)]
    if not rows:
        raise SystemExit(f"no checks matched --only {args.only!r}; refusing to report PASS")

    families = {}
    for r in rows:
        fam = r["invariant"].split("-")[0]
        d = families.setdefault(fam, [0, 0])
        d[0] += 1
        if not r["pass"]:
            d[1] += 1

    failures = [r for r in rows if not r["pass"]]

    print(f"binary:   {TELOSCOPE}")
    print(f"subjects: {len(subjects)}")
    for fam in sorted(families):
        total, bad = families[fam]
        print(f"  {fam}-*  {total - bad}/{total} pass")
    if failures:
        seen = {}
        for f in failures:
            seen.setdefault(f["invariant"], []).append(f)
        print(f"\n{len(failures)} failure(s) across {len(seen)} invariant(s):\n")
        for inv in sorted(seen):
            group = seen[inv]
            print(f"  {inv}  ({len(group)} failing)")
            for f in group[:3]:
                print(f"     {f['subject']}  {f['detail']}")
            if len(group) > 3:
                print(f"     ... and {len(group) - 3} more")
        return 1
    print("PASS: every invariant holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
