#!/usr/bin/env python3
"""Invariant checks for teloscope, with no recorded expected values anywhere.

Every check re-derives an output field via a different code path, compares two runs of the
tool against each other, or measures the output against the input FASTA directly.

Usage:
    python3 scripts/check_invariants.py [--only GLOB] [--quick]
    TELOSCOPE=/path/to/other/teloscope python3 scripts/check_invariants.py
"""

import argparse
import bisect
import fnmatch
import os
import pathlib
import shlex
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import teloscope_model as M  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_TELOSCOPE = ROOT / "build/bin" / ("teloscope.exe" if os.name == "nt" else "teloscope")
TELOSCOPE = pathlib.Path(os.environ.get("TELOSCOPE", DEFAULT_TELOSCOPE))
TESTFILES = ROOT / "testFiles"
MANIFEST = TESTFILES / "synthetic" / "manifest.tsv"
WAIVERS = ROOT / "validateFiles" / "invariant_waivers.tsv"
REGISTER = ROOT / "docs" / "conflicts.md"


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


# ---------------------------------------------------------------- families

def check_derivations(rec, subject, terminal, interstitial, gaps, rows, summary, fasta):
    lengths = {h: len(s) for h, s in fasta}
    by_chrom = {}
    for b in terminal:
        by_chrom.setdefault(b["chrom"], []).append(b)
    its_by_chrom = {}
    for b in interstitial:
        its_by_chrom.setdefault(b["chrom"], []).append(b)
    gaps_by_chrom = {}
    for c, s, e in gaps:
        gaps_by_chrom.setdefault(c, []).append((s, e))

    for b in terminal + interstitial:
        where = f"{subject}:{b['chrom']}:{b['start']}"
        total = b["forwardCount"] + b["reverseCount"]
        rec.eq("DER-01-strand-label-from-counts", where, b["label"],
               M.strand_label(b["forwardCount"], total),
               f"(fwd {b['forwardCount']} of {total})")
        overlaps = any(not (e <= b["start"] or s >= b["end"])
                       for s, e in gaps_by_chrom.get(b["chrom"], []))
        rec.eq("DER-02-gapstatus-from-gaps-bed", where, b["gapStatus"],
               "gapped" if overlaps else "contiguous")

    for chrom, row in rows.items():
        blocks = sorted(by_chrom.get(chrom, []), key=lambda b: b["start"])
        where = f"{subject}:{chrom}"
        rec.eq("DER-03-type-from-blocks", where, row.get("type"), M.scaffold_type(blocks))
        if "anomaly" in row:
            rec.eq("DER-04-anomaly-from-blocks", where, row["anomaly"], M.anomalies(blocks))
        rec.eq("DER-05-granular-from-blocks", where, row.get("granular", ""), M.granular(blocks))
        rec.eq("DER-06-labels-from-blocks", where, row.get("labels"), M.labels_column(blocks))
        # DER-03 tests the same predicate; this checks the granular string, built separately.
        rec.eq("DER-07-telomere-count-matches-granular", where,
               int(row.get("telomeres", -1)),
               sum(1 for c in row.get("granular", "") if c in "PQ"))
        rec.eq("DER-08-gap-count", where, int(row.get("gaps", -1)),
               len(gaps_by_chrom.get(chrom, [])))
        if "its" in row:
            rec.eq("DER-09-its-count", where, int(row["its"]),
                   len(its_by_chrom.get(chrom, [])))

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
        if "Scaffold N50" in summary:
            rec.eq("DER-13-scaffold-n50", subject, int(summary["Scaffold N50"]),
                   M.n50(list(lengths.values())))
        if "Contig N50" in summary:
            contigs = []
            for _, seq in fasta:
                contigs += M.contig_lengths(len(seq), M.gap_runs(seq))
            rec.eq("DER-14-contig-n50", subject, int(summary["Contig N50"]), M.n50(contigs))


FLAG_PARAMS = {"-t": "terminal_limit", "--terminal-limit": "terminal_limit",
               "-k": "max_match_dist", "--max-match-distance": "max_match_dist",
               "-d": "max_block_dist", "--max-block-distance": "max_block_dist",
               "-l": "min_block_len", "--min-block-length": "min_block_len",
               "-y": "min_block_density", "--min-block-density": "min_block_density",
               "-x": "edit_distance", "--edit-distance": "edit_distance",
               "-c": "canonical", "--canonical": "canonical",
               "--terminal-tolerance": "terminal_tolerance",
               "--min-its-length": "min_its_length",
               "--min-block-counts": "min_block_counts"}
FULL_SCAN_FLAGS = {"-r", "-g", "-e", "-m", "-i", "--out-win-repeats", "--out-gc",
                   "--out-entropy", "--out-matches", "--out-its"}


def params_from_flags(header_params, flags):
    """Header params overridden by the flags actually passed; older binaries write no header."""
    out = dict(header_params)
    it = iter(range(len(flags)))
    for i in it:
        f = flags[i]
        if f in FLAG_PARAMS and i + 1 < len(flags):
            out[FLAG_PARAMS[f]] = flags[i + 1]
            next(it, None)
        elif "=" in f and f.split("=")[0] in FLAG_PARAMS:
            out[FLAG_PARAMS[f.split("=")[0]]] = f.split("=", 1)[1]
        elif len(f) > 2 and f[:2] in FLAG_PARAMS and not f.startswith("--"):
            out[FLAG_PARAMS[f[:2]]] = f[2:]
    if any(f in FULL_SCAN_FLAGS for f in flags):
        out["ultra_fast"] = "false"
    elif "ultra_fast" not in out:
        out["ultra_fast"] = "true"
    return out


def check_oracle(rec, subject, terminal, interstitial, gaps, fasta, params):
    canonical = params.get("canonical", "CCCTAA/TTAGGG").split("/")[-1]
    edit = int(params.get("edit_distance", 1))
    min_len = int(params.get("min_block_len", 300))
    min_its = int(params.get("min_its_length", 100))
    # headLimit = min(tolerance, firstBase + terminalLimit) (src/teloscope.cpp:345-346)
    tolerance = min(int(params.get("terminal_tolerance", 2000)),
                    int(params.get("terminal_limit", 50000)))

    seqs = dict(fasta)

    # The gaps BED is a pure function of the input and nothing else reads it back.
    observed = {}
    for c, s, e in gaps:
        observed.setdefault(c, []).append((s, e))
    for header, seq in fasta:
        rec.eq("ORA-01-gaps-are-the-N-runs", f"{subject}:{header}",
               sorted(observed.get(header, [])), M.gap_runs(seq))

    density = float(params.get("min_block_density", 0.5))
    gaps_by_chrom = {}
    for c, gs, ge in gaps:
        gaps_by_chrom.setdefault(c, []).append((gs, ge))

    # terminal blocks: canonical coverage (src/teloscope.cpp:403); interstitial: all matches (:461)
    for b, canonical_only in [(x, True) for x in terminal] + [(x, False) for x in interstitial]:
        seq = seqs.get(b["chrom"])
        if seq is None:
            continue
        where = f"{subject}:{b['chrom']}:{b['start']}"
        sub = seq[b["start"]:b["end"]]
        canon_iv, all_iv = M.locate_matches(sub, canonical, edit)
        merged = M.merge_intervals(canon_iv if canonical_only else all_iv)
        covered = M.covered_bases(merged, 0, len(sub))
        gap_bases = sum(max(0, min(ge, b["end"]) - max(gs, b["start"]))
                        for gs, ge in gaps_by_chrom.get(b["chrom"], []))
        called = b["blockLen"] - gap_bases
        # the gate is union coverage over called bases (src/teloscope.cpp:403, :461), not a count
        rec.check("ORA-02-block-meets-canonical-density", where,
                  called <= 0 or covered >= density * called,
                  f"independently measured {covered} "
                  f"{'canonical' if canonical_only else 'matched'} bases over {called} "
                  f"called bases = {covered / called if called else 0:.3f}, -y is {density}")

    for b in terminal:
        seq = seqs.get(b["chrom"])
        if seq is None:
            continue
        where = f"{subject}:{b['chrom']}:{b['start']}"
        rec.check("DER-17-terminal-length-floor", where, b["blockLen"] >= min_len,
                  f"blockLen {b['blockLen']} against -l {min_len}")
        # tolerance counts called bases, so a leading N run hides nothing (docs/parameters.md:91)
        from_start = M.called_offset_from_start(seq, b["start"])
        from_end = M.called_offset_from_end(seq, b["end"])
        rec.check("ORA-03-terminal-block-is-terminal", where,
                  min(from_start, from_end) <= tolerance,
                  f"{from_start} called bases from start, {from_end} from end, "
                  f"--terminal-tolerance is {tolerance}")

    check_recall(rec, subject, terminal, interstitial, fasta, params)

    for b in interstitial:
        where = f"{subject}:{b['chrom']}:{b['start']}"
        rec.check("DER-18-its-length-floor", where, b["blockLen"] >= min_its,
                  f"blockLen {b['blockLen']} against --min-its-length {min_its}")

    # Coordinates must be disjoint between terminal and interstitial calls (DER-15/16).
    for chrom in {b["chrom"] for b in terminal}:
        t = [(b["start"], b["end"]) for b in terminal if b["chrom"] == chrom]
        i = [(b["start"], b["end"]) for b in interstitial if b["chrom"] == chrom]
        for ts, te in t:
            for is_, ie in i:
                rec.check("DER-15-terminal-and-its-disjoint", f"{subject}:{chrom}",
                          te <= is_ or ie <= ts, f"terminal {ts}-{te} overlaps ITS {is_}-{ie}")
        for a in range(len(t)):
            for b2 in range(a + 1, len(t)):
                rec.check("DER-16-terminal-blocks-disjoint", f"{subject}:{chrom}",
                          t[a][1] <= t[b2][0] or t[b2][1] <= t[a][0],
                          f"{t[a]} overlaps {t[b2]}")


def check_recall(rec, subject, terminal, interstitial, fasta, params):
    """Verify reported blocks cover dense canonical FASTA windows within -d + -k + motif slack."""
    canonical = params.get("canonical", "CCCTAA/TTAGGG").split("/")[-1]
    motif = len(canonical)
    min_counts = int(params.get("min_block_counts", 2))
    min_len = int(params.get("min_block_len", 300))
    min_its = int(params.get("min_its_length", 100))
    density = float(params.get("min_block_density", 0.5))
    max_dist = int(params.get("max_block_dist", 500))
    match_dist = int(params.get("max_match_dist", 50))
    terminal_limit = int(params.get("terminal_limit", 50000))
    tolerance = min(int(params.get("terminal_tolerance", 2000)), terminal_limit)
    full_scan = params.get("ultra_fast", "true") == "false"
    slack = max_dist + match_dist + motif

    t_by = {}
    for b in terminal:
        t_by.setdefault(b["chrom"], []).append(b)
    i_by = {}
    for b in interstitial:
        i_by.setdefault(b["chrom"], []).append(b)

    for header, seq in fasta:
        runs = M.gap_runs(seq)
        gaps = M.GapIndex(runs)
        n = len(seq)
        if gaps.called_before(n) == 0:
            continue
        first_base = runs[0][1] if runs and runs[0][0] == 0 else 0
        last_base = runs[-1][0] if runs and runs[-1][1] == n else n
        head_called, tail_called = gaps.called_before(first_base), gaps.called_before(last_base)
        t_blocks = t_by.get(header, [])
        i_blocks = i_by.get(header, [])
        t_iv = [(b["start"], b["end"]) for b in t_blocks]
        merged_t = M.merge_intervals(t_iv)
        all_iv = M.merge_intervals(t_iv + ([(b["start"], b["end"]) for b in i_blocks]
                                           if full_scan else []))

        def hit_of(iv):
            for b in t_blocks:
                if b["start"] < iv[1] and iv[0] < b["end"]:
                    return b
            return None

        matches = M.canonical_positions(seq, canonical)
        windows = M.dense_windows(matches, gaps, min_len, density, max_dist, min_counts)

        # Clipped to the -t extent the tool may report (src/teloscope.cpp:391, :428).
        starts = [s for s, _, _ in matches]

        def clip(w, lo, hi):
            a, b = bisect.bisect_left(starts, max(w[0], lo)), bisect.bisect_left(starts, min(w[1], hi))
            inside = [(s, e, o) for s, e, o in matches[a:b] if e <= min(w[1], hi) and o == w[2]]
            if len(inside) < min_counts:
                return None
            s, e = inside[0][0], inside[-1][1]
            if e - s < min_len or len(inside) * motif < density * (e - s):
                return None
            return (s, e, w[2], len(inside))

        ends, both_ends = {}, set()
        for w in windows:
            near_p = gaps.called_before(w[0]) - head_called <= tolerance
            near_q = tail_called - gaps.called_before(w[1]) <= tolerance
            if near_p and near_q:
                both_ends.add(w)
            if near_p:
                c = clip(w, 0, min(n, first_base + terminal_limit))
                if c:
                    ends.setdefault("p", []).append((c, w))
            if near_q:
                c = clip(w, max(0, last_base - terminal_limit), n)
                if c:
                    ends.setdefault("q", []).append((c, w))

        for end, cands in ends.items():
            cands.sort(key=lambda cw: cw[0][0] if end == "p" else -cw[0][1])
            outer, raw = cands[0]
            where = f"{subject}:{header}:{end}:{outer[0]}-{outer[1]}"
            hit = hit_of(outer)
            rec.check("ORA-06-terminal-recall", where, hit is not None,
                      f"{outer[3]} canonical {outer[2]} matches over {outer[1] - outer[0]} bp "
                      f"start within {tolerance} called bases of the {end} end, no terminal "
                      f"block reported over it")
            if hit is not None:
                label = hit.get("label")
                rec.check("ORA-06-terminal-recall", where + ":label", label in (outer[2], "b"),
                          f"window is {outer[2]}-oriented, overlapping block is labelled {label}")
                arm = hit.get("arm")
                if raw not in both_ends:
                    rec.check("ORA-06-terminal-recall", where + ":arm", arm == end,
                              f"window sits at the {end} end, overlapping block says arm {arm}")
            targets = cands if full_scan else cands[:1]
            for i, (c, _) in enumerate(targets):
                cov = M.covered_bases(merged_t if i == 0 else all_iv, c[0], c[1])
                length = c[1] - c[0]
                # an overlapping block is >= -l long and starts <= -k + motif outside the window
                need = max(length - slack, min(length, min_len) - match_dist - motif)
                rec.check("ORA-07-terminal-window-covered",
                          f"{subject}:{header}:{end}:{c[0]}-{c[1]}", cov >= need,
                          f"reported blocks cover {cov} of {length} bp, need {need} "
                          f"(slack -d+-k+motif = {slack})")

        if not full_scan:
            continue

        # Reported blocks must cover pure tandem clusters >= --min-its-length (REG-003).
        clusters = M.merge_runs(M.pure_tandem_runs(seq, canonical, 4), gaps, max_dist, density)
        arms = {b["arm"] for b in t_blocks}
        has_p, has_q = "p" in arms, "q" in arms
        for s, e, o in clusters:
            if e - s < min_its:
                continue
            if M.covered_bases(merged_t, s, e) == e - s:
                continue
            near_p = gaps.called_before(s) - head_called <= tolerance
            near_q = tail_called - gaps.called_before(e) <= tolerance
            if e - s < min_len and ((near_p and not has_p) or (near_q and not has_q)):
                # src/teloscope.cpp:465-467 drops end-adjacent arrays shorter than -l (REG-017)
                continue
            cov = M.covered_bases(all_iv, s, e)
            need = (e - s) - motif - match_dist
            rec.check("ORA-08-its-recall", f"{subject}:{header}:{s}-{e}", cov >= need,
                      f"pure {o} tandem cluster of {e - s} bp, reported blocks cover {cov}, "
                      f"need {need}")


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
        rows, summary, _ = M.read_report(p.stdout.decode("utf-8", "replace"), is_text=True)
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

        # terminal calls must not depend on scan mode: strip the full-scan flags for the ultra run
        if any(x in FULL_SCAN_FLAGS for x in base_flags):
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
            ub = {(b["chrom"], b["start"], b["end"])
                  for b in M.read_block_bed(outputs(ultra_dir, stem)["terminal"])}
            fb = {(b["chrom"], b["start"], b["end"])
                  for b in M.read_block_bed(outputs(full_dir, stem)["terminal"])}
            rec.check("MET-04-scan-mode-blocks", subject, ub == fb,
                      f"terminal blocks differ: only ultra {sorted(ub - fb)[:3]}, "
                      f"only full {sorted(fb - ub)[:3]}")

        # -x is monotonic at the match level, not at the block level (REG-012)
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
    finally:
        for d in dirs:
            shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------- driver

def load_waivers(path):
    out = {}
    if not path.exists():
        return out
    for raw in path.read_text().splitlines():
        if raw.startswith("#") or not raw.strip():
            continue
        f = raw.split("\t")
        if len(f) < 4:
            raise SystemExit(f"{path}: malformed waiver row: {raw!r}")
        out[(f[0], f[1])] = f[2]
    return out


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

    waivers = load_waivers(WAIVERS)
    if REGISTER.exists():
        import re as _re
        known = set(_re.findall(r"REG-\d{3}", REGISTER.read_text()))
        unknown = {v for v in waivers.values() if v not in known}
        if unknown:
            raise SystemExit(f"{WAIVERS} cites unknown register ids: {sorted(unknown)}")

    rec = Recorder()
    subjects = list(manifest_runs())
    if not subjects:
        raise SystemExit(
            f"no subjects: {MANIFEST} is missing or empty. Run testFiles/"
            "generate_synthetic.sh. A harness with nothing to check must not report PASS.")
    if not args.quick:
        for g in ["bTaeGut7_chr33_mat.fa.gz", "bTaeGut7_chr33_pat.fa.gz"]:
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
            fasta = M.read_fasta(fasta_path)
            terminal = M.read_block_bed(o["terminal"])
            interstitial = M.read_block_bed(o["interstitial"])
            gaps = M.read_gaps_bed(o["gaps"])
            rows, summary, params = M.read_report(
                proc.stdout.decode("utf-8", "replace"), is_text=True)
            if o["report"].exists():
                _, _, params = M.read_report(o["report"])

            check_derivations(rec, subject, terminal, interstitial, gaps, rows, summary, fasta)
            check_oracle(rec, subject, terminal, interstitial, gaps, fasta,
                         params_from_flags(params, base_flags))
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

    def waiver_for(r):
        return (waivers.get((r["invariant"], r["subject"]))
                or waivers.get((r["invariant"], "*")))

    failures = [r for r in rows if not r["pass"] and not waiver_for(r)]
    waived = [r for r in rows if not r["pass"] and waiver_for(r)]
    # A waiver that excuses nothing has to go with the fix that made it unnecessary.
    used = {(r["invariant"], r["subject"]) for r in rows if not r["pass"]}
    used |= {(r["invariant"], "*") for r in rows if not r["pass"]}
    stale_waivers = sorted(set(waivers) - used)

    print(f"binary:   {TELOSCOPE}")
    print(f"subjects: {len(subjects)}")
    for fam in sorted(families):
        total, bad = families[fam]
        print(f"  {fam}-*  {total - bad}/{total} pass")
    if waived:
        print(f"waived: {len(waived)}")
    if stale_waivers:
        print(f"\n{len(stale_waivers)} stale waiver(s): these no longer excuse anything and")
        print(f"must be removed from {WAIVERS} together with the fix that made them pass:")
        for inv, subj in stale_waivers:
            print(f"  {inv}  {subj}")
        return 1
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
