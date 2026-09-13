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


def block_spans(path):
    """(chrom, start, end) of every row in a block BED, read fresh from disk."""
    return {(b["chrom"], b["start"], b["end"]) for b in M.read_block_bed(path)}


# ---------------------------------------------------------------- families

def contig_bounds(gaps_sorted, pos, seq_len):
    """[start, end) of the contig containing pos, from a sorted list of gap runs."""
    cs = max((ge for gs, ge in gaps_sorted if ge <= pos), default=0)
    ce = min((gs for gs, ge in gaps_sorted if gs >= pos), default=seq_len)
    return cs, ce


def check_one_arm_per_end(rec, subject, terminal, gaps_by_chrom, lengths, manual_curation):
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
        gaps_sorted = sorted(gaps_by_chrom.get(chrom, []))
        by_contig = {}
        for b in blocks:
            bounds = contig_bounds(gaps_sorted, b["start"], lengths[chrom])
            by_contig.setdefault(bounds, []).append(b)
        for bounds, cblocks in by_contig.items():
            for end in ("p", "q"):
                n = sum(1 for b in cblocks if b["closestEnd"] == end)
                rec.check("DER-19-one-arm-per-end", f"{subject}:{chrom}:{bounds}:{end}", n <= 1,
                          f"{n} terminal rows with closestEnd {end} on contig {bounds}")


def check_junction_class(rec, subject, terminal, interstitial, gaps_by_chrom, link_distance):
    """Recompute each interstitial row's junction class from its nearest same-contig neighbour."""
    by_chrom = {}
    for b in terminal + interstitial:
        by_chrom.setdefault(b["chrom"], []).append(b)
    its_ids = {id(b) for b in interstitial}
    for chrom, blocks in by_chrom.items():
        ordered = sorted(blocks, key=lambda b: b["start"])
        gap_index = M.GapIndex(sorted(gaps_by_chrom.get(chrom, [])))
        for i, b in enumerate(ordered):
            if id(b) not in its_ids:
                continue
            candidates = []
            if i > 0 and not gap_index.has_gap_between(ordered[i - 1]["end"], b["start"]):
                candidates.append((b["start"] - ordered[i - 1]["end"], ordered[i - 1], b))
            if i + 1 < len(ordered) and not gap_index.has_gap_between(b["end"], ordered[i + 1]["start"]):
                candidates.append((ordered[i + 1]["start"] - b["end"], b, ordered[i + 1]))
            candidates = [c for c in candidates if c[0] <= link_distance]
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


def check_derivations(rec, subject, terminal, interstitial, gaps, rows, summary, fasta, params):
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

    threshold = float(params.get("label_threshold", 0.667))
    for b in terminal + interstitial:
        where = f"{subject}:{b['chrom']}:{b['start']}"
        forward = b["fwdCan"] + b["fwdNonCan"]
        total = forward + b["revCan"] + b["revNonCan"]
        rec.eq("DER-01-strand-label-from-counts", where, b["teloLabel"],
               M.strand_label(forward, total, threshold),
               f"(fwd {forward} of {total})")
    for b in terminal:
        where = f"{subject}:{b['chrom']}:{b['start']}"
        rec.check("DER-01-terminal-label-is-strand-pure", where, b["teloLabel"] != "b",
                  "a terminal row's teloLabel must be p or q, never b (R1)")

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

    manual_curation = params.get("manual_curation") == "true"
    check_one_arm_per_end(rec, subject, terminal, gaps_by_chrom, lengths, manual_curation)
    link_distance = int(params.get("link_distance", 1000))
    check_junction_class(rec, subject, terminal, interstitial, gaps_by_chrom, link_distance)

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
               "--link-distance": "link_distance",
               "--label-threshold": "label_threshold",
               "--min-block-counts": "min_block_counts"}
FULL_SCAN_FLAGS = {"-r", "-g", "-e", "-m", "-i", "-n", "--out-win-repeats", "--out-gc",
                   "--out-entropy", "--out-matches", "--manual-curation"}


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
    # zone = min(tolerance, -t) (R8): -t caps the start zone, never the extent
    tolerance = min(int(params.get("terminal_tolerance", 3000)),
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
        # a fragmented row's span includes an unheld --link-distance gap; gate on teloLen instead
        called = b["teloLen"] if canonical_only and b["teloLen"] < b["end"] - b["start"] \
            else b["end"] - b["start"]
        rec.check("ORA-02-block-meets-canonical-density", where,
                  called <= 0 or covered >= density * called,
                  f"independently measured {covered} "
                  f"{'canonical' if canonical_only else 'matched'} bases over {called} "
                  f"bases = {covered / called if called else 0:.3f}, -y is {density}")

    manual_curation = params.get("manual_curation") == "true"
    for b in terminal:
        seq = seqs.get(b["chrom"])
        if seq is None:
            continue
        where = f"{subject}:{b['chrom']}:{b['start']}"
        rec.check("DER-17-terminal-length-floor", where, b["teloLen"] >= min_len,
                  f"teloLen {b['teloLen']} against -l {min_len}")
        rec.check("DER-18-contig-row-needs-manual-curation", where,
                  b["teloType"] != "contig" or manual_curation,
                  "a contig row appeared without -n/--manual-curation")
        # tolerance counts called bases from its OWN contig's end named by closestEnd (R4.3)
        cs, ce = contig_bounds(M.gap_runs(seq), b["start"], len(seq)) \
            if b["teloType"] == "contig" else (0, len(seq))
        from_start = M.called_offset_from_start(seq, b["start"]) - M.called_offset_from_start(seq, cs)
        from_end = M.called_offset_from_end(seq, b["end"]) - M.called_offset_from_end(seq, ce)
        dist = from_start if b["closestEnd"] == "p" else from_end
        rec.check("ORA-03-terminal-block-is-terminal", where, dist <= tolerance,
                  f"{dist} called bases from the {b['closestEnd']} end of its contig, "
                  f"--terminal-tolerance is {tolerance}")

    check_recall(rec, subject, terminal, interstitial, fasta, params)

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
    """Verify reported blocks cover dense canonical FASTA windows within link+d+k+motif slack."""
    canonical = params.get("canonical", "CCCTAA/TTAGGG").split("/")[-1]
    motif = len(canonical)
    min_counts = int(params.get("min_block_counts", 2))
    min_len = int(params.get("min_block_len", 300))
    density = float(params.get("min_block_density", 0.5))
    max_dist = int(params.get("max_block_dist", 500))
    match_dist = int(params.get("max_match_dist", 50))
    link_distance = int(params.get("link_distance", 1000))
    terminal_limit = int(params.get("terminal_limit", 50000))
    tolerance = min(int(params.get("terminal_tolerance", 3000)), terminal_limit)
    full_scan = params.get("ultra_fast", "true") == "false"
    manual_curation = params.get("manual_curation") == "true"
    # pieces up to link_distance apart chain into one row (R6's D), so recall needs that slack too
    slack = link_distance + max_dist + match_dist + motif

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

        # No extent cap any more (R2): re-validate each window as found, unclipped.
        starts = [s for s, _, _ in matches]

        def validate(w):
            a, b = bisect.bisect_left(starts, w[0]), bisect.bisect_left(starts, w[1])
            inside = [(s, e, o) for s, e, o in matches[a:b] if e <= w[1] and o == w[2]]
            if len(inside) < min_counts:
                return None
            s, e = inside[0][0], inside[-1][1]
            if e - s < min_len or len(inside) * motif < density * (e - s):
                return None
            return (s, e, w[2], len(inside))

        ends, both_ends = {}, set()
        for w in windows:
            # a chain never crosses a contig boundary (R1/R4), however close in called bases
            near_p = (gaps.called_before(w[0]) - head_called <= tolerance
                      and not gaps.has_gap_between(first_base, w[0]))
            near_q = (tail_called - gaps.called_before(w[1]) <= tolerance
                      and not gaps.has_gap_between(w[1], last_base))
            if not (near_p or near_q):
                continue
            if near_p and near_q:
                both_ends.add(w)
            c = validate(w)
            if c is None:
                continue
            if near_p:
                ends.setdefault("p", []).append((c, w))
            if near_q:
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
                label = hit.get("teloLabel")
                rec.check("ORA-06-terminal-recall", where + ":label", label in (outer[2], "b"),
                          f"window is {outer[2]}-oriented, overlapping block is labelled {label}")
                closest_end = hit.get("closestEnd")
                # R3 can send a same-strand chain to its concordant end, not the nearer one
                r3_redirected = closest_end != end and label == closest_end
                if raw not in both_ends and not r3_redirected:
                    rec.check("ORA-06-terminal-recall", where + ":arm", closest_end == end,
                              f"window sits at the {end} end, overlapping block says "
                              f"closestEnd {closest_end}")
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

        # -n: every contig's own two ends count as ends too (R4.3), covered by any terminal row
        if manual_curation:
            bounds, prev = [], 0
            for gs, ge in runs:
                if gs > prev:
                    bounds.append((prev, gs))
                prev = ge
            if n > prev:
                bounds.append((prev, n))
            for cs, ce in bounds:
                c_head, c_tail = gaps.called_before(cs), gaps.called_before(ce)
                for edge, near, key in (
                    ("p", lambda w: gaps.called_before(w[0]) - c_head <= tolerance
                     and not gaps.has_gap_between(cs, w[0]), lambda w: w[0]),
                    ("q", lambda w: c_tail - gaps.called_before(w[1]) <= tolerance
                     and not gaps.has_gap_between(w[1], ce), lambda w: -w[1]),
                ):
                    cand = sorted((w for w in windows if cs <= w[0] and w[1] <= ce and near(w)), key=key)
                    c = validate(cand[0]) if cand else None
                    if c is None:
                        continue
                    rec.check("ORA-06-contig-recall", f"{subject}:{header}:{cs}-{ce}:{edge}",
                              hit_of(c) is not None,
                              f"contig [{cs},{ce}) {edge} end has an uncovered dense window {c}")

        if not full_scan:
            continue

        # Reported blocks must cover every pure tandem cluster (REG-003): no length floor (R5).
        clusters = M.merge_runs(M.pure_tandem_runs(seq, canonical, 4), gaps, max_dist, density)
        for s, e, o in clusters:
            if M.covered_bases(merged_t, s, e) == e - s:
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

        # terminal calls must not depend on scan mode; -n forces full scan by itself (R8), so
        # there is no ultra-fast run to compare it against.
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
            n_terminal = M.read_block_bed(outputs(d8, stem)["terminal"])
            n_set = {(b["chrom"], b["start"], b["end"]) for b in n_terminal}
            no_n_set = block_spans(outputs(d7, stem)["terminal"])
            contig_rows = {(b["chrom"], b["start"], b["end"])
                           for b in n_terminal if b["teloType"] == "contig"}
            rec.check("MET-07-manual-curation", subject,
                      n_set == no_n_set | contig_rows,
                      f"-n terminal BED should equal the -n-less full-scan terminal BED plus "
                      f"contig rows: missing {sorted((no_n_set | contig_rows) - n_set)[:3]}, "
                      f"extra {sorted(n_set - (no_n_set | contig_rows))[:3]}")
            its_set = block_spans(outputs(d8, stem)["interstitial"])
            overlap = contig_rows & its_set
            rec.check("MET-07-manual-curation", subject + ":its", not overlap,
                      f"contig rows {sorted(overlap)[:3]} also appear in the -n run's "
                      f"own interstitial BED")
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
            fasta = M.read_fasta(fasta_path)
            terminal = M.read_block_bed(o["terminal"])
            interstitial = M.read_block_bed(o["interstitial"])
            gaps = M.read_gaps_bed(o["gaps"])
            rows, summary, params = M.read_report(
                proc.stdout.decode("utf-8", "replace"), is_text=True)
            if o["report"].exists():
                _, _, params = M.read_report(o["report"])
            params = params_from_flags(params, base_flags)

            check_derivations(rec, subject, terminal, interstitial, gaps, rows, summary,
                               fasta, params)
            check_oracle(rec, subject, terminal, interstitial, gaps, fasta, params)
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
