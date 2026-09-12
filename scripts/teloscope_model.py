#!/usr/bin/env python3
"""Independent parsers and re-derivations for teloscope outputs, computed from the FASTA
or the tool's own primitive outputs so a check can disagree with teloscope.

Column numbering follows docs/outputs.md, 1-based in the prose and 0-based here.
"""

import gzip
import itertools
import re

# ---------------------------------------------------------------- block BED

BED_FIELDS = [
    "chrom", "start", "end", "blockLen", "label", "forwardCount", "reverseCount",
    "canonicalCount", "nonCanonicalCount", "chromSize", "blockType", "arm",
    "gapStatus", "fwdCan", "revCan", "fwdNonCan", "revNonCan", "status",
]
INT_FIELDS = {"start", "end", "blockLen", "forwardCount", "reverseCount", "canonicalCount",
              "nonCanonicalCount", "chromSize", "fwdCan", "revCan", "fwdNonCan", "revNonCan"}


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
    if not path or not str(path) or not _exists(path):
        return []
    out = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip() or s.lstrip().startswith(("#", "track", "browser")):
                continue
            out.append(s)
    return out


def _exists(path):
    import os
    return os.path.exists(path)


# ---------------------------------------------------------------- report

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


# ---------------------------------------------------------------- FASTA ground truth

def read_fasta(path):
    """-> list of (header, sequence). Header is the first whitespace token after '>'."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records, header, chunks = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif header is not None:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


GAP_CHARS = "NnXx"


def gap_runs(seq):
    """Maximal runs of N/n/X/x as half-open [start, end)."""
    return [(m.start(), m.end()) for m in re.finditer(f"[{GAP_CHARS}]+", seq)]


def called_offset_from_start(seq, pos):
    """Number of non-gap bases before pos."""
    return sum(1 for c in seq[:pos] if c not in GAP_CHARS)


def called_offset_from_end(seq, pos):
    return sum(1 for c in seq[pos:] if c not in GAP_CHARS)


# ---------------------------------------------------------------- motif oracle

def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))[::-1]


def canonical_pair(canonical="TTAGGG"):
    """Smaller of motif/revcomp lexicographically is forward (docs/outputs.md:102)."""
    rc = revcomp(canonical)
    return (canonical, rc) if canonical < rc else (rc, canonical)


def substitution_variants(motif, distance):
    """Every substitution variant within `distance` of motif, motif included."""
    out = {motif}
    alphabet = "ACGT"
    for d in range(1, distance + 1):
        for positions in itertools.combinations(range(len(motif)), d):
            for repl in itertools.product(alphabet, repeat=d):
                cand = list(motif)
                for p, c in zip(positions, repl):
                    cand[p] = c
                out.add("".join(cand))
    return out


def find_occurrences(seq, pattern):
    """All start offsets of pattern in seq, overlapping allowed."""
    out, i = [], seq.find(pattern)
    while i != -1:
        out.append(i)
        i = seq.find(pattern, i + 1)
    return out


def locate_matches(seq, canonical="TTAGGG", edit_distance=1):
    """Locate motif matches -> (canonical_intervals, all_intervals), sorted half-open."""
    seq = seq.upper()
    fwd, rev = canonical_pair(canonical)
    canon = {fwd, rev}
    patterns = set()
    for m in canon:
        patterns |= substitution_variants(m, edit_distance)

    canon_iv, all_iv = [], []
    for pat in patterns:
        for s in find_occurrences(seq, pat):
            iv = (s, s + len(pat))
            all_iv.append(iv)
            if pat in canon:
                canon_iv.append(iv)
    return sorted(set(canon_iv)), sorted(set(all_iv))


def merge_intervals(intervals):
    """Union of intervals; teloscope's coverage is a union too (getCoverRuns), not a sum."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def covered_bases(merged, start, end):
    return sum(max(0, min(e, end) - max(s, start)) for s, e in merged)


# ---------------------------------------------------------------- documented rules

FORWARD_LABEL_THRESHOLD = 666
REVERSE_LABEL_THRESHOLD = 333
LABEL_SCALE = 1000


def strand_label(forward_count, total_count):
    """Strand label from counts, integer arithmetic as in include/teloscope.h:236-242."""
    if total_count == 0:
        return "b"
    scaled = forward_count * LABEL_SCALE
    if scaled > total_count * FORWARD_LABEL_THRESHOLD:
        return "p"
    if scaled < total_count * REVERSE_LABEL_THRESHOLD:
        return "q"
    return "b"


def elect_longest(blocks):
    """Arm election by canonical count (src/teloscope.cpp:515-528; REG-005)."""
    # counters start at 0, not -1 (src/teloscope.cpp:512-513), matching the binary
    longest_p = longest_q = None
    max_p = max_q = 0
    for i, b in enumerate(blocks):
        canonical = b["fwdCan"] + b["revCan"]
        if b["arm"] == "p" and canonical > max_p:
            longest_p, max_p = i, canonical
        elif b["arm"] == "q" and canonical >= max_q:
            longest_q, max_q = i, canonical
    return longest_p, longest_q


def granular(blocks):
    """docs/classification.md:64-72. One token per terminal block, ascending by start."""
    blocks = sorted(blocks, key=lambda b: b["start"])
    lp, lq = elect_longest(blocks)
    out = []
    for i, b in enumerate(blocks):
        ch = b["arm"]
        if i in (lp, lq):
            ch = ch.upper()
        if b["label"] == "b":
            ch += "~"
        elif b["label"] != b["arm"]:
            ch += "*"
        out.append(ch)
    return "".join(out)


def scaffold_type(blocks):
    """docs/classification.md:41-43: both arms t2t, one incomplete, neither none."""
    lp, lq = elect_longest(sorted(blocks, key=lambda b: b["start"]))
    n = (lp is not None) + (lq is not None)
    return {2: "t2t", 1: "incomplete", 0: "none"}[n]


ANOMALY_ORDER = ["discordant_p", "discordant_q", "balanced_p", "balanced_q", "misassembly"]


def anomalies(blocks):
    """Anomaly flags per docs/classification.md:43-45; balanced is tested before discordant."""
    blocks = sorted(blocks, key=lambda b: b["start"])
    lp, lq = elect_longest(blocks)
    flags = set()
    for idx, suffix in ((lp, "p"), (lq, "q")):
        if idx is None:
            continue
        b = blocks[idx]
        if b["label"] == "b":
            flags.add(f"balanced_{suffix}")
        elif b["label"] != b["arm"]:
            flags.add(f"discordant_{suffix}")
    for i, _ in enumerate(blocks):
        if i not in (lp, lq):
            flags.add("misassembly")
            break
    ordered = [f for f in ANOMALY_ORDER if f in flags]
    return ",".join(ordered) if ordered else "."


def labels_column(blocks):
    """docs/outputs.md: lowercase arm chars of the elected blocks, or literal 'none'."""
    blocks = sorted(blocks, key=lambda b: b["start"])
    lp, lq = elect_longest(blocks)
    out = "".join(blocks[i]["arm"] for i in (lp, lq) if i is not None)
    return out or "none"


def n50(lengths):
    """include/teloscope.h:244-255: sort descending, first length where cum*2 >= total."""
    if not lengths:
        return 0
    total = sum(lengths)
    cum = 0
    for length in sorted(lengths, reverse=True):
        cum += length
        if cum * 2 >= total:
            return length
    return 0


def contig_lengths(seq_len, gaps):
    """The runs between gaps. A record with no gaps is one contig of its full length."""
    out, prev = [], 0
    for s, e in gaps:
        if s > prev:
            out.append(s - prev)
        prev = e
    if seq_len > prev:
        out.append(seq_len - prev)
    return out


# ---------------------------------------------------------------- recall oracle

def canonical_positions(seq, canonical="TTAGGG"):
    """Sorted exact canonical occurrences (start, end, orientation); 'p' is the smaller motif."""
    seq = seq.upper()
    fwd, rev = canonical_pair(canonical)
    out = [(s, s + len(fwd), "p") for s in find_occurrences(seq, fwd)]
    out += [(s, s + len(rev), "q") for s in find_occurrences(seq, rev)]
    return sorted(out)


class GapIndex:
    """Gap runs of one record, indexed once for O(log n) called-base arithmetic."""

    def __init__(self, gaps):
        self.gaps = list(gaps)
        self.starts = [s for s, _ in self.gaps]
        self.prefix, total = [], 0
        for s, e in self.gaps:
            self.prefix.append(total)
            total += e - s
        self.total = total

    def called_before(self, pos):
        """Called (non-gap) bases in [0, pos)."""
        import bisect
        i = bisect.bisect_left(self.starts, pos)
        gap_bases = self.prefix[i] if i < len(self.prefix) else self.total
        if i > 0 and self.gaps[i - 1][1] > pos:
            gap_bases -= self.gaps[i - 1][1] - pos
        return pos - gap_bases

    def has_gap_between(self, a, b):
        """True when any gap run overlaps [a, b)."""
        import bisect
        i = bisect.bisect_left(self.starts, b)
        return i > 0 and self.gaps[i - 1][1] > a


def dense_windows(matches, gaps, min_len, density, max_dist, min_counts):
    """Same-orientation chains within max_dist, shrunk to meet -l, -y and the count."""
    groups, cur = [], []
    for s, e, o in matches:
        if cur and (o != cur[-1][2] or s - cur[-1][1] > max_dist
                    or gaps.has_gap_between(cur[-1][1], s)):
            groups.append(cur)
            cur = []
        cur.append((s, e, o))
    if cur:
        groups.append(cur)

    out = []

    def shrink(g):
        merged = merge_intervals([(s, e) for s, e, _ in g])
        lo, hi = 0, len(g) - 1
        while hi - lo + 1 >= min_counts:
            span = g[hi][1] - g[lo][0]
            if span < min_len:
                return
            covered = covered_bases(merged, g[lo][0], g[hi][1])
            if covered >= density * span:
                out.append((g[lo][0], g[hi][1], g[lo][2], hi - lo + 1))
                shrink(g[:lo]) if lo > 0 else None
                shrink(g[hi + 1:]) if hi + 1 < len(g) else None
                return
            # drop the end match sitting behind the wider stretch of non-telomeric sequence
            if g[lo + 1][0] - g[lo][1] >= g[hi][0] - g[hi - 1][1]:
                lo += 1
            else:
                hi -= 1

    for g in groups:
        if g:
            shrink(g)
    return sorted(out)


def pure_tandem_runs(seq, canonical="TTAGGG", min_repeats=4):
    """Runs of at least min_repeats canonical motifs -> sorted (start, end, orientation)."""
    seq = seq.upper()
    fwd, rev = canonical_pair(canonical)
    out = []
    for motif, o in ((fwd, "p"), (rev, "q")):
        for m in re.finditer(f"(?:{motif}){{{min_repeats},}}", seq):
            out.append((m.start(), m.end(), o))
    return sorted(out)


def merge_runs(runs, gaps, max_dist, density):
    """Merge same-orientation runs within max_dist and no gap while coverage stays >= density."""
    out, covered = [], 0
    for s, e, o in runs:
        if out and out[-1][2] == o and s - out[-1][1] <= max_dist \
                and not gaps.has_gap_between(out[-1][1], s) \
                and covered + (e - s) >= density * (e - out[-1][0]):
            out[-1] = (out[-1][0], e, o)
            covered += e - s
        else:
            out.append((s, e, o))
            covered = e - s
    return out
