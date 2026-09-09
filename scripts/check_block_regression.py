#!/usr/bin/env python3
"""Extract per-scaffold block calls from a teloscope build, and compare two extractions.

extract  runs the binary over the .tst fixtures and writes one sorted row per scaffold
compare  joins two extractions, prints a scaffold-type transition matrix, and gates on a roster
"""
import argparse, os, re, shlex, subprocess, sys, tempfile

FIELDS = ["fixture", "scaffold", "scaffold_type", "telomere_count",
          "granular_label", "its_count", "terminal_blocks", "interstitial_blocks"]


def granular_tokens(label):
    """Split a granular label into one token per block: a letter plus an optional star."""
    out, i = [], 0
    while i < len(label):
        tok = label[i]
        i += 1
        if i < len(label) and label[i] == "*":
            tok += "*"
            i += 1
        out.append(tok)
    return out


def parse_report(path):
    """Return (ultra_fast, [row dicts]). Columns are read by name, never by position."""
    ultra, names, rows = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#params"):
                m = re.search(r"ultra_fast=(\w+)", line)
                if m:
                    ultra = m.group(1) == "true"
            elif line.startswith("#columns"):
                names = line.split("\t")[1:]
            elif not line or line.startswith("#") or line.startswith("+++"):
                if line.startswith("+++"):
                    break
            elif names and line.split("\t")[0] == "pos":
                continue
            elif names:
                parts = line.split("\t")
                if len(parts) < len(names):
                    continue
                rows.append(dict(zip(names, parts)))
    return ultra, rows


def parse_blocks(path, has_tag):
    """chrom -> list of (start, end, label, tag), sorted by start."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            tag = p[9] if has_tag and len(p) > 9 else "."
            out.setdefault(p[0], []).append((int(p[1]), int(p[2]), p[3], tag))
    for v in out.values():
        v.sort()
    return out


def extract(binary, fixture_dir, out_path):
    rows, skipped = [], 0
    for name in sorted(os.listdir(fixture_dir)):
        if not name.endswith(".tst"):
            continue
        with open(os.path.join(fixture_dir, name)) as fh:
            cmd = fh.readline().strip()
        if not cmd:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            cmd = cmd.replace("%OUTDIR%", tmp)
            argv = shlex.split(cmd)
            if "-o" not in argv and "--output" not in argv:
                argv += ["-o", tmp]
            if "-n" not in argv and "--manual-curation" not in argv:
                argv += ["-n"]          # keeps BED rows aligned with the granular label
            try:
                subprocess.run([binary] + argv, cwd=os.getcwd(), timeout=300,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.TimeoutExpired:
                skipped += 1
                continue
            reports = [f for f in os.listdir(tmp) if f.endswith("_report.tsv")]
            if not reports:
                skipped += 1        # gfa, read-subset and error-exit fixtures
                continue
            base = reports[0][: -len("_report.tsv")]
            ultra, rrows = parse_report(os.path.join(tmp, reports[0]))
            term = parse_blocks(os.path.join(tmp, base + "_terminal_telomeres.bed"), True)
            inter = parse_blocks(os.path.join(tmp, base + "_interstitial_telomeres.bed"), False)
            for r in rrows:
                chrom = r.get("header", "")
                toks = granular_tokens(r.get("granular", ""))
                tb = term.get(chrom, [])
                tstr = ";".join(
                    "%d-%d:%s:%s" % (s, e, toks[i] if i < len(toks) else lab, tag)
                    for i, (s, e, lab, tag) in enumerate(tb)) or "-"
                istr = ("NA" if ultra else
                        ";".join("%d-%d:%s" % (s, e, lab)
                                 for s, e, lab, _ in inter.get(chrom, [])) or "-")
                rows.append([name, chrom, r.get("type", ""), r.get("telomeres", ""),
                             r.get("granular", "") or "-",
                             "NA" if ultra else r.get("its", ""), tstr, istr])
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, "w") as fh:
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("extracted %d scaffold rows, skipped %d fixtures with no report" % (len(rows), skipped))


def load(path):
    with open(path) as fh:
        names = fh.readline().rstrip("\n").split("\t")
        return {(r[0], r[1]): dict(zip(names, r))
                for r in (l.rstrip("\n").split("\t") for l in fh) if len(r) == len(names)}


def load_roster(path):
    roster = {}
    if not path or not os.path.exists(path):
        return roster
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                roster[(p[0], p[1])] = (p[2], p[3])
    return roster


def compare(old_path, new_path, roster_path):
    old, new, roster = load(old_path), load(new_path), load_roster(roster_path)
    matrix, changed, failures = {}, [], []
    for key in sorted(set(old) | set(new)):
        o, n = old.get(key), new.get(key)
        if o and n:
            matrix.setdefault((o["scaffold_type"], n["scaffold_type"]), 0)
            matrix[(o["scaffold_type"], n["scaffold_type"])] += 1
        if o == n:
            continue
        if not o or not n:
            failures.append((key, "removed" if not n else "added", "", ""))
            continue
        diff = [f for f in FIELDS[2:] if o[f] != n[f]]
        changed.append((key, o["scaffold_type"], n["scaffold_type"], ",".join(diff)))
        allowed = roster.get(key)
        if not (allowed and allowed == (o["scaffold_type"], n["scaffold_type"])):
            failures.append((key, o["scaffold_type"], n["scaffold_type"], ",".join(diff)))

    print("\n=== scaffold type transitions ===")
    off = 0
    for (a, b), c in sorted(matrix.items()):
        mark = "" if a == b else "   <-- off diagonal"
        off += 0 if a == b else c
        print("  %-22s -> %-22s %5d%s" % (a, b, c, mark))
    print("  off-diagonal total: %d" % off)

    if changed:
        print("\n=== changed scaffolds ===")
        for (fx, sc), a, b, d in changed:
            print("  %s  %s  %s -> %s  [%s]" % (fx, sc, a, b, d))

    if failures:
        print("\n=== FAIL: %d change(s) not on the roster ===" % len(failures))
        for (fx, sc), a, b, d in failures:
            print("  %s\t%s\t%s\t%s\t<reason>\t<pr>" % (fx, sc, a, b))
        print("\nAdd a roster line only with a real reason naming a plan section.")
        return 1
    print("\nPASS: no unrostered changes (%d scaffolds compared)" % len(set(old) & set(new)))
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract"); e.add_argument("--binary", required=True)
    e.add_argument("--fixtures", default="validateFiles"); e.add_argument("--out", required=True)
    c = sub.add_parser("compare"); c.add_argument("--old", required=True)
    c.add_argument("--new", required=True); c.add_argument("--roster", default="")
    a = ap.parse_args()
    if a.cmd == "extract":
        extract(a.binary, a.fixtures, a.out); return 0
    return compare(a.old, a.new, a.roster)


if __name__ == "__main__":
    sys.exit(main())
