[Back to README](index.md)

# Documentation and behaviour conflict register

Where the documentation and the code disagree, the documentation is the specification; the
disagreement is decided here first, never by editing a test to match the binary. Each row
has an id, the claim, the behaviour, a decision, and a status: `open` (undecided), `decided`
(a resolution is chosen and the row names the pull request that carries it), `resolved-doc`
(the documentation was corrected), or `resolved-code` (the code was fixed).
`validateFiles/intent_waivers.tsv`, `invariant_waivers.tsv`, and `intent_blocked.tsv` each
cite a row's id.

Repair pull requests, in order: B1 docs-only resolutions; B2 the REG-010 banner; B3a block
engine (no gap bridging, interstitial gate, guard deletion, REG-012, strand-split
interstitial segments); B3b chained terminal rule; B3c per-contig scanning and `-n`; B4 the
12-column BED, junction classes and the report script; B5 `--label-threshold`. This pull
request carries B1 through B5. B6 (`-a` output and REG-011) remains open.

## Decided, awaiting their pull request

### REG-011 — reverse-oriented path components are silently skipped

`src/input.cpp:1013` has an empty `else` branch for path components in `-` orientation, so
such a segment contributes its length but none of its matches. Pre-existing, and there is
no fixture for it.

**Decision (B6):** resolved-code, with a GFA fixture carrying a reverse-oriented path
component.

## Resolved

### REG-001 — the terminal rule is named after two different flags

**Resolved (doc and code).** `-t` sets how far the fast scan reads and caps the start zone;
`--terminal-tolerance` (now `3000`, called bases) is where a telomere may start — the start
zone is the smaller of the two. Chaining pieces together is its own knob,
`--link-distance` (`1000`). Each of the four distances now does exactly one job.

### REG-002 — two pages give different minimum block lengths

**Resolved (doc).** `docs/algorithm.md` now says `42` bp for read subsets and `300` bp for
assembly annotation, matching the code.

### REG-003 — `minCanonicalCount` is undocumented and unreachable

**Resolved (doc and code).** The four-canonical-match floor for an interstitial row is
documented and written to the `#params` header as `min_canonical_count`; it stays a
built-in constant, not a flag.

### REG-004 — should a head-to-head fusion across a gap count as T2T?

**Resolved (code).** Arrays on either side of an `N` run sit on different contigs, since a
gap is a contig boundary, and never see each other: each is one interstitial row, classed
against its own contig's neighbours only, usually `single`. A true head-to-head pair on one
contig, with no gap between, is two rows with a junction class instead: `fusion` for a
reverse array (`TTAGGG`) then a forward one (`CCCTAA`) along the sequence, `tail_to_tail`
for forward then reverse.

### REG-005 — "longest arm" is elected by canonical count, not by length

**Resolved (code).** There is no election by canonical count any more. The p arm is
whichever chain anchors at the first contig's start, and the q arm whichever chain anchors
at the last contig's end; each qualifies on its own, and a strand that doesn't match its
end is flagged `discordant_p`/`discordant_q` rather than out-voted by a shorter block. When
a single contig's array runs end to end, both chains would claim it; it belongs to the end
its strand points to (forward to the p end, reverse to the q end), and the other end has no
telomere.

### REG-006 — the strand label is not symmetric under reverse complement

**Resolved (code).** The strand threshold is now exact thirds, symmetric by construction,
and tunable with one flag, `--label-threshold` (default `0.667`), used in two places: the
`p`/`b`/`q` label of an interstitial row (`p` above it, `q` below one minus it, `b`
between), and the strand-purity floor a terminal piece's own exact repeats must clear to
qualify at all. `closestEnd` at an exact midpoint reads `p`.

### REG-007 — two definitions of block density, five lines apart

**Resolved (code).** A block never crosses a gap — blocks are built per contig, and a
contig has no `N` inside — so density is always measured over the plain block length, with
no gap to discount. Density itself is coverage, not a count-times-motif-length formula: the
fraction of the block's own bases covered by repeats, canonical coverage for a terminal
piece and all-repeat coverage for an interstitial row.

### REG-008 — threshold comparisons are inclusive, and nothing says so

**Resolved (doc).** The docs state plainly that a piece or row exactly at a threshold
(`-l`, `-y`, `--min-block-counts`) is kept, not rejected.

### REG-009 — retired vocabulary inside the specification

**Resolved (doc).** `docs/classification.md` no longer calls `discordant` a scaffold type;
it is an anomaly flag on an arm.

### REG-010 — two summary headers are missing a space

**Resolved (code).** Both summary banners carry the missing space; the checked-in
snapshots were regenerated to match.

### REG-012 — block filters are not pure post-filters

**Resolved (code).** The old cross-arm guard is gone along with the code that produced it.
Each contig end builds its own chain independently, so a piece that fails `-l` or `-y` is
simply not part of any row — it can no longer loosen a filter elsewhere.

### REG-013 — a terminal array of near-canonical variants can never be called

**Resolved (doc).** Terminal density counts canonical coverage only, as it has since
v0.1.5; `-x`-derived variants affect the match BED files and interstitial rows, not
terminal density.

### REG-014 — abutting opposite orientations are split, never pooled as one balanced array

**Resolved (code).** The outermost qualifying match anchors the telomere; same-strand
pieces within `--link-distance` chain into one row (`fragmented_p`/`fragmented_q` when
there is more than one piece). The first real array of the other strand ends the chain and
becomes an interstitial row; its junction class follows the order of the two arrays along
the sequence, `fusion` or `tail_to_tail`, not always the same one. `balanced` is retired for
terminal rows — it survives only as the interstitial `b` label for an array that
alternates strands.

### REG-015 — `-n/--manual-curation` and `-a/--out-fasta` do nothing

**Resolved (code), in part.** Fast mode never looks at internal contig ends. With `-i` (the
full scan), an internal array is an ordinary interstitial row. `-n/--manual-curation`
implies the full scan and additionally builds a telomere chain at every contig end, writing
the internal ones to the terminal BED as `teloType=contig` rows instead of only accepting
the flag for compatibility. `-a/--out-fasta` still writes nothing; that half of this row
stays open for B6.

### REG-016 — the report script reads strand where it means arm

**Resolved (code).** `scripts/teloscope_report.py` reads `closestEnd` for arm decisions,
not `teloLabel`.

### REG-017 — a short end-adjacent canonical array vanishes

**Resolved (code).** The guard that dropped a short end-adjacent canonical array is
deleted. The array is now an interstitial row, or the outermost array at its end.

### REG-018 — sequence beyond `-t` inside an accepted seed is reported nowhere

**Resolved (code).** `-t` no longer caps a telomere's extent. It only caps the start zone
and the fast-mode window; a telomere reaches as far as its pieces chain, whatever `-t` is
set to.
