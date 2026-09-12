[Back to README](index.md)

# Documentation and behaviour conflict register

Where the documentation and the code disagree, the documentation is the specification; the
disagreement is decided here first, never by editing a test to match the binary. Each row
has an id, the claim, the behaviour, a decision, and a status: `open` (undecided), `decided`
(a resolution is chosen and the row names the pull request that carries it), `resolved-doc`
(the documentation was corrected), or `resolved-code` (the code was fixed).
`validateFiles/intent_waivers.tsv`, `invariant_waivers.tsv`, and `intent_blocked.tsv` each
cite a row's id. `ORA-02` measures density the way the current code does until B3a lands
(REG-007).

Repair pull requests, in order: B1 docs-only resolutions; B2 the REG-010 banner; B3a block
engine (no gap bridging, interstitial gate, guard deletion, REG-012, strand-split
interstitial segments); B3b chained terminal rule; B3c per-contig scanning and `-n`; B4 the
12-column BED, junction classes and the report script; B5 `--label-threshold`; B6 `-a`
output and REG-011.

## Decided, awaiting their pull request

### REG-013 — a terminal array of near-canonical variants can never be called

`docs/parameters.md:61` documents `-x/--edit-distance` as "allowed mismatches per repeat
unit", and `docs/algorithm.md:15` says pattern expansion includes "allowed edit-distance
variants". Both imply a terminal array built from one-substitution repeats is detectable
at `-x 1`.

It is not: the terminal density gate at `src/teloscope.cpp:403` measures
`getCoveredBases(canRuns, ...)` — **canonical** coverage only — against the whole block
length. An array of pure variants has zero canonical coverage, so it fails for any
`-y > 0`. Verified on `testFiles/edit_test.fa`, whose p arm is 100 × `CTCTAA`: all 100
matches are found and written to `*_noncanonical_matches.bed`, and the block is still
absent from the terminal BED at `-y 0.5`, `-y 0.1` and `-y 0.01`.

**Decision (B1):** resolved-doc. Terminal density counts canonical coverage only, as it has
since v0.1.5; `-x` is documented as affecting the match BED files and interstitial blocks,
not terminal density. The `edit_test` expectation becomes `none`.

*Fixtures:* `edit_test` · *Test:* `scripts/test_synthetic_intent.py`

### REG-005 — "longest arm" is elected by canonical count, not by length

`docs/classification.md:41-45` says completeness is decided by "both longest arms" and
`:72` says "uppercase marks the longest block for that arm". The election at
`src/teloscope.cpp:515-528` compares `fwdCanCount + revCanCount`, i.e. the number of exact
canonical matches, not `blockLen`. A long variant-rich block therefore loses to a short
canonical one, and the report calls the shorter block the arm.

Demonstrated by `mo_elect_longer_loses`, built so the two rules disagree: a 720 bp outer
array carrying 60 canonical matches against a 600 bp inner one carrying 100. The shorter
array wins the election. Under the documented rule the granular string is `Pp`; the code
emits `pP`.

**Decision (B1):** resolved-doc. The arm is the block with the most canonical repeats, as
`docs/classification.md:72` already says.

*Fixtures:* `mo_elect_longer_loses` · *Test:* `scripts/test_synthetic_intent.py`

### REG-014 — abutting opposite orientations are split, never pooled as one balanced array

`docs/classification.md:28` describes `balanced_p`/`balanced_q` as "that arm carries both
orientations in similar proportion, which is what a fusion or a collapsed repeat looks
like", which reads as though a head-to-head array at one end produces one balanced block.

It does not. The strand-inversion cut in `getSeeds` (`inversionRun = 3`,
`src/teloscope.cpp:121`) splits a run of three opposite-orientation matches, so
`mirror_both_start.fa` and `mirror_both_end.fa` yield two separate terminal blocks and are
reported as `misassembly`, not `balanced_p`/`balanced_q`. A balanced label arises only from
*interleaved* orientations within one seed, as in `balanced.fa`.

**Decision (B3b):** resolved-code. The outermost array anchors the telomere; same-orientation
pieces within the linking distance chain into one row whose `teloLen` is the sum of the
pieces (`fragmented_p`/`fragmented_q` replaces `misassembly`); the first qualifying
opposite-orientation array ends the chain and is an interstitial row with junction class
`tail_to_tail`. Balanced stays reserved for interleaved orientations inside one block, and
the docs say so. The linking distance is `--terminal-tolerance`: one knob for how far a
telomeric array may sit from a record end or from the previous array and still belong to
the same structure. Fast mode writes every interstitial block its tip scan finds, with its
junction class, so the array is never lost without `-i`.

*Fixtures:* `mirror_both_start`, `mirror_both_end`, `balanced`, `jn_p_then_q_k3`, `jn_p_then_q_k4` · *Test:* `scripts/test_synthetic_intent.py`

### REG-001 — the terminal rule is named after two different flags

`docs/classification.md:7` says a block counts as scaffold-terminal "when it falls within
`-t/--terminal-limit` of the sequence start or end". `docs/outputs.md:104` says it "has to
start within `--terminal-tolerance` called bases of a record end". These are different
flags with different defaults, 50000 and 2000.

The code uses both, for different things: `--terminal-tolerance` gates where a seed may
*start*, and `-t/--terminal-limit` bounds the block's *extent* and the ultra-fast scan
window, with tolerance capped at the limit. Both pages are incomplete rather than one
being wrong.

**Decision (B1, B3b):** resolved-doc for the two passages, resolved-code for the default.
Both flags stay: `-t` caps the extent and the fast-mode scan window, `--terminal-tolerance`
is where a block may start, in called bases. `docs/classification.md:7` is rewritten to say
both, because zebra finch telomeres reach 26 kb and one knob cannot serve both roles. The
tolerance's original name was "distance to end" (`d2end`) with a default of 3000, not 2000;
B3b restores the 3000 default, and the docs call the flag by both names. `-d` is what merges
small seeds into blocks, and 500 is the smallest value that still chains two minimal 300 bp
blocks (2×300 + 500 keeps repeat coverage above 0.5). `--terminal-tolerance` also serves as
the linking distance between pieces of one telomere and between abutting blocks (REG-014):
one knob, one meaning, "how far a telomeric array may sit from a record end or from the
previous array and still belong to the same structure".

### REG-002 — two pages give different minimum block lengths

`docs/algorithm.md:46` says read-subset mode defaults to 60 bp and assembly annotation to
500 bp. `docs/parameters.md:84` says 42 and 300. `include/input.h:44` says 300, and
`src/read-filter.cpp:12-14` says 42. `docs/algorithm.md` is stale.

**Decision (B1):** resolved-doc. `docs/algorithm.md:46` says 42/300.

### REG-003 — `minCanonicalCount` is undocumented and unreachable

`src/teloscope.cpp:90` declares `constexpr uint32_t minCanonicalCount = 4;` and applies it
at `:475` to every interstitial call. It appears in no documentation page, is not on the
command line, and is not written to the `#params` header, so a run cannot be reproduced
from its own metadata.

**Decision (B1, B3a):** resolved-doc for the docs, resolved-code for the `#params` header.
The four-canonical-repeat rule is documented and written to the header; it stays a
constant, not a flag.

### REG-004 — should a head-to-head fusion across a gap count as T2T?

`its_gap_headtohead.fa` and `its_headtohead.fa` place reverse repeats before a gap and
forward repeats after it. Whether that is a T2T chromosome, an assembly artefact, or an
interstitial array is not decided by the code. The manifest carries `?` for these fixtures
until it is.

**Decision (B3a, B4):** resolved-code. Two interstitial rows, one per orientation. The
interstitial `teloType` column carries a junction class judged against the nearest
abutting block in either file, only across called sequence (a gap is a contig boundary):
`fusion` for q→p (TTAGGG array then CCCTAA array, reading left to right along the
sequence), `tail_to_tail` for p→q, `fragmentation` for p→p or q→q, `single` otherwise. The
blocked `?` expectations stay until B4.

### REG-007 — two definitions of block density, five lines apart

`docs/outputs.md:98` defines canonical block density as
`(fwdCan + revCan) * motifLen / (chromEnd - chromStart)`. The code uses OR-merged covered
bases over a gap-discounted denominator: `src/teloscope.cpp:403` and `:461` both compute
`density * (blockLen - getGapBases(...))`. The two differ whenever matches overlap or a
block spans a gap.

The interstitial reclassification guard at `:467` compounds this: it uses bare
`density * blockLen` with no gap discount, so a gapped array can be dense enough to keep as
an interstitial block and, in the same breath, not dense enough to be recognised as a
would-be telomere.

**Decision (B3a):** resolved-code. Blocks never cross an N gap, so every density
denominator is the plain block length: canonical coverage for terminal blocks, all-match
coverage for interstitial ones. `docs/outputs.md:98` is rewritten to that rule. Until B3a
the oracle follows the current gap-discounted code.

### REG-008 — threshold comparisons are inclusive, and nothing says so

`-y` rejects on strict `<` (`src/teloscope.cpp:403`), as do `-l` (`:401`) and
`--min-block-counts` (`:322`). A block exactly at a threshold is therefore kept. The
documentation says "minimum", which is consistent but never states the tie direction, and
the shipped `density_edge.fa` fixture sits exactly on the `-y` tie without saying what
should happen.

**Decision (B1):** resolved-doc. The docs state that a block exactly at a threshold is kept.

### REG-006 — the strand label is not symmetric under reverse complement

`computeStrandLabel` (`include/teloscope.h:236-242`) returns `p` iff
`fwd*1000 > counts*666` and `q` iff `fwd*1000 < counts*333`, so the balanced band is
`[0.333, 0.666]` while its mirror image is `[0.334, 0.667]`. A block with 1000 matches of
which 333 are forward is labelled `b`; its reverse complement has 667 forward and is
labelled `p`. A sequence and its own reverse complement get different answers.

The symmetric pair would be 666/334 or 667/333. The interstitial arm rule
`(start + len <= firstBase + lastBase) ? 'p' : 'q'` (`src/teloscope.cpp:482`) has the same
problem at the exact midpoint, where both a block and its mirror resolve to `p`.

Note that the neighbouring asymmetry is *not* a defect: the arm election uses `>` for p at
`:518` and `>=` for q at `:524`, and because blocks ascend by start that is exactly what
makes both arms elect their outermost block. It is correct as written.

**Decision (B5):** resolved-code. Exact thirds, symmetric by construction, tunable with one
fraction `--label-threshold` (default 0.667): `p` when the forward share exceeds it, `q`
when the forward share is below one minus it, `b` between. `closestEnd` at an exact
midpoint reads `p`, and the docs say so.

### REG-011 — reverse-oriented path components are silently skipped

`src/input.cpp:1013` has an empty `else` branch for path components in `-` orientation, so
such a segment contributes its length but none of its matches. Pre-existing, and there is
no fixture for it.

**Decision (B6):** resolved-code, with a GFA fixture carrying a reverse-oriented path
component.

### REG-010 — two summary headers are missing a space

`src/teloscope.cpp:1212` and `:1218` print `Chromosome Telomere Counts+++` and
`Chromosome Telomere/Gap Completeness+++` without the space every other banner has.

**Decision (B2):** resolved-code. Both banners gain the missing space; the checked-in
snapshots are regenerated to match.

### REG-012 — block filters are not pure post-filters

Rejecting a p-arm block on `-l` or `-y` leaves `pTo == 0`, which relaxes the q-arm guard
at `src/teloscope.cpp:419` and the interstitial exclusions at `:452-454`. Tightening a
filter can therefore *add* a q-arm block or an interstitial block.

The invariant harness therefore asserts no monotonicity in `-t`, `-l`, or `-y` at all.

**Decision (B3a):** resolved-code. A rejected p block keeps its extent for the q-arm guard.

### REG-009 — retired vocabulary inside the specification

`docs/classification.md:62` still says "which is what the scaffold type `discordant`
reports", after `discordant` stopped being a scaffold type and became an anomaly flag.

**Decision (B1):** resolved-doc.

### REG-015 — `-n/--manual-curation` and `-a/--out-fasta` do nothing

`-n` gated contig-terminal blocks into the BED from cc242d4 (2026-03) until d3207f8 removed
the gate with the rest of the per-contig scan; `docs/parameters.md:103` now calls it
"accepted for compatibility" while `docs/troubleshooting.md:160` and the `--help` text
(`src/main.cpp:582`) still describe the gate. `-a` sets `outFasta` and nothing reads it, at
v0.1.5 and now.

**Decision (B3c, B6):** resolved-code. Per-contig scanning returns: every contig end is
scanned in fast mode, a telomere at an internal contig end is a `teloType=contig` row, and
`-n` includes those rows. `-a` writes `<input>_telomeres.fa`, one record per reported
block, header `>chr:start-end teloLabel=… closestEnd=… teloType=… teloLen=…`. Between B3a
(no gap bridging) and B3c, an array at an internal contig end is an interstitial row.

*Fixtures:* three-contig scaffold with an internal contig-end telomere, run with and without `-n`

### REG-016 — the report script reads strand where it means arm

`scripts/teloscope_report.py` reads `["label"]` (strand composition) wherever it decides
which arm a block sits on and never reads `["arm"]`. Every arm decision in the PDF
report is therefore a strand decision, which is wrong on every discordant block.

**Decision (B4):** resolved-code. The script reads `closestEnd`; a failing case is added to
`scripts/test_teloscope_report.py` first.

### REG-017 — a short end-adjacent canonical array vanishes

`src/teloscope.cpp:465-467` drops an interstitial segment that sits in the terminal zone
and is canonical-dense enough to have been a telomere, whatever the reason the arm loops
rejected it: too short for `-l`, or refused by the "never the same physical array" rule at
`:425-426`. On `jn_short_q_then_p.fa` (two abutting 360 bp arrays of opposite orientation,
each within tolerance of its nearer end) the second array is refused as the p arm's own
array and then dropped by the guard, so it is reported nowhere. `docs/algorithm.md:19` says
everything that is not terminal is scored as interstitial, and `docs/outputs.md:104` makes
both arrays terminal.

**Decision (B3a):** resolved-code. The guard is deleted and the array becomes an
interstitial row, or, under REG-014, the outermost array at its end. `ORA-08-its-recall`
skips the too-short case until B3a.

*Fixtures:* `jn_short_q_then_p`, `jn_short_p_then_q` · *Test:* both harnesses

### REG-018 — sequence beyond `-t` inside an accepted seed is reported nowhere

On `boundary_extend_its.fa` at `-i -t 300`, each 600 bp terminal array is capped at 300 bp
by `-t`, but the accepted seed's extent (`pTo`, `qFrom` at `src/teloscope.cpp:411`, `:444`)
still spans all 600, and the interstitial scan starts beyond it (`:452-454`). The other
300 bp of each array are neither terminal nor interstitial. Found by `ORA-08-its-recall`.

**Decision (B3b):** resolved-code. The terminal chain stops at `-t` and the remainder is an
interstitial row; its junction class (`fragmentation`, since it abuts a same-orientation
block) arrives with B4.

*Fixtures:* `boundary_extend_its` · *Test:* `scripts/check_invariants.py`

## Resolved

*(none yet)*
