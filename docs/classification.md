[Back to README](index.md)

# Chromosome classification

This page applies to FASTA mode. GFA mode annotates graph segments but does not emit scaffold classes.

Teloscope classifies each sequence from its two arms: the p arm is the telomere chain anchored at the first contig's start, and the q arm is the chain anchored at the last contig's end. These are the scaffold's own ends, so leading and trailing runs of `N` don't count. Each arm has to anchor within the start zone — the smaller of `--terminal-tolerance` and `-t` — of its contig end. An arm's own strand can be `p` or `q` regardless of which end it sits at.

When one array runs from a contig's start to its own end, both chains would claim it. It belongs to the end its strand points to: forward (the `CCCTAA` family) claims the p end, reverse (`TTAGGG`) the q end, and the other end has no telomere. A scaffold that is `TTAGGG` from end to end is therefore `incomplete`, `q`, with no anomaly.

The same rule covers a short record: one no longer than about twice the distance to end (`--terminal-tolerance`; `6 kb` at the defaults) has its two start zones overlap, so both ends reach its only array. The strand still decides which end it belongs to, and no `discordant` flag is raised. On chromosome-length records the two zones never overlap, and position decides as usual.

## Output classes

How many telomeres a sequence has, and whether they are plausible, are two separate
questions, so they live in two separate columns. An anomaly never overwrites the count.

The `type` column answers completeness only:

| Type | Meaning |
| --- | --- |
| `t2t` | both arms present |
| `incomplete` | one arm present |
| `none` | neither arm present |

The `anomaly` column answers plausibility. It reads `.` when there is nothing to report,
or one or more of the following, comma separated:

| Anomaly | Meaning |
| --- | --- |
| `discordant_p`, `discordant_q` | that arm's strand points the wrong way for the end it sits at, which is an inverted terminal repeat or a fusion, not an ordinary chromosome end |
| `fragmented_p`, `fragmented_q` | that arm is built from more than one piece (`teloLen` is less than `end - start`) |

Neither anomaly is expected at an ordinary chromosome end. Both do occur in real biology:
`discordant` at genuine head-to-head fusions and inverted terminal repeats, `fragmented`
where a short non-telomeric stretch splits an otherwise continuous array. They are flagged
for a curator to judge rather than asserted as errors. A terminal row is always
strand-pure; an array that alternates strands has no qualifying strand-pure piece and
becomes an interstitial `b` row instead.

Gappedness is not part of either column. The `gaps` column already carries it, and the
assembly summary splits each completeness class on it.

## Decision order

The two columns are computed independently from the two arms. Completeness is three
cases: both arms present is `t2t`, exactly one is `incomplete`, neither is `none`. The
anomaly set is then read off each arm: a strand that disagrees with the end it sits at
sets `discordant_p`/`discordant_q`; more than one piece in the chain sets
`fragmented_p`/`fragmented_q`.

Contig-terminal rows, written only with `-n`, never count toward `type`, `telomeres`, the
anomalies, or the length statistics. They appear lowercase in `granular` and nowhere else.

## Block labels

Every row carries two labels, in two separate columns of both BED files.

Column 5, `teloLabel`, is the strand.

A terminal row is strand-pure by construction, so its `teloLabel` is simply the strand of its repeats: `p` (forward) or `q` (reverse); it never reads `b`.

An interstitial row's matches can mix both strands, so its `teloLabel` comes from `--label-threshold` (default `0.667`) instead:

- `p`: forward-strand share above the threshold
- `q`: forward-strand share below one minus the threshold
- `b`: between the two

Column 6, `closestEnd`, is the end a row belongs to:

- for a telomere, the end it is the arm of: `p` at the start of the record, `q` at the end; when one array reaches both ends of a contig it is the end its strand points to
- for an interstitial row, the nearer record end; an exact midpoint reads `p`

For an ordinary telomere the two agree, since a p arm carries the forward motif and sits near the start. They disagree on a discordant arm, which the `anomaly` column flags as `discordant_p` or `discordant_q`.

The `granular` column in `*_report.tsv` shows one character per terminal row, in position order along the sequence. An arm is uppercase (`P` or `Q`). A contig-terminal row, written only with `-n`, is lowercase (`p` or `q`). A discordant row is followed by `*`.

Examples:

- `PQ`: a p arm and a q arm
- `P`: a p arm only, no q arm
- `P*`: a p arm whose strand disagrees with the end it sits at
- `Pq`: a p arm, plus a contig-terminal telomere elsewhere on the sequence (`-n`)
