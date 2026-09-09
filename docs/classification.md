[Back to README](index.md)

# Chromosome classification

This page applies to FASTA mode. GFA mode annotates graph segments but does not emit scaffold classes.

Teloscope classifies each sequence from scaffold-terminal telomere blocks only. A block counts as scaffold-terminal when it falls within `-t/--terminal-limit` of the sequence start or end.

## Output classes

How many telomeres a sequence has, and whether they are plausible, are two separate
questions, so they live in two separate columns. An anomaly never overwrites the count.

The `type` column answers completeness only:

| Type | Meaning |
| --- | --- |
| `t2t` | one arm at each end |
| `incomplete` | one terminal arm |
| `none` | no scaffold-terminal telomere block was detected |

The `anomaly` column answers plausibility. It reads `.` when there is nothing to report,
or one or more of the following, comma separated:

| Anomaly | Meaning |
| --- | --- |
| `discordant_p`, `discordant_q` | that arm points the wrong way for the end it sits at, which is an inverted terminal repeat or a fusion, not an ordinary chromosome end |
| `balanced_p`, `balanced_q` | that arm carries both orientations in similar proportion, which is what a fusion or a collapsed repeat looks like |
| `misassembly` | a further terminal block sits beside the two longest arms, so one end carries more than one array |

Neither anomaly is expected at an ordinary chromosome end. Both do occur in real biology,
`discordant` at genuine head-to-head fusions and in inverted terminal repeats, `balanced`
in ALT telomeres, hairpin telomeres and subtelomeric variant zones, so they are flagged
for a curator to judge rather than asserted as errors.

Gappedness is not part of either column. The `gaps` column already carries it, and the
assembly summary splits each completeness class on it.

## Decision order

The two columns are computed independently. Completeness is three cases: both longest
arms present is `t2t`, exactly one is `incomplete`, neither is `none`. The anomaly set is
then accumulated over the same blocks: for each longest arm, both orientations mixed sets
`balanced_`, otherwise a strand disagreeing with the arm sets `discordant_`; and any
terminal block that is not one of the two longest arms sets `misassembly`.

## Block labels

Every block carries two labels, in two separate columns of both BED files.

Column 5, `label`, is the strand composition, and it means the same thing in both files:

- `p`: forward-strand dominant, more than 66.6% of matches
- `q`: reverse-strand dominant, less than 33.3%
- `b`: balanced or mixed

Column 12, `arm`, is the position:

- `p`: the block sits nearer the start of the called sequence
- `q`: the block sits nearer the end

For an ordinary telomere the two agree, since a p arm carries the forward motif and a q arm the reverse. They disagree on an inverted terminal repeat, which is what the scaffold type `discordant` reports.

The `granular` column in `*_report.tsv` shows the arm pattern for each sequence. An arm may report more than one block, which is how a duplicated arm is detected, so the column holds one letter per scaffold-terminal block. A block whose orientation is mixed is marked `~`, and one whose strand disagrees with its arm is marked `*`. The two are different findings and are kept apart.

Examples:

- `PQ`: a p arm and a q arm
- `P`: a p arm only
- `P*`: a p arm whose strand composition disagrees with its position

Uppercase marks the longest block for that arm. `*` marks arm and strand disagreement.
