[Back to README](index.md)

# Chromosome classification

This page applies to FASTA mode. GFA mode annotates graph segments but does not emit scaffold classes.

Teloscope classifies each sequence from scaffold-terminal telomere blocks only. A block counts as scaffold-terminal when it falls within `-t/--terminal-limit` of the sequence start or end.

## Output classes

| Type | Gapped variant | Meaning |
| --- | --- | --- |
| `t2t` | `gapped_t2t` | one `p` block at the left end and one `q` block at the right end |
| `incomplete` | `gapped_incomplete` | only one terminal arm is present |
| `misassembly` | `gapped_misassembly` | terminal arms exist but the arrangement is wrong, or the same arm appears twice |
| `discordant` | `gapped_discordant` | a terminal array points the wrong way for the end it sits at, which is an inverted terminal repeat or a fusion, not an ordinary chromosome end |
| `balanced` | `gapped_balanced` | a terminal arm carries both orientations in similar proportion, which is what a fusion or an inverted repeat looks like, not an ordinary chromosome end |
| `none` | `gapped_none` | no scaffold-terminal telomere block was detected |

The `gapped_` prefix is added when the sequence contains assembly gaps.

## Decision order

Teloscope applies the rules in this order:

1. No scaffold-terminal blocks: `none`
2. Either longest arm carries both strands in similar proportion: `balanced`
3. Either longest arm has the wrong positional orientation: `discordant`
4. Any further scaffold-terminal block beside the two longest arms: `misassembly`
5. Both arms present: `t2t`
6. Neither arm present: `none`
7. One terminal arm only: `incomplete`

Rule 4 runs before the `t2t` check, so a duplicated arm is reported as `misassembly` whether or not the opposite arm is present.

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

The `granular` column in `*_report.tsv` shows the arm pattern for each sequence. An arm may report more than one block, which is how a duplicated arm is detected, so the column holds one letter per scaffold-terminal block plus a `*` for each block whose strand disagrees with its arm.

Examples:

- `PQ`: a p arm and a q arm
- `P`: a p arm only
- `P*`: a p arm whose strand composition disagrees with its position

Uppercase marks the longest block for that arm. `*` marks arm and strand disagreement.
