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
| `discordant` | `gapped_discordant` | a terminal block is present but it sits closer to the opposite end of the sequence than the end it was scanned from |
| `none` | `gapped_none` | no scaffold-terminal telomere block was detected |

The `gapped_` prefix is added when the sequence contains assembly gaps.

## Decision order

Teloscope applies the rules in this order:

1. No scaffold-terminal blocks: `none`
2. Any longest terminal arm with the wrong positional orientation: `discordant`
3. One `p` block left of one `q` block: `t2t`
4. Both arms present but reversed: `misassembly`
5. One arm present, plus a second scaffold-terminal block of the same arm: `misassembly`
6. One terminal arm only: `incomplete`

Rule 5 only fires when the opposite arm is absent. When both arms are present, rule 3 or 4 already returns a result first, so a duplicated arm next to a normal opposite arm is not flagged as `misassembly`.

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

The `granular` column in `*_report.tsv` shows the arm pattern for each sequence. At most one p block and one q block are kept per sequence, so the column holds at most two letters.

Examples:

- `PQ`: a p arm and a q arm
- `P`: a p arm only
- `P*`: a p arm whose strand composition disagrees with its position

Uppercase marks the longest block for that arm. `*` marks arm and strand disagreement.
