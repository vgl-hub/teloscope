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

Terminal blocks are labeled from the end they were scanned from, not from strand balance:

- `p`: found scanning from the start of the sequence
- `q`: found scanning from the end of the sequence

Interstitial blocks are labeled from strand balance instead:

- `p`: forward-strand dominant
- `q`: reverse-strand dominant
- `b`: balanced or mixed

The `granular` column in `*_report.tsv` shows the block pattern for each sequence. It includes every terminal block along the sequence, not just the two at the scaffold ends, so a gapped scaffold with contig-internal terminal blocks can produce a longer string than the examples below.

Examples:

- `PQ`: one dominant `p` block and one dominant `q` block
- `P`: one dominant `p` block only
- `Pq`: one large `p` block and one smaller `q` block
- `P*`: a `p` block in a discordant position

Uppercase marks the longest block for that arm. `*` marks positional discordance.
