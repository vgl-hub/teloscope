[Back to README](index.md)

# Classification

FASTA mode classifies every sequence from its two arms. GFA mode has no classes.

## Arms

The p arm is the telomere anchored at the start of the first contig, and the q arm the one anchored at the end of the last contig; leading and trailing runs of `N` don't count. An arm must start inside the start zone, the smaller of `--terminal-tolerance` and `-t`. Its strand can be `p` or `q` at either end.

When one array runs from a contig's start to its end, both arms would claim it. It belongs to the end its strand points to: forward (`CCCTAA`) to the p end, reverse (`TTAGGG`) to the q end. A sequence that is `TTAGGG` from end to end is therefore `incomplete`, `q`, with no anomaly. The same holds for a record shorter than about twice `--terminal-tolerance` (6 kb by default), whose two start zones overlap. With `-n`, a scaffold end keeps such an array over the internal end of its contig.

## Classes

`type` counts the arms, and `anomaly` says whether they look like ordinary chromosome ends. An anomaly never changes the type.

| `type` | Arms |
| --- | --- |
| `t2t` | both |
| `incomplete` | one |
| `none` | neither |

| `anomaly` | Meaning |
| --- | --- |
| `.` | nothing to report |
| `discordant_p`, `discordant_q` | that arm's strand points the wrong way for its end, as at a head-to-head fusion or an inverted terminal repeat |
| `fragmented_p`, `fragmented_q` | that arm is built from more than one piece (`teloLen` below `end - start`): something breaks its `-k` chain, and the next piece starts within `-d` and outweighs the gap |

Flags are comma-separated, and are for a curator to judge, not errors.

Gaps are in neither column: `gaps` counts them, and the summary splits each type by them. Contig rows (`-n`) never count toward `type`, `telomeres`, `anomaly`, or the statistics; they appear in the terminal BED, not in the report.

## Labels

Every block has two labels, BED columns 5 and 6.

`teloLabel` is the strand. A terminal row is strand-pure, so it is `p` (forward) or `q` (reverse), never `b`; an array that alternates strands has no qualifying piece and becomes an interstitial `b` row. An interstitial row is labeled from its forward-strand share and `--label-threshold` (default `0.667`):

- `p`: above the threshold
- `q`: below one minus the threshold
- `b`: in between

`closestEnd` is the end a row belongs to: for a telomere, the end it is anchored at ([Arms](#arms)); for an interstitial row, the nearer end, `p` at the exact midpoint.

On an ordinary telomere the two agree, since a p arm carries the forward motif and sits at the start. They disagree on a discordant arm.
