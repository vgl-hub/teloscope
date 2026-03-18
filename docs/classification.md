[← Back to README](../README.md)

Chromosome classification
============

Teloscope classifies each chromosome based on its scaffold-terminal telomere blocks and whether the sequence contains assembly gaps (Ns). Only blocks within `-t`/`--terminal-limit` of a scaffold end are considered.

| Type | Gapped variant | Description |
|------|---------------|-------------|
| `t2t` | `gapped_t2t` | Both a p-arm and q-arm telomere, in the expected orientation (p left, q right) |
| `incomplete` | `gapped_incomplete` | Only one arm (p or q) has a scaffold-terminal telomere |
| `none` | `gapped_none` | No scaffold-terminal telomeres detected |
| `misassembly` | `gapped_misassembly` | Telomeres in wrong arrangement: p right of q, or the same arm appears twice at scaffold ends |
| `discordant` | `gapped_discordant` | A terminal block's strand composition contradicts its position (e.g., forward-dominant block near the 3' end). Not biological; cannot be resolved by simple reorientation |

**Classification logic.** For each chromosome, Teloscope finds the longest p-arm and q-arm blocks among scaffold-terminal blocks (by canonical coverage). It then applies these rules in order:

1. If no scaffold-terminal blocks exist → `none`
2. If the longest p or q block has **invalid orientation** (p closer to the 3' end, or q closer to the 5' end) → `discordant`
3. If both p and q exist and p is left of q → `t2t`
4. If both p and q exist but q is left of p → `misassembly`
5. If only one arm exists but a second block of the same arm is also scaffold-terminal → `misassembly`
6. If only one arm exists with no duplicates → `incomplete`

The `granular` column in the path summary shows the full block-by-block label string. Uppercase marks the longest block on each arm (e.g., `Pq` = one large p-arm block + one smaller q-arm block). An asterisk (`*`) after a label marks positional discordance: the block's strand identity does not match its chromosomal position.
