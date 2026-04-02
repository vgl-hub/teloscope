# teloscope validation

`validateFiles/` contains `.tst` manifests for `build/bin/teloscope-validate`.

## Running the suite

Run one file:

```sh
build/bin/teloscope-validate validateFiles/balanced.fa.0.tst
```

Run the full directory:

```sh
build/bin/teloscope-validate validateFiles
```

Use `-c` to print the exact Teloscope command for each test and `-v` to print output differences.

## Legacy `.tst` format

A legacy `.tst` file has:

1. one line of Teloscope arguments
2. either `embedded` or the path to an expected stdout file
3. embedded expected stdout lines when `embedded` is used

Example:

```text
-f testFiles/t2t.fa
embedded
...
```

## Directive-mode `.tst` format

Directive mode is used when the first non-empty line after the command starts with `expect_` or `gfa_`.

Supported directives:

- `expect_exit <code>`
- `expect_stdout ignore|<path>`
- `expect_stderr_substr <text>`
- `expect_output_name <basename>`
- `expect_gfa_header <version>`
- `gfa_expect <path-to-tsv>`
- `gfa_preserve_input strict|subset|segments_only|skip`

Example:

```text
-f testFiles/gfa_pathless_small.gfa -o %OUTDIR% -j 1
expect_exit 0
expect_stdout ignore
expect_output_name gfa_pathless_small.gfa.telo.annotated.gfa
expect_gfa_header 1.2
gfa_expect testFiles/expected/gfa/gfa_pathless_small.tsv
gfa_preserve_input strict
```

`%OUTDIR%` is replaced with a temporary output directory created for that test.

## GFA semantic checks

`gfa_expect` points to a TSV with these columns:

1. `path`
2. `segment`
3. `terminal_role`
4. `path_orient`
5. `node_name`
6. `seg_edge_orient`
7. `connector_type` (optional, defaults to `L`)
8. `connector_value` (optional, defaults to `0M`)
9. `tl_bp`

The validator parses the annotated GFA and checks:

- expected telomere segments exist
- expected telomere links exist
- `LN:i:6`, `RC:i:6000`, and `TL:i:<bp>` tags are correct
- no unexpected telomere segments or telomere links were added
- non-telomere graph content is preserved under the requested mode

Preservation modes:

- `strict`: non-telomere lines must match exactly as a multiset
- `subset`: every input non-telomere line must still exist in output
- `segments_only`: segment content must survive semantic normalization
- `skip`: do not compare the original graph payload

## Regenerating tests

`build/bin/teloscope-generate-tests` can regenerate legacy stdout-based expected outputs. Use it only when the current behavior is accepted.

Directive-mode GFA tests are checked in directly and should be edited by hand.
