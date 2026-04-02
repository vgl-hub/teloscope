# GFA Tests Apr 1

This file records the GFA coverage added on April 1 for Teloscope's path-aware annotation mode.

## Scope

The new suite adds tests only. It does not change existing fixtures, existing `.tst` manifests, or production GFA code.

The new tests live in:

- `testFiles/*.gfa`
- `testFiles/expected/gfa/*.tsv`
- `validateFiles/*.tst`

## Validator Contract

Legacy `.tst` files still work exactly as before.

New GFA-focused `.tst` files use directive mode:

- `expect_exit <code>`
- `expect_stdout ignore`
- `expect_output_name <basename>`
- `expect_gfa_header <version>`
- `gfa_expect <tsv>`
- `gfa_preserve_input <strict|subset|segments_only|skip>`

`gfa_expect` is a semantic oracle. It does not diff the whole annotated GFA byte-for-byte. Instead it checks:

- the exact synthetic telomere node set
- `LN:i:6`, `RC:i:6000`, and `TL:i:<len>` on each telomere node
- exactly one telomere link per expected node
- `L telomere_* + <segment> <orient> 0M`
- `RC:i:0` on telomere links
- no unexpected telomere nodes or telomere links

`gfa_preserve_input` controls how strictly the non-telomere graph is compared:

- `strict`: non-telomere input lines must match the non-telomere output lines exactly as a multiset
- `subset`: every non-telomere input line must still be present in the output, but duplicates in output are tolerated
- `segments_only`: non-telomere segments must be preserved semantically even if format changes, used for GFA2-to-GFA1 normalization
- `skip`: do not compare the original graph payload

## Fixture Catalog

### Pathless fixtures

- `gfa_pathless_small.gfa`
  - baseline pathless annotation
  - covers start-only, end-only, both-ends, and no-telomere segments

- `gfa_multiblock_same_end_small.gfa`
  - two terminal blocks on the same physical end
  - locks down current one-node dedup behavior for that end

- `gfa_noseq_small.gfa`
  - no-sequence `S ... * LN:i:...`
  - asserts that no telomere nodes are emitted

- `gfa_headerless_small.gfa`
  - headerless GFA1 input
  - asserts successful normalization to a `VN:Z:1.2` output header

- `gfa_internal_only_small.gfa`
  - repeat block away from the scanned tips
  - asserts no false terminal telomere node in tips-only mode

- `gfa_plant_small.gfa`
  - non-default canonical motif
  - verifies GFA mode respects `-c CCCTAAA`

### Path-aware fixtures

- `gfa_path_orient_pairs_small.gfa`
  - two multi-segment paths
  - covers all four start/end orientation rules:
    - first `+`
    - last `+`
    - first `-`
    - last `-`

- `gfa_single_seg_paths_small.gfa`
  - same terminal segment reused in `+` and `-` paths
  - verifies path orientation is encoded into node names and prevents over-dedup

### Format-normalization fixture

- `gfa2_small.gfa`
  - minimal GFA2 segment input
  - asserts semantic preservation of segments plus normalized `VN:Z:1.2` output

## Maintenance Notes

- Every new GFA case runs with `-j 1`.
  - This keeps the test surface deterministic and avoids coupling the regression suite to current thread-order sensitivity.

- The expectation files are intentionally small and semantic.
  - If node order changes, the tests still pass.
  - If telomere naming, orientation, or tag semantics change, the tests fail immediately.

- Existing FASTA and stdout-based validation stays untouched.
