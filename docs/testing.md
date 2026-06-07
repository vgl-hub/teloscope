[Back to README](../README.md)

# Testing

Teloscope has three main test layers:

- `.tst` manifests in `validateFiles/` for the main binary
- focused shell checks in `scripts/`
- Python regression checks for the plotting code

## Build the helper binaries

```sh
make all -j
```

If you only need the validator:

```sh
make validate
```

## Main validation suite

Run the checked-in `.tst` suite:

```sh
build/bin/teloscope-validate validateFiles
```

This is the same validator used by the CI workflow in `.github/workflows/validate.yml`.

## `.tst` formats

`teloscope-validate` supports two styles.

Legacy mode compares stdout against embedded text or an expected file.

Directive mode is used for GFA and file-oriented checks. Supported directives are:

- `expect_exit`
- `expect_stdout`
- `expect_stderr_substr`
- `expect_output_name`
- `expect_gfa_header`
- `gfa_expect`
- `gfa_preserve_input`

Minimal directive-mode example:

```text
-f testFiles/gfa_pathless_small.gfa -o %OUTDIR% -j 1
expect_exit 0
expect_stdout ignore
expect_output_name gfa_pathless_small.gfa.telo.annotated.gfa
expect_gfa_header 1.2
gfa_expect testFiles/expected/gfa/gfa_pathless_small.tsv
gfa_preserve_input strict
```

`%OUTDIR%` is replaced by a per-test temporary directory. GFA expectations are semantic, not raw file diffs. See [validateFiles/README.md](../validateFiles/README.md) for the full format.

## Regenerate legacy expected outputs

Only regenerate expected outputs when the current behavior is accepted:

```sh
make regenerate
build/bin/teloscope-generate-tests
```

Directive-mode GFA manifests and their TSV expectations are hand-authored and checked in directly.

## Report regression script

The plotting regression script now lives in `scripts/`:

```sh
python3 scripts/test_teloscope_report.py
```

It exercises `scripts/teloscope_report.py` directly with synthetic in-memory data.

## Gap BED regression script

```sh
bash scripts/test_gaps_bed.sh
```

This script compares generated `*_gaps.bed` files against the checked-in expected files in `testFiles/expected/`.

## BAM subset regression script

```sh
python3 scripts/test_bam_subset.py
```

The script uses only the Python standard library. It generates BAM/BGZF fixtures in a temporary directory and checks record preservation, scoring parity, malformed input handling, batching, thread determinism, and deterministic mutations.

## BAM hardening

The Linux hardening target adds an independent `samtools` oracle and zlib fault injection:

```sh
make test-bam-hardening
```

`samtools` generates the input BAM from SAM, validates the output with `quickcheck`, and checks headers and records independently. GNU linker wrappers force zlib initialization, compression, finalization, size, and output failures without production test hooks.

The strict coverage target uses GNU gcov:

```sh
make test-bam-coverage
```

It requires 100% of executable lines in `src/bam.cpp`, `src/bgzf.cpp`, and `src/read-filter.cpp`. The checker reads gcov JSON and counts a line as executable when gcov reports execution or an unexecuted block. Compiler-only cleanup braces are not counted as source executable lines. There are no source exclusions or manual coverage allowlists. Set `COVERAGE_BRANCH_DETAILS=1` to show raw gcov branch counts for diagnostics; they are not a gate because exceptions and standard-library templates add implementation-dependent branches.

The sanitizer target builds a temporary binary with AddressSanitizer and UndefinedBehaviorSanitizer:

```sh
make test-bam-sanitize
```

Pull requests run the standard integration suite, 512 deterministic mutations, strict coverage, sanitizers, zlib fault injection, and the `samtools` oracle in independent jobs. A weekly workflow repeats normal and sanitized mutation testing with 4096 cases.

## Useful local runs

Quick validator pass:

```sh
bash .github/workflows/val.sh
```

Single `.tst` file:

```sh
build/bin/teloscope-validate -c validateFiles/gfa_pathless_small.tst
```
