[Back to README](index.md)

# Testing

Teloscope has four test layers:

- invariant and declared-intent checks over generated fixtures (`make test-synthetic`)
- `.tst` manifests in `validateFiles/` for the main binary
- focused shell checks in `scripts/`
- Python regression checks for the plotting code

## Invariants and declared intent

```sh
make test-synthetic     # fixture check, declared intent, and invariants
make test-intent        # declared intent alone
make test-invariants    # invariants alone
```

`TELOSCOPE=/path/to/binary` overrides the binary for either script.

`scripts/test_synthetic_intent.py` asserts `testFiles/synthetic/manifest.tsv` against what
the binary reports. `scripts/check_invariants.py`, with `scripts/teloscope_model.py`, needs
no recorded expected values: derivation checks re-derive a report field from the BED files,
oracle checks measure telomere recall directly against the input FASTA, and cross-run
checks compare two runs against each other — determinism, thread count, fast mode against
full scan, `-x` monotonicity at match level, a tip-window check (a small `-t` gives the same
terminal BED as the default), and a manual-curation check (`-n` only adds contig rows, and
only to the terminal BED).

A known deviation is waived in `validateFiles/intent_waivers.tsv` or
`validateFiles/invariant_waivers.tsv` against a [conflict register](conflicts.md) id, never
by editing the manifest to match the binary; a waiver that stops failing fails the run.

### Fixtures

```sh
bash testFiles/generate_synthetic.sh          # write fixtures and the manifest
bash testFiles/generate_synthetic.sh --check  # regeneration must be a no-op
bash testFiles/generate_synthetic.sh --list   # fixture ids
```

`make fixtures` and `make fixtures-check` wrap the same script. It reads its thresholds
from `include/input.h`, `include/teloscope.h`, and `src/teloscope.cpp`, and refuses to run
when one has moved.

Declaration format, as used by `fx` in `testFiles/synthetic_fixtures.sh`:

- `fx <id> <path> <record_spec> <flags> <expect> <intent>` declares one fixture; `xfx` declares a checked-in file the script does not write.
- `<expect>` is `key=value` pairs joined by `;`; for a multi-record file, one `<expect>` per record joined by `|`, or a single one for all.
- `type` t2t, incomplete or none; `anom` `.` or comma-joined anomaly flags; `telo` terminal row count; `labels` lowercase arm letters or none; `gaps` gap rows.
- `gran` one token per terminal row by start (uppercase an arm, lowercase a contig row, `*` after a discordant row); empty with no row.
- `its` interstitial block count in every full-scan run (`-i`, `-n`, `-r`, `-g`, `-e`, or `-m`), `-` otherwise.
- A key missing from `<expect>` reads `?`: blocked on a `docs/conflicts.md` entry named in `validateFiles/intent_blocked.tsv`.

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
- `expect_file`
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

`%OUTDIR%` is replaced by a per-test temporary directory. GFA expectations are semantic, not raw file diffs. See [validateFiles/README.md](https://github.com/vgl-hub/teloscope/blob/main/validateFiles/README.md) for the full format.

Use `expect_file <output-basename> <golden-path>` for generated BED, BEDgraph, or report files. The validator compares the file under `%OUTDIR%` with the golden after dropping blank and `#`-prefixed lines. The directive is repeatable, so one manifest can check several companion files.

## Regenerate legacy expected outputs

Only regenerate expected outputs when the current behavior is accepted:

```sh
make regenerate
build/bin/teloscope-generate-tests
```

Directive-mode manifests and their golden files are hand-authored and checked in directly.

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

## Assembly record filter regression script

```sh
make test-filters
```

This builds temporary FASTA and GFA fixtures and checks exact-ID and prefix selection, include/exclude precedence, selector validation, compressed input, line endings, unsupported modes, and output isolation. The CI workflow runs the same script on Linux, macOS, and Windows.

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
