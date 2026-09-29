[Back to README](index.md)

# Testing

`make all -j` builds the binary and its helpers: `teloscope-validate`, `teloscope-generate-tests`, and `teloscope-simulate`.

## Suites

| Command | Checks |
| --- | --- |
| `make test-synthetic` | fixtures regenerate unchanged, declared intent, and invariants |
| `build/bin/teloscope-validate validateFiles` | the `.tst` manifests |
| `bash .github/workflows/val.sh` | the `.tst` manifests, output file counts, and classification spot checks |
| `make test-filters` | record filters on FASTA and GFA |
| `make test-gaps` | gap BEDs against `testFiles/expected/` |
| `make test-n50` | scaffold and contig N50 |
| `make test-bam` | kept BAM records and malformed BAM |
| `make test-read-tl` | read telomere rows and report |
| `make test-bam-hardening` | a samtools oracle, zlib fault injection, and 512 BAM mutations |
| `make test-bam-coverage` | every executable line of the BAM code |
| `make test-bam-sanitize` | the BAM tests under AddressSanitizer and UndefinedBehaviorSanitizer |
| `python3 scripts/test_teloscope_report.py` | the report plots, on synthetic data |

`TELOSCOPE=/path/to/binary` points the Python checks at another binary.

CI runs the manifests and the synthetic, filter, BAM, and read checks on Linux, macOS, and Windows (`val.sh` outside Windows), the report check on two matplotlib versions, and the three BAM hardening targets on Linux. A weekly job repeats hardening and sanitizing with 4096 mutations.

## Synthetic fixtures

`scripts/test_synthetic_intent.py` checks the binary against `testFiles/synthetic/manifest.tsv`. `scripts/check_invariants.py` needs no expected values: it re-derives report fields from the BED files and compares runs with each other: repeated runs, thread counts, fast against full scan, a higher `-x` never losing matches, a small `-t` giving the same terminal BED, and `-n` adding only contig rows.

```sh
bash testFiles/generate_synthetic.sh          # write the fixtures and the manifest
bash testFiles/generate_synthetic.sh --check  # regeneration must change nothing
bash testFiles/generate_synthetic.sh --list   # fixture ids
```

`make fixtures` and `make fixtures-check` wrap the first two. The script reads its thresholds from `include/input.h`, `include/teloscope.h`, and `src/teloscope.cpp`, and refuses to run when one has moved.

Fixtures are declared in `testFiles/synthetic_fixtures.sh` and `testFiles/synthetic_axis_*.sh`:

- `fx <id> <path> <record_spec> <flags> <expect> <intent>` declares one fixture; `xfx` declares a checked-in file the script does not write.
- `<expect>` holds `key=value` pairs joined by `;`, one set per record joined by `|`, or one set for all records.
- Keys: `type`; `anom`, `.` or comma-joined flags; `telo`, the arm count; `labels`, the arm ends or `none`; `gaps`; `its`, the interstitial count in any full scan; `telolen`, an arm's length.

## `.tst` manifests

`teloscope-validate` runs the manifests in `validateFiles/`. A legacy manifest compares stdout with embedded text. A directive manifest checks exit codes, named output files, and GFA graphs semantically. `-c` prints each test's command:

```sh
build/bin/teloscope-validate -c validateFiles/gfa_pathless_small.tst
```

`make regenerate` builds `teloscope-generate-tests`, which rewrites the legacy manifests from the current binary; run it only when the new behavior is accepted. Directive manifests and their goldens are edited by hand. [validateFiles/README.md](https://github.com/vgl-hub/teloscope/blob/main/validateFiles/README.md) has the full format.

## BAM hardening

`make test-bam-hardening` builds its BAM from SAM with samtools, validates the output with `samtools quickcheck`, compares headers and records independently, and forces zlib failures through GNU linker wrappers. `make test-bam-coverage` reads gcov JSON and requires every executable line of `src/bam.cpp`, `src/bgzf.cpp`, and `src/read-filter.cpp`, with no exclusions; `COVERAGE_BRANCH_DETAILS=1` also prints branch counts. These targets need Linux.

## Repo layout

- `src/`, `include/`: C++ sources and headers
- `scripts/`: the report and plotting scripts, and the Python checks
- `testFiles/`: FASTA and GFA fixtures, with expected outputs in `testFiles/expected/`
- `validateFiles/`: `.tst` manifests
- `gfalibs/`: the GFA and I/O submodule
- `docs/`: this documentation
