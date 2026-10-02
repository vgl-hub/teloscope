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

To test another binary, run a Python check directly with `TELOSCOPE=/path/to/binary`; the make targets set their own.

| Workflow | When | Runs |
| --- | --- | --- |
| `validate.yml` | every push, and pull requests to `main` | the manifests, `test-synthetic`, and the filter, BAM, and read checks on Linux, macOS, and Windows (`val.bat` on Windows); the three BAM hardening targets on Linux |
| `report.yml` | when the report scripts change | the report check on matplotlib 3.8.4 and 3.10.8, on Linux |
| `bam_hardening.yml` | weekly | hardening and sanitizing with 4096 mutations |

CI never runs `test-gaps` or `test-n50`; run them before a release.

## Synthetic fixtures

`scripts/test_synthetic_intent.py` checks the binary against `testFiles/synthetic/manifest.tsv`. `scripts/check_invariants.py` needs no expected values. It re-derives report fields from the BED files, checks the builder's guarantees (span at least `-l`, exact-repeat coverage at least `-y`, `teloLen` within the span, edges on repeat edges), and compares runs: repeats, thread counts, fast against full scan, a higher `-x` never losing matches, a small `-t` giving the same terminal BED, the reverse complement giving the mirror image, and `-n` adding only contig rows.

```sh
bash testFiles/generate_synthetic.sh          # write the fixtures and the manifest
bash testFiles/generate_synthetic.sh --check  # regeneration must change nothing
bash testFiles/generate_synthetic.sh --list   # fixture ids
```

`make fixtures` and `make fixtures-check` wrap the first two. The script reads its thresholds from `include/input.h` and `src/teloscope.cpp`, and refuses to run when one has moved.

Fixtures are declared in `testFiles/synthetic_fixtures.sh` and `testFiles/synthetic_axis_*.sh`:

- `fx <id> <path> <record_spec> <flags> <expect> <intent>` declares one fixture; `xfx <id> <path> <scaffold> <flags> <expect> <intent>` declares a checked-in file the script does not write.
- `<expect>` holds `key=value` pairs joined by `;`, one set per record joined by `|`, or one set for all records.

| Key | Expected value | When absent |
| --- | --- | --- |
| `type` | `t2t`, `incomplete`, or `none` | fails as unstated |
| `anom` | `.` or comma-joined flags | fails as unstated |
| `telo` | the arm count | fails as unstated |
| `labels` | the arm ends, or `none` | fails as unstated |
| `gaps` | the gap count | fails as unstated |
| `its` | the interstitial count, in any full scan | not checked |
| `telolen` | the `teloLen` of one arm | not checked |

## `.tst` manifests

`teloscope-validate` runs the manifests in `validateFiles/`. A legacy manifest compares stdout with embedded text. A directive manifest checks exit codes, named output files, and GFA graphs semantically. `-c` prints each test's command:

```sh
build/bin/teloscope-validate -c validateFiles/gfa_pathless_small.tst
```

`make regenerate` builds `teloscope-generate-tests`. Run from the repo root, it deletes and rewrites the legacy manifests from the current binary, adds the two `chr_only` ones, and leaves outputs in `testFiles/`, so review the diff. Directive manifests and their goldens are edited by hand; a new BED golden needs `git add -f`, since `.gitignore` matches it. [validateFiles/README.md](https://github.com/vgl-hub/teloscope/blob/main/validateFiles/README.md) has the full format.

## BAM hardening

`make test-bam-hardening` builds its BAM from SAM with samtools, validates the output with `samtools quickcheck`, compares headers and records independently, and forces zlib failures through GNU linker wrappers. `make test-bam-coverage` reads gcov JSON and requires every executable line of `src/bam.cpp`, `src/bgzf.cpp`, and `src/read-filter.cpp`, with no exclusions; `COVERAGE_BRANCH_DETAILS=1` also prints branch counts. These targets need Linux.

## Repo layout

- `src/`, `include/`: C++ sources and headers
- `scripts/`: the report and plotting scripts, and the Python checks
- `testFiles/`: FASTA and GFA fixtures, with expected outputs in `testFiles/expected/`
- `validateFiles/`: `.tst` manifests
- `tests/`: the BGZF fault-injection test
- `gfalibs/`: the GFA and I/O submodule
- `docs/`: this documentation
