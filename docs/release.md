[Back to README](index.md)

# Release checklist

Run through this once the release commit is on `main`.

## Before tagging

- Confirm `src/main.cpp` reports the new version.
- Keep `testFiles/` in the repo, since the manifests and CI read it. Expected outputs go in `testFiles/expected/`, never as root-level `testFiles/*_report.tsv` or `testFiles/*_gaps.bed`.
- Run the checks:

```sh
make all -j
bash .github/workflows/val.sh
make test-synthetic test-filters test-gaps
python3 scripts/test_teloscope_report.py
```

## GitHub release

```sh
git tag v0.1.6
git push origin v0.1.6
```

The `Create Release` workflow publishes:

- `teloscope.v0.1.6-linux.zip`
- `teloscope.v0.1.6-macOS.zip`
- `teloscope.v0.1.6-win.zip`
- `teloscope.v0.1.6-with_submodules.zip`, the source archive Bioconda builds from

## Bioconda

Bioconda does not follow this repository. Once the assets exist, open a PR on [`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/teloscope) that updates `recipes/teloscope/meta.yaml`:

- `version`: `0.1.6`
- `source.url`: `https://github.com/vgl-hub/teloscope/releases/download/v{{version}}/teloscope.v{{version}}-with_submodules.zip`
- `source.sha256`: the SHA-256 of the new `with_submodules` asset
- `doc_url`: pointing at `v{{ version }}`

The recipe installs only the binary. To ship the plotting scripts too, extend its `build.sh` and add the Python dependencies in the same PR.

## Zenodo

With the [Zenodo GitHub integration](https://zenodo.org/account/settings/github/) enabled for the repository, which an organization owner may need to approve, every GitHub release is archived with a DOI. A DOI needed before the release can be reserved from a Zenodo draft.
