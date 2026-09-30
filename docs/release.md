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
make test-synthetic test-filters test-gaps test-n50
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

Bioconda's autobump bot opens the version bump on [`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/teloscope) once the release is published, which the workflow does only after every asset is uploaded. Check that PR, or open one that updates `recipes/teloscope/meta.yaml`:

- `version`: `0.1.6`
- `source.url`: `https://github.com/vgl-hub/teloscope/releases/download/v{{version}}/teloscope.v{{version}}-with_submodules.zip`
- `source.sha256`: the SHA-256 of the new `with_submodules` asset
- `doc_url`: pointing at `v{{ version }}`

The recipe also installs `scripts/teloscope_report.py` in `bin/scripts/`, with Python, NumPy, pandas, and matplotlib as run dependencies; keep them when bumping.

## Galaxy

The [tools-iuc wrapper](https://github.com/galaxyproject/tools-iuc/tree/main/tools/teloscope) pins the version in `macros.xml`. Its tests check the report section titles and two terminal BED rows, so update them with the bump.

## Zenodo

Once the [Zenodo GitHub integration](https://zenodo.org/account/settings/github/) is on (an organization owner may need to approve it), every GitHub release gets a DOI. For a DOI before the release, reserve one from a Zenodo draft.
