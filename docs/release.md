[Back to README](../README.md)

# Release checklist

Use this checklist after the release commit is on `main`.

## Repository contents

- Keep `testFiles/` on GitHub. The validator manifests in `validateFiles/` and the CI workflow call files from `testFiles/` directly.
- Do not add generated root-level `testFiles/*_report.tsv` or `testFiles/*_gaps.bed` outputs. Expected regression outputs live under `testFiles/expected/`.
- Confirm `src/main.cpp` reports the release version.
- Run the normal validation checks before tagging.

```sh
make all -j
bash .github/workflows/val.sh
python3 scripts/test_teloscope_report.py
bash scripts/test_gaps_bed.sh
```

## GitHub release

Push a tag from the release commit:

```sh
git tag v0.1.5
git push origin v0.1.5
```

The `Create Release` workflow should publish these assets:

- `teloscope.v0.1.5-linux.zip`
- `teloscope.v0.1.5-macOS.zip`
- `teloscope.v0.1.5-win.zip`
- `teloscope.v0.1.5-with_submodules.zip`

The `with_submodules` asset is the source archive used by the Bioconda recipe.

## Bioconda

Bioconda does not update from this repository automatically. After the GitHub release assets exist, open a PR against [`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/teloscope) that updates `recipes/teloscope/meta.yaml`:

- set `version` to `0.1.5`
- set `source.url` to `https://github.com/vgl-hub/teloscope/releases/download/v{{version}}/teloscope.v{{version}}-with_submodules.zip`
- replace `source.sha256` with the SHA-256 digest of the new `with_submodules` asset
- update `doc_url` to point at `v{{ version }}`

The current recipe installs only the `teloscope` binary. If Python plotting helpers should be shipped by Conda, update the Bioconda `build.sh` in the same PR and add any required runtime Python dependencies there.

## Zenodo

No repository file is required just to generate a Zenodo DOI URL. Zenodo must be connected to the public GitHub repository from the [Zenodo GitHub settings](https://zenodo.org/account/settings/github/), and for an organization repository an owner may need to approve the Zenodo GitHub app. Once enabled, Zenodo archives each new GitHub release and issues a DOI.

`CITATION.cff` or `.zenodo.json` can improve citation metadata, but they are not required for DOI minting. If a DOI is needed before the release is published, reserve it from a Zenodo draft record.
