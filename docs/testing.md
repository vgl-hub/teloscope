[← Back to README](../README.md)

Testing
============

Tests live under `testFiles/` and `validateFiles/`. Build and run with:
```sh
make validate
build/bin/teloscope-validate validateFiles
```

To regenerate expected test outputs from a known-good build:
```sh
make regenerate
build/bin/teloscope-generate-tests
```
