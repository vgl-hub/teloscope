[← Back to README](../README.md)

Troubleshooting
============

### Compressed input on stdin is not supported

```
Error: Compressed input on stdin is not supported. Decompress first:
  zcat file.fa.gz | teloscope -o results/
```

Teloscope can read `.fa.gz` files directly when passed as a file argument (`teloscope asm.fa.gz`). However, piping gzipped data through stdin is not supported. Decompress first with `zcat` or `gunzip -c`:

```sh
# This works
teloscope asm.fa.gz

# This does NOT work
cat asm.fa.gz | teloscope -o results/

# This works (decompress before piping)
zcat asm.fa.gz | teloscope -o results/
```

### Unusually high pattern count

```
Warning: 2400 patterns is unusually high and may be slow on large genomes.
  Consider fewer IUPAC wildcards or a lower -x value.
```

Pattern count grows multiplicatively: IUPAC wildcards expand first (e.g., `N` = 4 bases), then edit distance generates substitution variants on top. Some combinations:

| Input | `-x` | Patterns (after dedup) |
|-------|------|----------------------|
| `TTAGGG` | 1 | 38 |
| `TTAGGG` | 2 | 308 |
| `NNNGGG` | 0 | 127 |
| `NNNGGG` | 1 | 1,180 |
| `NNNGGG` | 2 | 3,367 |
| `NNNNNN` | 2 | 4,096 |

The trie handles large pattern sets correctly, but scanning a multi-gigabase genome against thousands of patterns generates proportionally more matches, increasing memory use and runtime. If you are using IUPAC wildcards to approximate fuzzy matching, use `-x 1` with concrete patterns instead.

### Output directory is not writable

```
Error: Output directory '/path/to/dir' is not writable.
```

Teloscope checks that the output directory is writable before starting analysis. Check that the directory exists, has write permissions, and that the filesystem has free space.
