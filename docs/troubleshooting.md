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

### Step size cannot be larger than window size

```
Error: Step size (2000) cannot be larger than window size (1000).
```

The step (`-s`) must be less than or equal to the window (`-w`). Note that if you set `-s` before `-w` on the command line, Teloscope validates both after all flags are parsed, so flag order does not matter:

```sh
# Both are valid
teloscope asm.fa -w 2000 -s 500
teloscope asm.fa -s 500 -w 2000
```

### Output directory is not writable

```
Error: Output directory '/path/to/dir' is not writable.
```

Teloscope checks that the output directory is writable before starting analysis. This catches permission errors and full filesystems early, before spending time on a genome scan that would fail at the output stage.

Fix: check that the directory exists, has write permissions, and that the filesystem has free space.

### Performance on large genomes

**Use ultra-fast mode (default).** Unless you need genome-wide window metrics (`-r`, `-g`, `-e`) or ITS detection (`-i`), the default mode only scans near contig ends and skips the interior. This is orders of magnitude faster on large assemblies.

**Thread count on shared systems.** By default, Teloscope uses all available cores. On shared HPC nodes, set `-j` to your allocation (e.g., `-j 16`) to avoid contention with other jobs.

**Window size on large genomes.** Small windows (`-w 100 -s 50`) on a 3 Gbp genome produce ~60 million windows per BEDgraph file. The default 1000 bp window is appropriate for most use cases.

**Output file sizes.** Each enabled BEDgraph flag (`-r`, `-g`, `-e`) produces one file per metric. On a 3 Gbp genome with 1 kbp windows, each BEDgraph is ~100 MB. The `-m` flag (all matches) can produce much larger files on repeat-rich genomes. Plan disk space accordingly.
