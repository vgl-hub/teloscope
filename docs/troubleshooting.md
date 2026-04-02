[Back to README](../README.md)

# Troubleshooting

## Build fails because `gfalibs` is missing

Teloscope needs the `gfalibs` submodule for GFA support.

Fix:

```sh
git submodule update --init --recursive
make -j
```

## `--plot-report` fails because Python packages are missing

The report generator needs:

- `matplotlib`
- `numpy`
- `pandas`

Install them in the Python environment that will run `scripts/teloscope_report.py`.

## Compressed stdin is not supported

This fails:

```sh
cat asm.fa.gz | teloscope -o results/
```

Use one of these instead:

```sh
teloscope asm.fa.gz
zcat asm.fa.gz | teloscope -o results/
```

## GFA run writes no telomere nodes

Common reasons:

- the graph segment sequence is `*`, so there are no bases to scan
- the detected repeat block is shorter than `-l`
- the repeat density is below `-y`
- paths are present and the segment end is not path-terminal

If you want to inspect the graph result, make sure the output file exists:

```sh
ls results/*.telo.annotated.gfa
```

## Very large pattern sets are slow

Pattern count grows quickly when IUPAC expansion and edit distance are combined.

Examples:

| Input | `-x` | Concrete patterns after expansion |
| --- | --- | --- |
| `TTAGGG` | `1` | `38` |
| `TTAGGG` | `2` | `308` |
| `NNNGGG` | `0` | `127` |
| `NNNGGG` | `1` | `1180` |
| `NNNGGG` | `2` | `3367` |
| `NNNNNN` | `2` | `4096` |

If runtime is high, reduce ambiguity in `-p` or lower `-x`.

## Output directory is not writable

If you see:

```text
Error: Output directory '/path/to/dir' is not writable.
```

Check:

- the directory exists
- you have write permission
- the filesystem has free space
- another process is not writing conflicting outputs into the same directory
