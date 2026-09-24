# Benchmark results

Reproduce with:

```bash
cargo build --release
./benchmark/fetch_data.sh nextseq  # the two NextSeq runs used below (~0.9 GB)
./benchmark/run_benchmark.sh       # picks up fastqc/falco if on PATH
```

## Setup

| | |
|---|---|
| Machine | Intel Core i9-9900K @ 3.6 GHz, 8 cores / 16 threads, 32 GB RAM |
| OS | macOS 15 (Darwin 24.6.0), x86-64 |
| RastQC | `main` at `faa547b` (unreleased) and 0.2.0 — release builds, `lto = true` |
| Falco | 1.2.5 (Bioconda) |
| FastQC | not run here — no JRE on this machine; see the FastQC table in the README |
| Method | 3 repetitions per row, **fastest** wall time, peak RSS across the same runs via `/usr/bin/time -l` |
| Measured | 2026-09-24 — every Falco, 0.2.0 and `main` row below in one sitting |
| Datasets | Provenance and measured read counts in [`DATA.md`](DATA.md) |

**Which statistic.** The harness reports the *fastest* of its repetitions, not
the median: timing interference on this machine is one-sided and occasionally
outlasts two of three repetitions, which a median cannot survive (one FastQC
measurement came out 47.7 s / 242 s / 47.7 s). Every Falco, 0.2.0 and `main`
row on this page was measured under that statistic, in one sitting, so the
rows are comparable to each other and to the paper's figures.

**The RastQC 0.1.0 rows are the exception**, marked †. They have not been
re-measured and are still the 0.2.0-era median of 3. A median can only be
equal to or slower than the minimum of the same runs, so those rows are if
anything slightly pessimistic, and the `vs Falco` ratios computed against them
slightly overstate the gap.

**On thread counts.** Falco is single-threaded: its `-t/--threads` flag is
documented in its own `--help` as *"NOT YET IMPLEMENTED IN FALCO"*. The
`rastqc -t 1` rows are therefore the like-for-like, one-core comparison; the
`-t 4` rows show what RastQC does when given more cores.

## Short-read (Illumina NextSeq 500)

### DRR045135_1 — 18.7M reads, 72 bp, 828 MB gzipped

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 30.55 s | 86 MB | — |
| **RastQC `main` `-t 1`** | **9.83 s** | 73 MB | **3.1× faster** |
| **RastQC `main` `-t 4`** | **6.66 s** | 70 MB | **4.6× faster** |
| RastQC 0.2.0 `-t 1` | 10.43 s | 118 MB | 2.9× faster |
| RastQC 0.2.0 `-t 4` | 7.10 s | 150 MB | 4.3× faster |
| RastQC 0.1.0 `-t 1` † | 64.9 s | 114 MB | 2.1× *slower* |
| RastQC 0.1.0 `-t 4` † | 20.9 s | 322 MB | 1.5× faster |

### DRR048760 — 1.2M reads, 67 bp, 54 MB gzipped

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 2.34 s | 88 MB | — |
| **RastQC `main` `-t 1`** | **0.78 s** | 74 MB | **3.0× faster** |
| **RastQC `main` `-t 4`** | **0.66 s** | 78 MB | **3.5× faster** |
| RastQC 0.2.0 `-t 1` | 0.80 s | 104 MB | 2.9× faster |
| RastQC 0.2.0 `-t 4` | 0.70 s | 126 MB | 3.3× faster |
| RastQC 0.1.0 `-t 1` † | 4.37 s | 107 MB | 1.9× *slower* |

### Both files in one invocation

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 32.80 s | 94 MB | — |
| **RastQC `main` `-t 4`** | **7.71 s** | 133 MB | **4.3× faster** |
| RastQC 0.2.0 `-t 4` | 7.60 s | 211 MB | 4.3× faster |

This is the one row where `main` does not beat 0.2.0 on wall time — 7.71 s
against 7.60 s, a 1.4% difference that is inside this machine's run-to-run
spread. Peak RSS drops from 211 MB to 133 MB over the same pair.

## What changed between 0.2.0 and `main`

Measured on the same machine and in the same sitting as the tables above.
Reported QC values are unchanged: `fastqc_data.txt` from a 0.2.0 `-t 1` run
matches `main` at `-t 4` and `-t 16` byte for byte, on both a short-read and a
long-read file. The changes themselves are in the
[CHANGELOG](../CHANGELOG.md#unreleased).

### ERR5897746_1 — 4.3M reads, 126 bp, 320 MB gzipped

| | 0.2.0 | `main` |
|---|---|---|
| wall, `-t 1` | 4.14 s | **3.74 s** |
| wall, `-t 4` | 2.53 s | **2.28 s** |
| wall, `-t 16` | 2.49 s | **2.31 s** |
| system time, `-t 4` | 0.27 s | **0.09 s** |
| peak RSS, `-t 1` | 126 MB | **72 MB** |
| peak RSS, `-t 4` | 149 MB | **73 MB** |
| peak RSS, `-t 16` | 147 MB | **75 MB** |

### DRR242198_1 — 75.8k ONT reads, 5.9 kb mean, 100.6 kb max, 406 MB

| | 0.2.0 | `main` |
|---|---|---|
| wall, `-t 1` | 2.53 s | **2.42 s** |
| wall, `-t 4` | 2.98 s | **2.39 s** |
| wall, `-t 16` | 2.74 s | **2.40 s** |
| peak RSS, `-t 1` | 460 MB | **422 MB** |
| peak RSS, `-t 4` | 1045 MB | **410 MB** |
| peak RSS, `-t 16` | 1021 MB | **415 MB** |

0.2.0 got *slower* on this long-read file as `-t` rose — 2.98 s at `-t 4`
against 2.53 s at `-t 1`, for 2.3× the memory — because its per-file worker
cap keyed off the filename extension, so an ONT run still got four workers and
each built its own set of per-length GC models. `main` samples mean read
length from the first block instead and gives such files one worker.

## What changed between 0.1.0 and 0.2.0

[Issue #12](https://github.com/Huang-lab/RastQC/issues/12) reported RastQC
running ~2× slower than Falco on NextSeq runs. That reproduced exactly: 0.1.0
at `-t 1` was 2.1× slower than Falco on DRR045135_1 and 1.9× slower on
DRR048760. Profiling found the causes were not where the tool's own
documentation assumed:

| Fix | Effect |
|-----|--------|
| Adapter search: brute-force scan of every start position for every adapter → one SIMD Aho-Corasick pass that rejects the (overwhelmingly common) no-adapter read outright | Was **53% of all CPU**, now 6% |
| Gzip: flate2's default `miniz_oxide` backend → the pure-Rust `zlib-rs` backend | Inflate was ~2.6× slower than zlib and runs on the serial reader thread, capping every gzipped input |
| FASTQ parsing moved off the single reader thread into the workers, on borrowed slices | Removed ~7 heap allocations and a UTF-8 validation per read |
| Per-base loops (`basic_stats`, `per_base_sequence_content`, `per_sequence_gc`) made branch-free | ~57% of the remaining per-read work, roughly halved |
| `-t` made a total budget rather than a per-file multiplier | `rastqc *.fastq.gz -t 16` on 6 files went from 5.1 GB and 7.0 s to 0.6 GB and 1.3 s |

## Memory and thread count

0.1.0 gave every worker thread a full copy of every module's state, so peak
memory grew with `-t` — and because file-level and within-file parallelism
both took the full `-t`, it grew with the *product*. 0.2.0 keeps the
whole-file modules on one instance and treats `-t` as a budget:

| Run | 0.1.0 † | 0.2.0 † |
|---|---|---|
| 1 file (828 MB gz), `-t 4` | 322 MB / 20.9 s | **163 MB / 7.4 s** |
| 1 file (828 MB gz), `-t 16` | 1119 MB / 22.3 s | **146 MB / 6.9 s** |
| 6 files (240 MB each), `-t 8` | 4052 MB / 6.3 s | **546 MB / 1.6 s** |
| 6 files (240 MB each), `-t 16` | 5148 MB / 7.0 s | **586 MB / 1.3 s** |

Note that 0.1.0 got *slower* as `-t` rose past 4 while its memory kept
climbing; 0.2.0 does not. This table is 0.2.0-era and has not been re-measured
under the current harness; for `main`'s memory against 0.2.0's, see the tables
above, where a `-t 4` short-read run drops from 150 MB to 70 MB.

## Why these datasets, and not the small sample

`benchmark/data/yeast_200k.fastq.gz` is 200k reads and finishes in under half
a second for every tool tested. That measures process startup, not throughput,
and it is why 0.1.0's published results looked healthy while real NextSeq runs
did not. Run `benchmark/fetch_data.sh` before drawing any conclusion about
performance.

Two things that file's name implies are both false: it is a subsample of
SARS-CoV-2 amplicon data rather than yeast, and it has never been committed to
this repository in any revision. [`DATA.md`](DATA.md) records what every
benchmark and validation dataset actually is, measured from the files rather
than taken from metadata, including three files in the 0.1.0 paper's set whose
names do not match their contents.
