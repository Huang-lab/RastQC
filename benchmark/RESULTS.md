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
| RastQC | 0.2.0 (release build, `lto = true`) |
| Falco | 1.2.5 (Bioconda) |
| FastQC | not run here — no JRE on this machine; see the FastQC table in the README |
| Method | 3 repetitions per row, median wall time, peak RSS via `/usr/bin/time -l` |
| Datasets | Provenance and measured read counts in [`DATA.md`](DATA.md) |

**The harness no longer reports the median.** Every wall time on this page was
measured with a median of 3, which is what `run_benchmark.sh` reported when
0.2.0 was released. It now reports the *fastest* of the repetitions instead:
timing interference on this machine is one-sided and occasionally outlasts two
of three repetitions, which a median cannot survive (one FastQC measurement
came out 47.7 s / 242 s / 47.7 s). A fresh run is therefore not strictly
comparable to the tables below — it can only come out equal or faster — and
the paper's figures, which are regenerated from the CSV, already use the
fastest. Re-measuring this page is tracked for the next release.

**On thread counts.** Falco is single-threaded: its `-t/--threads` flag is
documented in its own `--help` as *"NOT YET IMPLEMENTED IN FALCO"*. The
`rastqc -t 1` rows are therefore the like-for-like, one-core comparison; the
`-t 4` rows show what RastQC does when given more cores.

## Short-read (Illumina NextSeq 500)

### DRR045135_1 — 18.7M reads, 72 bp, 828 MB gzipped

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 32.6 s | 86 MB | — |
| **RastQC 0.2.0 `-t 1`** | **11.3 s** | 126 MB | **2.9× faster** |
| **RastQC 0.2.0 `-t 4`** | **6.9 s** | 158 MB | **4.7× faster** |
| RastQC 0.1.0 `-t 1` | 64.9 s | 114 MB | 2.0× *slower* |
| RastQC 0.1.0 `-t 4` | 20.9 s | 322 MB | 1.6× faster |

### DRR048760 — 1.2M reads, 67 bp, 54 MB gzipped

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 2.32 s | 88 MB | — |
| **RastQC 0.2.0 `-t 1`** | **0.83 s** | 95 MB | **2.8× faster** |
| **RastQC 0.2.0 `-t 4`** | **0.67 s** | 128 MB | **3.5× faster** |
| RastQC 0.1.0 `-t 1` | 4.37 s | 107 MB | 1.9× *slower* |

### Both files in one invocation

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 35.1 s | 94 MB | — |
| **RastQC 0.2.0 `-t 4`** | **8.6 s** | 208 MB | **4.1× faster** |

## What changed between 0.1.0 and 0.2.0

[Issue #12](https://github.com/Huang-lab/RastQC/issues/12) reported RastQC
running ~2× slower than Falco on NextSeq runs. That reproduced exactly: 0.1.0
at `-t 1` was 2.2× slower than Falco on DRR045135_1 and 1.9× slower on
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

| Run | 0.1.0 | 0.2.0 |
|---|---|---|
| 1 file (828 MB gz), `-t 4` | 322 MB / 20.9 s | **163 MB / 7.4 s** |
| 1 file (828 MB gz), `-t 16` | 1119 MB / 22.3 s | **146 MB / 6.9 s** |
| 6 files (240 MB each), `-t 8` | 4052 MB / 6.3 s | **546 MB / 1.6 s** |
| 6 files (240 MB each), `-t 16` | 5148 MB / 7.0 s | **586 MB / 1.3 s** |

Note that 0.1.0 got *slower* as `-t` rose past 4 while its memory kept
climbing; 0.2.0 does not.

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
