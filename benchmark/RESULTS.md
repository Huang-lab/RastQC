# Benchmark results

Reproduce with:

```bash
cargo build --release
./benchmark/fetch_data.sh          # downloads the NextSeq run used below
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

**On thread counts.** Falco is single-threaded: its `-t/--threads` flag is
documented in its own `--help` as *"NOT YET IMPLEMENTED IN FALCO"*. The
`rastqc -t 1` rows are therefore the like-for-like, one-core comparison; the
`-t 4` rows show what RastQC does when given more cores.

## Short-read (Illumina NextSeq 500)

### DRR045135_1 — 18.7M reads, 144 bp, 828 MB gzipped

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 30.0 s | 86 MB | — |
| **RastQC 0.2.0 `-t 1`** | **10.3 s** | 115 MB | **2.9× faster** |
| **RastQC 0.2.0 `-t 4`** | **7.0 s** | 149 MB | **4.3× faster** |
| RastQC 0.1.0 `-t 1` | 64.9 s | 114 MB | 2.2× *slower* |
| RastQC 0.1.0 `-t 4` | 21.2 s | 330 MB | 1.4× faster |

### DRR048760 — 1.2M reads, 76 bp, 54 MB gzipped

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 2.30 s | 88 MB | — |
| **RastQC 0.2.0 `-t 1`** | **0.83 s** | 97 MB | **2.8× faster** |
| **RastQC 0.2.0 `-t 4`** | **0.68 s** | 121 MB | **3.4× faster** |
| RastQC 0.1.0 `-t 1` | 4.37 s | 107 MB | 1.9× *slower* |

### Both files in one invocation

| Tool | Wall | Peak RSS | vs Falco |
|------|------|----------|----------|
| Falco 1.2.5 | 32.5 s | 93 MB | — |
| **RastQC 0.2.0 `-t 4`** | **8.2 s** | 194 MB | **4.0× faster** |

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
| `-t` made a total budget rather than a per-file multiplier | `rastqc *.fastq.gz -t 16` on 6 files went from 4.6 GB and 4.2 s to 0.7 GB and 1.7 s |

## Memory and thread count

0.1.0 gave every worker thread a full copy of every module's state, so peak
memory grew with `-t` — and because file-level and within-file parallelism
both took the full `-t`, it grew with the *product*. 0.2.0 keeps the
whole-file modules on one instance and treats `-t` as a budget:

| Files × `-t` | 0.1.0 peak RSS | 0.2.0 peak RSS |
|---|---|---|
| 1 file, `-t 4` | 330 MB | 149 MB |
| 1 file, `-t 16` | 998 MB | 161 MB |
| 6 files, `-t 8` | 3476 MB | 319 MB |
| 6 files, `-t 16` | 4618 MB | 714 MB |

## A note on the committed sample data

`benchmark/data/yeast_200k.fastq.gz` is 200k reads and finishes in under half
a second for every tool tested — that measures process startup, not
throughput, and it is why 0.1.0's committed results looked healthy while real
NextSeq runs did not. Use `benchmark/fetch_data.sh` before drawing any
conclusion about performance.
