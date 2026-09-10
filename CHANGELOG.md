# Changelog

All notable changes to RastQC are documented here.

## [0.2.0] — unreleased

Performance, memory and reproducibility release. Reported QC values are
unchanged from a 0.1.0 **sequential** run; see "Output" below for the two
places where 0.1.0's *parallel* runs disagreed with it.

### Performance

Prompted by [#12](https://github.com/Huang-lab/RastQC/issues/12), which
reported RastQC running ~2× slower than Falco on NextSeq runs. That
reproduced: on a public 18.7M-read NextSeq 500 run, 0.1.0 at `-t 1` took
64.9 s against Falco's 30.0 s. 0.2.0 takes 10.3 s at `-t 1` and 7.0 s at
`-t 4`. Full numbers and method in [`benchmark/RESULTS.md`](benchmark/RESULTS.md).

- **Adapter Content** searched every start position of every read for every
  adapter, byte by byte — 53% of all CPU time. Replaced with a single
  SIMD Aho-Corasick pass that rejects reads containing no adapter outright,
  falling back to per-adapter search only for reads that do contain one.
- **Gzip decoding** switched from flate2's default `miniz_oxide` backend to
  the pure-Rust `zlib-rs` backend, which decodes ~2.6× faster. Inflate runs
  on the serial reader thread, so it had been the ceiling on every gzipped
  input regardless of `-t`. Still no C toolchain or cmake required.
- **FASTQ parsing** moved off the reader thread into the worker threads,
  which now parse records as slices borrowed from a shared block. Removes
  roughly seven heap allocations and a UTF-8 validation per read.
- **Per-base loops** in Basic Statistics, Per base sequence content and Per
  sequence GC content rewritten branch-free (lookup table / histogram).
- **Kmer Content** no longer allocates a key for every kmer *occurrence*, and
  its per-position counters are `u32` rather than `u64`.
- Hot lookup tables now use an FxHash-style hasher instead of SipHash.

### Memory

- **`-t` is now a total thread budget, not a per-file multiplier.** Files are
  analyzed concurrently *and* each file runs its own worker pool; both levels
  were given the full `-t`, so `rastqc *.fastq.gz -t 16` spawned up to ~256
  workers. On 6 files that meant 4.6 GB resident — and it ran *slower* than
  `-t 4`. The budget is now split across the two levels.
- **Per-file worker count is capped** at where a single file stops benefiting
  (4 for compressed input, 8 for uncompressed). Past that the reader thread
  is the ceiling and extra workers only add module state to allocate and
  merge.
- Modules that need the whole read stream now run on **one** instance instead
  of one per worker, so their state no longer scales with `-t`.

Together: 1 file at `-t 16` went from 998 MB to 161 MB; 6 files at `-t 16`
from 4.6 GB to 714 MB.

### Output

Both of these made a QC number depend on how the run was parallelized:

- **Sequence Duplication Levels and Overrepresented sequences are now
  thread-independent.** Each worker kept its own 100k-sequence observation
  table; once a worker froze its table, occurrences of any sequence it had
  never admitted were lost, and each worker froze on a different subset. One
  real NextSeq run reported 75.2% deduplicated at `-t 1` (matching FastQC) but
  82.7% at `-t 4`. These modules now run on a single instance fed every read
  in file order, so results match a sequential run exactly at any `-t`.
- **Kmer Content and Per tile sequence quality are now thread-independent**
  for the same reason — both sample off a per-instance read counter, which
  only means "every Nth read of the file" if one instance sees the whole file.

### Fixed

- **Reports are now reproducible.** Kmer Content and Overrepresented
  sequences sorted by score alone, and the tied rows came out of a
  randomly-seeded `HashMap`, so two runs over the same file produced
  different output. Kmer Content is also truncated to 20 rows, so ties
  straddling the cutoff changed *which* kmers were reported at all. Both
  sorts now break ties on the sequence.
- Kmer Content no longer panics on a NaN p-value.
- `rastqc --version` read from a hardcoded string that could drift from the
  package version; it now reports the real one.
- **A file that fails to process now exits non-zero (3).** RastQC printed
  `Error processing <file>` and then exited 0, so an unreadable or malformed
  input passed a Nextflow/Snakemake gate silently — the opposite of what
  `--exit-code` is for.
- The declared MSRV was 1.70 while a dependency already required 1.75. It is
  now 1.75, and CI builds on exactly the declared version so the two can't
  drift again.

### Added

- `benchmark/fetch_data.sh` downloads a real public NextSeq run. The
  committed 200k-read sample finishes in under half a second for every tool,
  which measures startup rather than throughput.
- `benchmark/run_benchmark.sh` now benchmarks **Falco** alongside FastQC
  (both optional, auto-detected), reports peak memory, repeats each
  measurement and reports the median, and runs on Linux as well as macOS.
- [`benchmark/RESULTS.md`](benchmark/RESULTS.md) records the results and how
  to reproduce them.

## [0.1.0] — 2026-04-29

Initial release: 12 FastQC modules plus 3 long-read modules, FASTQ/SAM/BAM/
Fast5/POD5 input, HTML/text/ZIP/MultiQC-JSON output, and a built-in report
browser.
