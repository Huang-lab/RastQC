# Changelog

All notable changes to RastQC are documented here.

## [Unreleased]

### Performance

Measured on a public 4.3M-read HiSeq run (ERR5897746_1, 320 MB gzipped),
median of 6 runs, against 0.2.0 on the same machine. Output is byte-identical
to 0.2.0 at every thread count — verified by diffing `fastqc_data.txt` between
a 0.2.0 `-t 1` run and a `-t 4` run of this build.

| | 0.2.0 | this build | |
|---|---|---|---|
| `-t 1` | 4.39 s | **3.96 s** | 1.11x |
| `-t 4` | 2.90 s | **2.54 s** | 1.14x |
| system time, `-t 4` | 0.43 s | **0.12 s** | |

- **The block reader now recycles its buffers.** It allocated a fresh 1 MB
  `Vec` for every block and dropped it once both consumers were done, so each
  megabyte of input cost an `mmap` and a `madvise` teardown — which profiled as
  the single largest cost in a run, above every QC module. Finished blocks are
  now returned to the reader and refilled. This is where the 3.5x drop in
  system time comes from.
- **The block buffer no longer reallocates on every block.** `next_block`
  fills in 128 KB chunks and stops once it has passed `BLOCK_SIZE`, so the last
  chunk pushed the length just past a capacity of exactly `BLOCK_SIZE` — a
  realloc and a 1 MB copy per block. Buffers are now allocated with room for
  that overshoot.
- **Per sequence quality scores** walked every quality string twice, once for
  the run's minimum character and once for the read's sum. The two are now one
  pass, and its per-read counter table is a flat array indexed by the mean
  quality character (which cannot leave `0..=255`) rather than a hash map.
- **Sequence Length Distribution** used a SipHash map for its once-per-read
  probe; it now uses the same FxHash the other hot tables use. It stays a map
  rather than an array because read length is unbounded — ONT reads reach
  megabases, and an array indexed by length would reintroduce the memory
  blowup 0.2.0 fixed.

### Added

- `benchmark/check_concordance.sh` — runs FastQC as the reference and compares
  RastQC's per-module PASS/WARN/FAIL calls against it, and Falco's too when
  Falco is installed. The paper's central correctness claim had no script in
  the repository, so the one number that most matters for a drop-in
  replacement could not be reproduced. It exits non-zero on any disagreement,
  so it also works as a gate.

### Changed

- `benchmark/fetch_data.sh` now fetches every dataset the benchmarks and the
  paper use, verifies each against the byte size ENA reports, and downloads
  large files as concurrent byte ranges. ENA throttles a single connection to
  roughly 200 KB/s, which is over five hours for the full set; twelve ranges
  measured ~2 MB/s.
- `benchmark/run_benchmark.sh` classifies long-read inputs by their measured
  mean read length instead of by filename. The previous `*_ont_*`/`*_pacbio_*`
  convention silently benchmarked real ENA files as short-read, and it now
  emits stable tool ids that `paper/analyze_benchmarks.py` matches on — the two
  had drifted apart, so the figures were being generated from no data.
- Both benchmark scripts now run on bash 3.2, which is what macOS ships.

## [0.2.0] — 2026-09-10

Performance, memory and reproducibility release. Reported QC values are
unchanged from a 0.1.0 **sequential** run; see "Output" below for the two
places where 0.1.0's *parallel* runs disagreed with it.

### Performance

Prompted by [#12](https://github.com/Huang-lab/RastQC/issues/12), which
reported RastQC running ~2× slower than Falco on NextSeq runs. That
reproduced: on a public 18.7M-read NextSeq 500 run, 0.1.0 at `-t 1` took
64.9 s against Falco's 32.6 s. 0.2.0 takes 11.3 s at `-t 1` and 6.9 s at
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
  workers. On 6 files that meant 5.1 GB resident — and it ran *slower* than
  `-t 4`. The budget is now split across the two levels.
- **Per-file worker count is capped** at where a single file stops benefiting
  (4 for compressed input, 8 for uncompressed). Past that the reader thread
  is the ceiling and extra workers only add module state to allocate and
  merge.
- Modules that need the whole read stream now run on **one** instance instead
  of one per worker, so their state no longer scales with `-t`.

Together: one 828 MB file at `-t 16` went from 1119 MB to 146 MB, and six
240 MB files at `-t 16` from 5148 MB to 586 MB — while getting 3× and 5×
faster respectively.

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
- The declared MSRV was 1.70 while dependencies already required 1.85. It is
  now 1.85, and CI builds on exactly the declared version so the two can't
  drift again.
- A FASTQ whose final record has no trailing newline no longer loses that
  record on the parallel path — EOF terminates the last line just as a
  newline does.
- `--quiet` now also suppresses the note printed when an input falls back to
  the tolerant per-record FASTQ reader.
- Adapter Content falls back to a direct scan if the pattern matcher cannot
  be built (reachable with a very large custom `--adapters` list) rather than
  silently reporting 0% adapter content and a PASS.

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
