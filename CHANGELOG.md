# Changelog

All notable changes to RastQC are documented here.

## [Unreleased]

### Performance

Measured against 0.2.0 on the same machine (Intel Core i9-9900K, 8C/16T,
macOS 15 x86-64) on 2026-09-24: 3 repetitions per row, fastest wall time, peak
RSS across the same runs via `/usr/bin/time -l`. The harness changed which
statistic it reports partway through this cycle (see Changed); every row below,
and every row in [`benchmark/RESULTS.md`](benchmark/RESULTS.md) that is not
marked 0.1.0, was re-measured in one sitting under the current one. Reported
values are unchanged: `fastqc_data.txt` from a 0.2.0 `-t 1` run matches this
build at `-t 4` and `-t 16` byte for byte, on both a short-read and a long-read
file.

Short-read, ERR5897746_1 (4.3M reads, 126 bp, 320 MB gzipped):

| | 0.2.0 | this build |
|---|---|---|
| wall, `-t 1` | 4.14 s | **3.74 s** |
| wall, `-t 4` | 2.53 s | **2.28 s** |
| wall, `-t 16` | 2.49 s | **2.31 s** |
| system time, `-t 4` | 0.27 s | **0.09 s** |
| peak RSS, `-t 1` | 126 MB | **72 MB** |
| peak RSS, `-t 4` | 149 MB | **73 MB** |
| peak RSS, `-t 16` | 147 MB | **75 MB** |

Long-read, DRR242198_1 (75.8k ONT reads, 5.9 kb mean, 100.6 kb max, 406 MB):

| | 0.2.0 | this build |
|---|---|---|
| wall, `-t 4` | 2.98 s | **2.39 s** |
| peak RSS, `-t 1` | 460 MB | **422 MB** |
| peak RSS, `-t 4` | 1045 MB | **410 MB** |
| peak RSS, `-t 16` | 1021 MB | **415 MB** |

- **The block reader now recycles its buffers.** It allocated a fresh 1 MB
  `Vec` for every block and dropped it once both consumers were done, so each
  megabyte of input cost an `mmap` and a `madvise` teardown — which profiled
  as the single largest cost in a run, above every QC module. Finished blocks
  are now returned to the reader and refilled. This is where the 3x drop in
  system time and roughly half the peak memory come from.
- **The block buffer no longer reallocates on every block.** `next_block`
  fills in 128 KB chunks and stops once it has passed `BLOCK_SIZE`, so the
  last chunk pushed the length just past a capacity of exactly `BLOCK_SIZE` —
  a realloc and a 1 MB copy per block. Buffers now have room for the
  overshoot.
- **Long-read runs no longer get worse as `-t` rises.** The per-file worker
  cap 0.2.0 introduced keyed off the filename extension alone, so long-read
  input still got 4 workers (gzipped) or 8 (plain). Per sequence GC content
  keeps one GC model per distinct read length and a model for length *L*
  holds *L+1* inner vectors, so every worker built its own set — Illumina has
  one length and one model, an ONT run buckets to ~100 lengths and millions
  of vectors — while the reader thread stayed the ceiling. At `-t 4` that was
  *slower* than `-t 1` (2.98 s vs 2.53 s) for 2.3x the memory. The cap now
  samples mean read length from the first block and gives such files one
  worker; the rest of the budget goes to other files. A 282 MB PacBio run
  drops from 465 MB to 198 MB at `-t 4`. Short-read files are untouched and
  still scale with `-t`. The deeper fix is to share those models between
  workers rather than rebuild them.
- **Quality Stratified Length** widened every base to `f64` and accumulated
  the read's quality sum in floating point. It now sums as an integer and
  converts once — exact rather than accumulating rounding error over hundreds
  of thousands of bases, and faster: `--long-read` on a 282 MB PacBio run
  drops from 3.84 s to 3.38 s, and on a 406 MB ONT run from 3.08 s to 2.81 s.
  This loop matters more than most per-base arithmetic because it walks the
  whole read, where the per-position modules stop at 1000 bases.
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
- A CI job that syntax-checks and shellchecks the benchmark scripts. They
  produce the numbers in `RESULTS.md` and the paper, but nothing checked them,
  and two breakages this cycle would have been caught by it.
- [`benchmark/DATA.md`](benchmark/DATA.md) — what every benchmark and
  validation dataset actually is: accession, organism, platform, and read
  counts and lengths *measured from the files* rather than taken from ENA
  metadata, whose `base_count` is per spot for some runs and per file for
  others. Dividing those two fields is how the README came to describe an
  18.7M-read NextSeq run as 144 bp when it is 72 bp. It also records that
  three files in the 0.1.0 paper's dataset are not what their names say — a
  second mouse run named `fly`, SARS-CoV-2 amplicon data named `yeast`, and a
  196-byte failed download named `zebrafish` — so the paper's "five model
  organisms" were three, with mouse counted twice. Which module calls match
  FastQC does not depend on the species, so the concordance result stands;
  only the description of the validation set was wrong.

### Changed

- `benchmark/fetch_data.sh` now fetches every dataset the benchmarks and the
  paper cite, verifies each against the byte size ENA reports, and downloads
  large files as concurrent byte ranges that resume from their partial offset.
  ENA throttles a single connection to roughly 200 KB/s, which is over five
  hours for the full set; twelve ranges measured ~2 MB/s.
- `benchmark/run_benchmark.sh` classifies long-read inputs by their measured
  mean read length instead of by filename. The previous `*_ont_*`/`*_pacbio_*`
  convention silently benchmarked real ENA files as short-read, and the script
  now emits stable tool ids that `paper/analyze_benchmarks.py` matches on —
  the two had drifted apart, so the figures were being generated from no data.
- Both benchmark scripts now run on bash 3.2, which is what macOS ships.
- `benchmark/run_benchmark.sh` reports the **fastest** of its repetitions
  rather than the median. Timing interference on this machine is one-sided and
  occasionally outlasts two of three repetitions — one FastQC measurement came
  out 47.7 s / 242 s / 47.7 s, and isolated repeats of the 242 s outlier gave
  47.7 s — which a median cannot survive and a minimum can. Peak RSS over the
  same runs agreed within 15% on 42 of 45 measurements, so only wall time was
  affected. The paper's figures are regenerated from the CSV and so already use
  it; `benchmark/RESULTS.md` has now been re-measured under it too, and marks
  the RastQC 0.1.0 rows it did not re-run.
- `benchmark/fetch_data.sh` gained a `human` group: one full 30X human WGS run
  (ERR3239334, 12 GB). It is the only group that takes hours rather than
  minutes, so `short long` remains the way to cover every platform and three
  orders of magnitude of size without it. `all` now includes it, and is ~16 GB.
- The two Bioconda recipes are now one multi-output recipe at
  [`recipes/rastqc-meta/`](recipes/rastqc-meta/), mirroring
  [bioconda-recipes#69452](https://github.com/bioconda/bioconda-recipes/pull/69452).
  `rastqc` and `rastqc-nanopore` are the same program from the same tag with
  different cargo features and are always released together, so two recipes
  meant two version bumps and two `sha256` edits per release — 0.2.0 produced
  two autobump PRs plus one maintainer PR for one tag. The staging copy here is
  byte-identical to what is live, which is the check that 0.2.0's dropped
  `additional-platforms` taught us to run first.

### Fixed

- The agent skill (`RastQC.md`) told users to install with "requires Rust
  1.70+", and its Nextflow gate example tested `[ $? -eq 2 ]` — which passes
  an unreadable or malformed file through the gate, the exact case exit code 3
  exists to stop.
- `benchmark/run_benchmark.sh` looked up per-file stats with a tab-delimited
  parse of a space-delimited table, so every lookup returned nothing: the CSV
  that feeds `RESULTS.md` and the paper's figures got empty `size_mb`, `reads`
  and `mean_read_len` columns, the per-file headers printed blanks, and
  `$((total_mb + $(stat_of ...)))` became an arithmetic syntax error that
  killed the aggregate-group rows outright under `set -e`. Both the table and
  the lookup are now tab-delimited, which was the point of the change — a
  filename containing a space stays in one field.
- `benchmark/run_benchmark.sh` skips an input it cannot read as FASTQ instead
  of benchmarking it. A truncated or non-gzip file still lets `awk` reach `END`
  and report zero reads, which memoized a 0 bp mean under a cache key of
  `basename:size` — permanently classifying the file as short-read, so its
  long-read rows silently vanished — and timed each tool's error path as if it
  were a wall time.
- `benchmark/fetch_data.sh` records the chunk layout a partial download was
  written under, and refuses to resume chunks that carry no recorded layout
  rather than guessing. Chunk boundaries derive from `$PARALLEL`; resuming
  under a different value appends bytes from the new layout onto bytes from the
  old one, and every chunk still reaches its expected length, so the total size
  check passes and the file is silently corrupt in the middle. Which layout
  produced an unlabelled chunk dir cannot be recovered from the chunks — a
  chunk that was complete under a larger `$PARALLEL` is shorter than its range
  under a smaller one, and so is indistinguishable from a partial one.
- Read lengths in the documentation, all now measured from the files:
  DRR045135_1 is 72 bp (was 144), DRR048760 is 67 bp (was 76, a
  transposition), and the paper's dataset tables described DRR609229 as a fixed
  76 bp when it is variable (83/84 bp mean, 35–151), the ONT run as 5,347 bp
  when it is 5,920, and the PacBio run as 18,814 bp when it is 17,609 — that
  column had been filled in with the N50. The N50 and median figures in the
  paper's text were checked against the tool's own output and are correct.

## [0.2.0] — 2026-09-10

Performance, memory and reproducibility release. Reported QC values are
unchanged from a 0.1.0 **sequential** run; see "Output" below for the two
places where 0.1.0's *parallel* runs disagreed with it.

### Performance

Prompted by [#12](https://github.com/Huang-lab/RastQC/issues/12), which
reported RastQC running ~2× slower than Falco on NextSeq runs. That
reproduced: on a public 18.7M-read NextSeq 500 run, 0.1.0 at `-t 1` took
64.9 s against Falco's 32.6 s. 0.2.0 takes 11.3 s at `-t 1` and 6.9 s at
`-t 4`. Those are the median of 3 that the harness reported at the time;
re-measured as the fastest of 3, Falco is 30.55 s and 0.2.0 is 10.43 s and
7.10 s, and 0.1.0 has not been re-run. Full numbers and method in
[`benchmark/RESULTS.md`](benchmark/RESULTS.md).

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
