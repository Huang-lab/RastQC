# RastQC

A fast quality control tool for high-throughput sequencing data, written in Rust. Drop-in replacement for [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with identical QC modules, matching algorithms, and compatible output formats.

## Features

- **15 QC modules**: all 12 FastQC modules + 3 long-read QC modules
- **Fast**: 2.9x faster than Falco single-threaded, 4.7x with 4 threads, on a real 18.7M-read NextSeq run ([benchmarks](benchmark/RESULTS.md))
- **Portable**: single 2.6 MB static binary, no Java runtime needed
- **Compatible output**: HTML reports, tab-separated data files, ZIP archives, native MultiQC JSON
- **Multi-file summary**: overview dashboard when processing many files
- **Web GUI**: built-in report browser (`--serve`)
- **Input formats**: FASTQ, gzip, bzip2, BAM, SAM, SOLiD colorspace, Fast5/POD5 (optional), stdin
- **Pipeline integration**: QC-aware exit codes (`--exit-code`) for Nextflow/Snakemake gates

## Installation

### Via conda (Bioconda)

```bash
# Core build (short-read QC)
conda install -c bioconda rastqc

# With Fast5/POD5 support (Oxford Nanopore)
conda install -c bioconda rastqc-nanopore
```

Both packages are built from one multi-output Bioconda recipe, which lives in
[`recipes/rastqc-meta/`](recipes/rastqc-meta/).

### From source

```bash
# Requires Rust 1.85+
cargo install --path .
```

### Build manually

```bash
git clone https://github.com/Huang-lab/RastQC.git
cd RastQC
cargo build --release
# Binary at ./target/release/rastqc
```

### With Nanopore format support

```bash
cargo build --release --features nanopore
```

## Quick start

```bash
# Single file
rastqc sample.fastq.gz

# Multiple files (processed in parallel)
rastqc *.fastq.gz

# Specify output directory
rastqc -o results/ sample_R1.fastq.gz sample_R2.fastq.gz

# HTML only (no ZIP)
rastqc --nozip -o results/ sample.fastq.gz

# Stream from stdin (gzip/bzip2 auto-detected)
samtools fastq aligned.bam | rastqc --stdin -o results/
zcat sample.fastq.gz | rastqc --stdin -o results/

# Use 8 threads
rastqc -t 8 -o results/ *.fastq.gz

# Pipeline QC gate (exit 2 if any module fails)
rastqc --exit-code sample.fastq.gz || echo "QC failed"

# Browse reports in browser
rastqc -o results/ *.fastq.gz --serve

# Native MultiQC JSON output
rastqc --multiqc-json -o results/ sample.fastq.gz
```

## Usage

```
rastqc [OPTIONS] [FILES]...

Arguments:
  [FILES]...  Input files (FASTQ, FASTA, BAM, SAM, Fast5, POD5). Use "-" for stdin (gzip/bzip2 auto-detected).

Options:
  -o, --outdir <DIR>            Output directory [default: current directory]
  -t, --threads <N>             Number of threads [default: all CPUs]
  -c, --contaminants <FILE>     Custom contaminant list (tab-separated: name\tsequence)
  -a, --adapters <FILE>         Custom adapter list (tab-separated: name\tsequence)
  -l, --limits <FILE>           Custom pass/warn/fail thresholds
  -k, --kmer-size <N>           Kmer size for enrichment analysis [default: 7]
      --stdin                   Read FASTQ from standard input (gzip/bzip2 auto-detected)
      --nofilter                Include all reads (don't skip QC-failed reads)
      --extract                 Extract ZIP contents after creation
      --nozip                   Write HTML report only, skip ZIP archive
      --summary                 Write multi-file summary report
      --multiqc-json            Output native MultiQC JSON alongside standard reports
      --exit-code               Return QC-aware exit codes: 0=pass, 1=warn, 2=fail
      --serve                   Start web server to browse reports
      --port <N>                Web server port [default: 8080]
      --long-read               Enable long-read QC modules (auto-enabled for Fast5/POD5 inputs)
      --time                    Show per-file and per-step timing breakdown
      --no-parallel             Disable streaming intra-file parallelism (on by default for >50MB files)
  -q, --quiet                   Suppress progress output
      --dup-length <N>          Truncation length for duplication detection [default: 50]
  -h, --help                    Print help
  -V, --version                 Print version
```

## Architecture

```
rastqc/
├── src/
│   ├── main.rs              # CLI entry point, file dispatch, exit codes
│   ├── config.rs            # Adapters, contaminants, limits, thresholds
│   ├── gui.rs               # Built-in HTTP server for report browsing
│   ├── parallel.rs          # Streaming parallel pipeline (reader → channel → workers → merge)
│   ├── io/
│   │   ├── mod.rs           # SequenceReader enum (unified format dispatch)
│   │   ├── fastq.rs         # FASTQ/gz/bz2 streaming reader + stdin
│   │   ├── bam.rs           # BAM/SAM reader via noodles
│   │   ├── colorspace.rs    # SOLiD di-base → basespace decoder
│   │   ├── fast5.rs         # Oxford Nanopore Fast5 (HDF5) reader
│   │   └── pod5.rs          # Oxford Nanopore POD5 (Arrow IPC) reader
│   ├── modules/
│   │   ├── mod.rs           # QCModule trait, merge support, factory
│   │   ├── basic_stats.rs   # Sequence count, length, %GC, encoding
│   │   ├── per_base_quality.rs
│   │   ├── per_tile_quality.rs
│   │   ├── per_sequence_quality.rs
│   │   ├── per_base_content.rs
│   │   ├── per_sequence_gc.rs
│   │   ├── n_content.rs
│   │   ├── sequence_length.rs
│   │   ├── duplication.rs
│   │   ├── overrepresented.rs
│   │   ├── adapter_content.rs
│   │   ├── kmer_content.rs
│   │   └── long_read_quality.rs  # N50, quality-stratified length, homopolymer
│   └── report/
│       └── mod.rs           # HTML, text, JSON, ZIP, summary generation
├── tests/
│   └── integration_test.rs  # 16 integration tests
├── paper/                   # Manuscript, benchmarks, figures
└── FastQC/                  # Reference FastQC for concordance testing
```

**Data flow**: Files → `SequenceReader` → streaming `Sequence` records → each record passed to all `QCModule` instances → `calculate_results()` → report generation (HTML/text/JSON/ZIP).

**Streaming parallel pipeline** (default for files >50MB): a reader thread decompresses the file and cuts it into record-aligned blocks, which it hands to two consumers. A pool of N worker threads takes whichever block is next and parses records as slices borrowed from it, so parsing scales with the pool and the hot loop allocates nothing. A single in-order consumer runs the modules that need the file's whole read stream — Sequence Duplication Levels, Overrepresented sequences, Kmer Content and Per tile sequence quality, each of which either keeps a capped observation table or samples off a read counter, and so cannot be reconstructed by merging per-worker partials. Worker states are merged via `merge_from()` at the end. Nothing buffers the whole file, output is identical at any `-t`, and per-worker state no longer scales with thread count.

**Threads**: `-t` is a budget for the whole run. Files are analyzed concurrently and each file runs its own worker pool; the budget is split across the two levels rather than applied to both, and each file's pool is capped at the point where a single file stops benefiting (the reader thread becomes the ceiling). Passing a large `-t` therefore costs neither the thread explosion nor the memory it used to.

All 15 modules implement the `QCModule` trait with `process_sequence()`, `calculate_results()`, `merge_from()` (for parallel chunk merging), and output methods. Modules are created by `ModuleFactory` based on the limits configuration.

## Output files

For each input file `sample.fastq.gz`, RastQC produces:

| File | Description |
|------|-------------|
| `sample_fastqc.zip` | ZIP archive containing all outputs below |
| `sample_fastqc/fastqc_report.html` | Self-contained HTML report with SVG charts |
| `sample_fastqc/fastqc_data.txt` | Tab-separated data for each module |
| `sample_fastqc/summary.txt` | One-line PASS/WARN/FAIL per module |
| `sample_multiqc.json` | Native MultiQC JSON (with `--multiqc-json`) |

When processing multiple files with `--summary`:

| File | Description |
|------|-------------|
| `summary.tsv` | Tab-separated matrix: rows = files, columns = modules |
| `summary.html` | Overview dashboard linking to all individual reports |

## QC modules

| # | Module | What it checks | Pass/Warn/Fail criteria |
|---|--------|---------------|------------------------|
| 1 | **Basic Statistics** | Sequence count, length, %GC, encoding | Informational only |
| 2 | **Per Base Sequence Quality** | Quality score distribution at each position | Median < 25 (warn) / < 20 (fail) |
| 3 | **Per Tile Sequence Quality** | Quality variation between flowcell tiles | Max deviation > 5 (warn) / > 10 (fail) |
| 4 | **Per Sequence Quality Scores** | Distribution of mean quality per read | Mode <= 27 (warn) / <= 20 (fail) |
| 5 | **Per Base Sequence Content** | A/T/G/C proportions at each position | |A-T| or |G-C| > 10% (warn) / > 20% (fail) |
| 6 | **Per Sequence GC Content** | GC% distribution vs theoretical normal | Deviation > 15% (warn) / > 30% (fail) |
| 7 | **Per Base N Content** | Unknown base (N) frequency per position | N% > 5 (warn) / > 20 (fail) |
| 8 | **Sequence Length Distribution** | Read length variability | Variable lengths (warn) |
| 9 | **Sequence Duplication Levels** | Library complexity estimate | < 70% unique (warn) / < 50% unique (fail) |
| 10 | **Overrepresented Sequences** | Frequently occurring sequences + contaminant matching | Any seq > 0.1% (warn) / > 1% (fail) |
| 11 | **Adapter Content** | Known adapter sequence contamination | > 5% (warn) / > 10% (fail) |
| 12 | **Kmer Content** | Positionally biased k-mers | -log10(p) > 2 (warn) / > 5 (fail) |
| 13 | **Read Length N50** (Long Read) | N50, N90, mean, median, min, max lengths | Informational only |
| 14 | **Quality Stratified Length** (Long Read) | Length distribution by quality tier (Q<10 to Q40+) | >50% below Q20 (warn) |
| 15 | **Homopolymer Content** (Long Read) | Homopolymer run frequency by base and length | >5% bases in runs (warn) / >10% (fail) |

Modules 13--15 are RastQC-exclusive, designed for long-read sequencing data (PacBio HiFi, Oxford Nanopore). These modules are **disabled by default** and enabled with `--long-read` or automatically when processing Fast5/POD5 files. Their thresholds are calibrated for long-read error profiles and would produce false positives on short-read Illumina data.

## Working with many files

### Batch processing

```bash
# Process all FASTQ files in a directory
rastqc -o qc_results/ data/*.fastq.gz

# Process with summary dashboard
rastqc -o qc_results/ --summary data/*.fastq.gz

# Use find for recursive discovery
find data/ -name "*.fastq.gz" | xargs rastqc -o qc_results/ --summary
```

### Summary report

The `--summary` flag generates two files for multi-file review:

**`summary.tsv`** -- machine-readable matrix for scripting:
```
Sample	Basic Statistics	Per Base Quality	...	Adapter Content
sample_A	PASS	PASS	...	WARN
sample_B	PASS	FAIL	...	PASS
```

**`summary.html`** -- browser-friendly dashboard with color-coded PASS/WARN/FAIL table.

### Filtering results

```bash
# Find all failing samples
grep "FAIL" qc_results/summary.tsv

# Count warnings per sample
awk -F'\t' '{n=0; for(i=2;i<=NF;i++) if($i=="WARN") n++; print $1, n}' qc_results/summary.tsv
```

## Custom configuration

### Adapter list

Tab-separated file with adapter name and 12bp sequence:

```
My Custom Adapter	AGATCGGAAGAG
Another Adapter		CTGTCTCTTATA
```

### Contaminant list

Tab-separated file with contaminant name and full sequence:

```
PhiX Control	GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACT
Custom Primer	AATGATACGGCGACCACCGA
```

### Limits file

Controls pass/warn/fail thresholds and which modules run:

```
# Disable a module
kmer    ignore  1

# Adjust thresholds
quality_base_lower  warn    10
quality_base_lower  error   5
adapter             warn    5
adapter             error   10
```

## Compatibility with FastQC

RastQC produces output compatible with tools that consume FastQC results:

- **MultiQC**: `fastqc_data.txt` files are compatible with MultiQC's FastQC module
- **Native JSON**: `--multiqc-json` provides structured output without parsing
- **summary.txt**: same PASS/WARN/FAIL format per module
- **Identical module names and data headers** in text output
- **100% concordance**: 55/55 module calls identical across 5 model organisms

## Performance

Two comparisons, measured differently — read both.

### vs Falco (measured for this release)

Intel Core i9-9900K (8C/16T), macOS x86-64, median of 3 runs, peak RSS via
`/usr/bin/time -l`. Falco is single-threaded — its `-t` flag is documented in
its own help as "NOT YET IMPLEMENTED" — so `rastqc -t 1` is the like-for-like
row. Reproduce with
`./benchmark/fetch_data.sh nextseq && ./benchmark/run_benchmark.sh`.

| Dataset | Falco 1.2.5 | RastQC `-t 1` | RastQC `-t 4` |
|---------|-------------|---------------|---------------|
| DRR045135_1 — 18.7M reads, 72 bp, 828 MB | 32.6 s / 86 MB | **11.3 s** / 126 MB | **6.9 s** / 158 MB |
| DRR048760 — 1.2M reads, 67 bp, 54 MB | 2.32 s / 88 MB | **0.83 s** / 95 MB | **0.67 s** / 128 MB |
| Both files, one invocation | 35.1 s / 94 MB | — | **8.6 s** / 208 MB |

Full method, and what changed since 0.1.0, in
[`benchmark/RESULTS.md`](benchmark/RESULTS.md); what each dataset actually is,
in [`benchmark/DATA.md`](benchmark/DATA.md).

### vs FastQC

Measured on macOS ARM64 with 4 threads against RastQC 0.1.0 — different
hardware and an older RastQC than the table above, so the two are not directly
comparable. 0.2.0 is substantially faster than the RastQC column here.

#### Short-read (Illumina)

| File | Size | Reads | FastQC 0.12.1 | RastQC 0.1.0 | Speedup |
|------|------|-------|---------------|--------------|---------|
| DRR609229 R1 | 22 MB | 720K | 3.5s | **2.0s** | 1.8x |
| DRR609229 R2 | 23 MB | 720K | 3.5s | **2.0s** | 1.7x |
| ERR5897746 R1 | 320 MB | 4.3M | 15.6s | **4.8s** | 3.2x |
| ERR5897746 R2 | 327 MB | 4.3M | 15.6s | **4.8s** | 3.2x |
| DRR013000 R1 | 1.4 GB | 24.8M | 51.8s | **19.6s** | 2.6x |
| All 5 files | 2.1 GB | 34.7M | 55.7s | **22.3s** | 2.5x |

#### Long-read (ONT / PacBio)

| File | Platform | Size | Reads | Mean Length | FastQC | RastQC 0.1.0 | Speedup |
|------|----------|------|-------|-------------|--------|--------------|---------|
| DRR242198 | ONT MinION | 406 MB | 76K | 5.3 kb | 14.6s | **3.1s** | 4.7x |
| DRR723651 | PacBio Revio | 281 MB | 42K | 18.8 kb | 17.6s | **2.7s** | 6.5x |

The `--long-read` flag enables 3 additional QC modules with negligible overhead.

### Memory

Peak resident memory no longer grows with `-t`. Measured on the 18.7M-read
NextSeq run above, and on six 240 MB uncompressed FASTQs passed in one
invocation:

| Run | RastQC 0.1.0 | RastQC 0.2.0 |
|-----|--------------|--------------|
| 1 file (828 MB gz), `-t 4` | 322 MB / 20.9 s | **163 MB / 7.4 s** |
| 1 file (828 MB gz), `-t 16` | 1119 MB / 22.3 s | **146 MB / 6.9 s** |
| 6 files (240 MB each), `-t 8` | 4052 MB / 6.3 s | **546 MB / 1.6 s** |
| 6 files (240 MB each), `-t 16` | 5148 MB / 7.0 s | **586 MB / 1.3 s** |

0.1.0 got *slower* as `-t` rose past 4 while its memory kept climbing.

### Resource comparison

| Metric | RastQC | FastQC (Java) |
|--------|--------|---------------|
| Binary size | 2.6 MB | ~215 MB (with JRE) |
| Startup time | <5 ms | ~2.5 s JVM warmup |
| Peak memory (small files) | 49-50 MB | 424-425 MB |
| Threading | streaming intra-file + multi-file parallel, single `-t` budget | per-file parallel |
| Modules | 12 core + 3 long-read | 11 |

### Reproducibility

Output is byte-identical at any `-t`, and identical to a sequential run —
verified across five datasets including two real NextSeq runs. Reports are
also stable run to run: sorts that previously left tied rows in
`HashMap` order now break ties deterministically.

---

## Citation

If you use RastQC in your research, please cite:

> Huang KL. RastQC: A fast, Rust-based quality control tool for high-throughput sequencing data. *bioRxiv* (2026). [https://www.biorxiv.org/content/10.64898/2026.03.31.715630v2]([https://www.biorxiv.org/content/10.64898/2026.03.31.71563](https://www.biorxiv.org/content/10.64898/2026.03.31.715630v2))

## Acknowledgments

RastQC is a reimplementation inspired by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) by Simon Andrews at the Babraham Institute. FastQC has served as the gold standard for sequencing quality control for over a decade, and its elegant module design, diagnostic algorithms, and output formats are the foundation upon which RastQC is built. We are grateful to the FastQC team for creating and maintaining such an essential tool for the genomics community.

## License

MIT License. See [LICENSE](LICENSE) for details.

Contributions are welcome! Please open an issue or pull request on [GitHub](https://github.com/Huang-lab/RastQC).

## Author

Written by **Kuan-Lin Huang** at [PrecisionOmics.org](https://PrecisionOmics.org)
