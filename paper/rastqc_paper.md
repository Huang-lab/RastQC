# RastQC: A High-Performance Sequencing Quality Control Tool Written in Rust

**Kuan-Lin Huang**

## Abstract

Quality control (QC) of high-throughput sequencing data is a critical first step in genomics analysis pipelines. FastQC has served as the de facto standard for sequencing QC for over a decade, but its Java runtime dependency introduces startup overhead, elevated memory consumption, and deployment complexity. Meanwhile, the growing adoption of long-read sequencing platforms (Oxford Nanopore, PacBio) has created demand for QC tools that handle both short and long reads, yet existing solutions require separate tools for each data type and a third tool (MultiQC) to aggregate results.

Here we present RastQC, a unified sequencing QC tool written in Rust that combines FastQC-compatible short-read QC, long-read-specific metrics, built-in multi-sample summary, native MultiQC JSON export, and a web-based report viewer in a single 2.1 MB static binary. RastQC implements all 12 standard FastQC modules with matching algorithms, plus 3 long-read modules (Read Length N50, Quality Stratified Length, Homopolymer Content), achieving 100% module-level concordance with FastQC (55/55 calls identical across five model organisms). RastQC's streaming parallel pipeline with adaptive batch sizing delivers 1.8--3.2x speedup on short-read Illumina data and 4.7--6.5x speedup on long-read ONT/PacBio data, while using 8--9x less memory on small files and comparable memory on large files. RastQC is freely available at https://github.com/Huang-lab/RastQC under the MIT license.

## Introduction

Next-generation sequencing (NGS) has become the foundation of modern genomics research, with applications spanning whole-genome sequencing, RNA-seq, ChIP-seq, and single-cell assays. Before downstream analysis, quality control of raw sequencing data is essential to identify technical artifacts including base-calling errors, adapter contamination, GC bias, sequence duplication, and position-dependent quality degradation (Andrews, 2010; Patel & Jain, 2012). Failure to detect these issues can propagate systematic errors into variant calls, expression quantification, and other analyses.

FastQC (Andrews, 2010) has been the most widely used sequencing QC tool for over a decade, offering 12 diagnostic modules covering base quality, sequence content, duplication, adapter contamination, and k-mer enrichment. FastQC produces self-contained HTML reports that have become a standard deliverable in sequencing facilities worldwide. Tools such as MultiQC (Ewels et al., 2016) aggregate FastQC outputs across samples for project-level QC review.

The sequencing landscape has evolved significantly since FastQC's introduction. Long-read platforms from Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) now routinely generate reads spanning thousands to tens of thousands of bases, with distinct error profiles that require specialized QC metrics. Existing long-read QC tools such as NanoPlot (De Coster et al., 2018), LongQC (Fukasawa et al., 2020), PycoQC (Leger & Leonardi, 2019), and MinIONQC (Lanfear et al., 2019) address this need but are platform-specific and require separate installation, execution, and result interpretation alongside FastQC for short-read data. The recently published Sequali (Vorderman, 2025) represents progress toward unified short/long-read QC but lacks FastQC output compatibility and built-in multi-sample aggregation.

Meanwhile, FastQC itself has practical limitations rooted in its Java implementation. The Java Virtual Machine (JVM) imposes a 2--3 second startup penalty per invocation, a minimum memory footprint of ~300 MB regardless of input size, and a runtime dependency that complicates deployment in minimal container images and heterogeneous HPC environments. Alternative implementations such as falco (de Sena Brandine & Smith, 2021) address performance but do not extend functionality beyond FastQC's original 11 modules.

Here we present RastQC, a complete reimplementation and extension of FastQC in Rust. To our knowledge, RastQC is the first single-binary tool that unifies FastQC-compatible short-read QC modules, long-read-specific metrics, built-in multi-sample summary, native MultiQC JSON export, and a web-based report viewer. We demonstrate that RastQC achieves 1.8--6.5x speedup over FastQC across both short-read and long-read datasets while maintaining 100% concordance with FastQC on all shared modules.

## Implementation

### Architecture

RastQC is organized into six components reflecting the analysis pipeline:

1. **I/O layer** (`io/`): Streaming parsers for FASTQ (plain, gzip, bzip2-compressed), BAM/SAM formats via the noodles library, Oxford Nanopore Fast5 (HDF5) and POD5 (Apache Arrow IPC) formats (feature-gated), SOLiD colorspace auto-detection and decoding, and standard input streaming. Sequences are processed one at a time to minimize memory footprint.

2. **QC modules** (`modules/`): Fifteen independent analysis modules, each implementing a common `QCModule` trait that defines `process_sequence()`, `calculate_results()`, merge support for parallel processing, and output generation methods. Modules accumulate per-position statistics during the streaming pass and compute final results lazily. Per-position arrays are capped at 1,000 bases to bound memory usage on long reads while preserving resolution at the diagnostically important 5' end. All modules support accumulator-state merging for intra-file parallelism.

3. **Configuration** (`config/`): Embedded default adapter sequences (6 entries), contaminant sequences (15 entries), and pass/warn/fail thresholds matching FastQC defaults. Long-read modules (Read Length N50, Quality Stratified Length, Homopolymer Content) are disabled by default to avoid false positives on short-read data, and enabled via `--long-read` or automatically when processing Fast5/POD5 files. All configurations are overridable via command-line flags.

4. **Report generation** (`report/`): Self-contained HTML reports with inline SVG charts featuring complete axis labels and tick marks, tab-separated data files compatible with MultiQC, native MultiQC JSON output (`--multiqc-json`), per-file `summary.txt` with module-level status, and ZIP archives. A multi-file summary dashboard (HTML + TSV) is generated when processing multiple samples, providing MultiQC-like functionality without additional software.

5. **Web GUI** (`gui/`): A built-in HTTP server (`--serve`) provides a browser-based interface for navigating and viewing reports, with auto-browser launch and multi-sample summary access.

6. **Streaming parallel pipeline** (`parallel/`): For files >50 MB, a dedicated reader thread streams sequence batches through a bounded crossbeam channel to N worker threads, each with independent module instances. After the file is fully read, worker states are merged via `merge_from()`. Batch size is computed adaptively from early read lengths, targeting ~4 MB per batch. This prevents memory blowup on long reads (where fixed 16K-read batches would consume hundreds of megabytes) while maintaining high throughput on short reads. QC-aware exit codes (`--exit-code`) support automated pipeline gates.

### QC Modules

RastQC implements all 12 FastQC modules with matching algorithms, plus 3 long-read QC modules (Table 1).

**Table 1.** QC modules implemented in RastQC. Modules 13--15 are RastQC-exclusive long-read modules, disabled by default and enabled via `--long-read`.

| # | Module | Description | Key Algorithm |
|---|--------|-------------|---------------|
| 1 | Basic Statistics | Sequence count, length, %GC, total bases, encoding | Streaming counters; Phred encoding auto-detection |
| 2 | Per Base Sequence Quality | Quality distribution at each read position | Per-position quality histograms with adaptive base grouping |
| 3 | Per Tile Sequence Quality | Tile-specific quality deviations (Illumina) | Per-tile quality accumulation with 10% sampling after 10K reads |
| 4 | Per Sequence Quality Scores | Distribution of mean quality per read | Per-read average quality histogram |
| 5 | Per Base Sequence Content | A/T/G/C proportions per position | Per-position base counters with adaptive grouping |
| 6 | Per Sequence GC Content | GC% distribution vs. theoretical normal | GCModel fractional binning; bidirectional mode averaging; normal distribution fit |
| 7 | Per Base N Content | Unknown base frequency per position | Per-position N counters |
| 8 | Sequence Length Distribution | Read length variability | Adaptive binning maintaining ~50 display categories |
| 9 | Sequence Duplication Levels | Library complexity estimate | String-based tracking (first 50 bp); 100K unique cutoff with iterative binomial correction |
| 10 | Overrepresented Sequences | High-frequency sequences with contaminant matching | Exact/approximate string matching; 1-mismatch tolerance for >20 bp |
| 11 | Adapter Content | Known adapter contamination by position | Substring matching of 12-mer adapters with cumulative position counting |
| 12 | K-mer Content | Positionally biased k-mers | 7-mer tracking with 2% sampling; binomial enrichment test |
| 13 | Read Length N50 (Long Read) | Assembly-style length statistics | Sorted-length N50/N90; streaming min/max/mean/median |
| 14 | Quality Stratified Length (Long Read) | Length distribution by quality tier | 5-tier quality binning (Q<10 through Q40+); per-tier base counting |
| 15 | Homopolymer Content (Long Read) | Systematic homopolymer error detection | Run-length encoding; per-base run tracking (3--22 bp) |

Long-read module thresholds are calibrated for long-read error profiles. Homopolymer Content uses 5%/10% warn/fail thresholds for bases in runs of 3+, which is appropriate for long reads but would produce false positives on short-read Illumina data where ~18--20% of bases naturally occur in homopolymer runs. This motivated the default-off design.

### Unified QC Workflow

A typical multi-platform sequencing project requires running FastQC on short reads, a separate tool (e.g., NanoPlot, LongQC, or PycoQC) on long reads, and MultiQC to aggregate the results. RastQC consolidates this into a single command:

```bash
# Short-read QC (12 modules, FastQC-compatible output)
rastqc -o results/ *.fastq.gz

# Long-read QC (all 15 modules, auto-detected for Fast5/POD5)
rastqc --long-read -o results/ ont_reads.fastq.gz pacbio_reads.fastq.gz

# Multi-sample summary generated automatically
# Native MultiQC JSON available via --multiqc-json
```

The built-in summary dashboard provides immediate cross-sample comparison without requiring MultiQC installation, while the FastQC-compatible output format ensures full backward compatibility with existing MultiQC-based pipelines.

### MultiQC Compatibility

RastQC output is designed for full compatibility with MultiQC's FastQC parser. The `fastqc_data.txt` file uses identical module names, column headers, and data formats. Additionally, the `--multiqc-json` flag produces structured JSON output that eliminates parsing overhead for richer metadata transfer. Key compatibility features include:

- **Filename field** populated in Basic Statistics for correct sample identification
- **Total Bases** metric included for modern MultiQC versions
- **summary.txt** populated with per-module PASS/WARN/FAIL status and filename
- **ZIP archive structure** matching the expected `{sample}_fastqc/fastqc_data.txt` path
- **Module names** exactly matching FastQC conventions

## Methods

### Benchmarking Environment

All benchmarks were performed on a single workstation:

- **CPU**: Apple M1 Ultra (20 cores)
- **Memory**: 128 GB
- **Storage**: SSD
- **OS**: macOS 15.7.4 (Darwin 24.6.0, ARM64)
- **RastQC**: v0.1.0, compiled with Rust (release profile, LTO enabled, opt-level 3)
- **FastQC**: v0.12.1, OpenJDK
- **Threads**: 4

### Short-Read Datasets

We benchmarked on real human whole-exome sequencing data from the European Nucleotide Archive (Table 2).

**Table 2.** Short-read datasets used for benchmarking.

| File | Accession | Reads | Read Length | File Size |
|------|-----------|-------|-------------|-----------|
| DRR609229 R1 | DRR609229 | 719,828 | 76 bp | 22 MB |
| DRR609229 R2 | DRR609229 | 719,828 | 76 bp | 23 MB |
| ERR5897746 R1 | ERR5897746 | 4,252,217 | 126 bp | 320 MB |
| ERR5897746 R2 | ERR5897746 | 4,252,217 | 126 bp | 327 MB |
| DRR013000 R1 | DRR013000 | 24,778,423 | 76 bp | 1,430 MB |

### Long-Read Datasets

We additionally benchmarked on long-read bacterial sequencing data to evaluate performance on ONT and PacBio platforms (Table 3).

**Table 3.** Long-read datasets used for benchmarking.

| File | Platform | Accession | Reads | Mean Length | File Size |
|------|----------|-----------|-------|-------------|-----------|
| E. coli ONT | MinION | DRR242198 | 75,766 | 5,347 bp | 406 MB |
| E. coli PacBio | Revio | DRR723651 | 41,996 | 18,814 bp | 281 MB |

### Metrics

- **Wall-clock time**: Elapsed time measured by `/usr/bin/time -l`
- **Peak memory (RSS)**: Maximum resident set size reported by `/usr/bin/time -l`
- **Output concordance**: Module-level PASS/WARN/FAIL agreement on identical input
- **Per-step timing**: RastQC's `--time` flag provides QC, report generation, and I/O write timing breakdown

All benchmarks used 4 threads. RastQC's streaming parallel pipeline activates automatically for files >50 MB.

## Results

### Short-Read Performance

RastQC consistently outperformed FastQC across all short-read datasets (Table 4).

**Table 4.** Short-read performance comparison (4 threads). RastQC runs 12 core modules; FastQC runs 11.

| File | Size | Reads | FastQC Time | RastQC Time | Speedup | FastQC RSS | RastQC RSS |
|------|------|-------|-------------|-------------|---------|------------|------------|
| DRR609229 R1 | 22 MB | 720K | 3.50 s | **1.99 s** | 1.8x | 425 MB | 49 MB |
| DRR609229 R2 | 23 MB | 720K | 3.51 s | **2.04 s** | 1.7x | 424 MB | 50 MB |
| ERR5897746 R1 | 320 MB | 4.3M | 15.56 s | **4.81 s** | 3.2x | 446 MB | 332 MB |
| ERR5897746 R2 | 327 MB | 4.3M | 15.57 s | **4.81 s** | 3.2x | 442 MB | 330 MB |
| DRR013000 R1 | 1.4 GB | 24.8M | 51.75 s | **19.62 s** | 2.6x | 434 MB | 315 MB |
| All 5 files | 2.1 GB | 34.7M | 55.74 s | **22.25 s** | 2.5x | 1,068 MB | 825 MB |

RastQC achieved 1.7--3.2x speedup across all file sizes. Speedup is greatest on medium-sized files (320 MB, 3.2x) where the streaming parallel pipeline provides maximum benefit. On small files (22 MB), the 1.8x speedup reflects eliminated JVM startup overhead. On the largest file (1.4 GB), RastQC is 2.6x faster.

For small files processed sequentially (< 50 MB), RastQC uses 8--9x less memory (49 MB vs. 425 MB). For larger files processed in parallel, memory usage is comparable to FastQC due to per-worker module state duplication. Across all sizes, RastQC memory remains at or below FastQC.

### Long-Read Performance

RastQC showed even greater speedup on long-read data (Table 5).

**Table 5.** Long-read performance comparison (4 threads).

| File | Platform | Size | Reads | Mean Len | FastQC | RastQC | RastQC --long-read | Speedup |
|------|----------|------|-------|----------|--------|--------|--------------------|---------|
| E. coli ONT | MinION | 406 MB | 76K | 5.3 kb | 14.56 s | **3.12 s** | **3.13 s** | 4.7x |
| E. coli PacBio | Revio | 281 MB | 42K | 18.8 kb | 17.57 s | **2.72 s** | **2.71 s** | 6.5x |

RastQC achieved 4.7x speedup on ONT and 6.5x on PacBio data. The `--long-read` flag, which enables 3 additional QC modules, adds negligible overhead (<1%). The greater speedup on long reads compared to short reads reflects RastQC's adaptive batch sizing, which automatically reduces batch size for long reads to prevent memory blowup while maintaining throughput.

On PacBio data, RastQC uses slightly less memory than FastQC (670 MB vs. 702 MB). On ONT data with very long reads (max 100 kb), RastQC uses more memory (1,257 MB vs. 854 MB) due to the parallel pipeline's per-worker module state, though per-position array capping at 1,000 bases prevents this from growing further with read length.

### Long-Read QC Module Results

The long-read modules produced biologically meaningful results on both platforms:

- **ONT MinION**: N50 = 10,807 bp, median = 3,347 bp. 100% of reads below Q20, with 84% in the Q10--19 tier, consistent with expected ONT accuracy. Quality Stratified Length flagged a warning. Homopolymer Content detected 19.6% of bases in runs of 3+.

- **PacBio Revio**: N50 = 18,830 bp, median = 16,663 bp. 99% of reads in the Q30--39 tier, reflecting PacBio HiFi accuracy. Quality Stratified Length passed. Homopolymer Content detected systematic runs consistent with genomic content.

### Resource Comparison

**Table 6.** Deployment and resource comparison.

| Metric | RastQC | FastQC (Java) |
|--------|--------|---------------|
| Binary size | 2.1 MB | ~15 MB (JARs) |
| Runtime dependency | None (static binary) | Java 11+ (~200 MB) |
| Total deployment size | 2.1 MB | ~215 MB |
| Startup time | <5 ms | ~2.5 s |
| Peak memory (small files) | 49--50 MB | 424--425 MB |
| Peak memory (1.4 GB short-read) | 315 MB | 434 MB |
| Peak memory (long reads) | 670--1,257 MB | 702--854 MB |
| Modules | 12 core + 3 long-read | 11 |
| Multi-sample summary | Built-in (HTML + TSV) | Requires MultiQC |
| Native MultiQC JSON | Yes | No |
| Web report viewer | Built-in (`--serve`) | No |
| Long-read QC | Yes (`--long-read`) | No |
| Input formats | FASTQ, BAM, SAM, Fast5, POD5, stdin | FASTQ, BAM, SAM |

### Output Concordance

We systematically compared module-level PASS/WARN/FAIL calls between RastQC and FastQC on all short-read benchmark files plus five model organism datasets from the European Nucleotide Archive (*E. coli*, *S. cerevisiae*, *D. melanogaster*, *M. musculus*, *H. sapiens*).

All 11 shared modules produced identical PASS/WARN/FAIL calls across all organisms, achieving **100% concordance** (55/55 module comparisons). Key algorithmic details ensuring concordance:

1. **Per Sequence GC Content**: RastQC implements FastQC's GCModel for fractional bin smoothing and uses `sum(gcDistribution)` (not raw sequence count) as the total count for standard deviation and theoretical curve computation, matching FastQC's behavior.

2. **Sequence Duplication Levels**: RastQC uses FastQC's exact iterative binomial correction formula for estimating true duplication counts beyond the 100K unique sequence observation window.

3. **Overrepresented Sequences**: RastQC uses the total sequence count (not count-at-limit) as the denominator for percentage calculations, matching FastQC's behavior.

## Discussion

### Comparison with Existing Tools

Several tools address sequencing QC, each with different scope and trade-offs (Table 7).

**Table 7.** Feature comparison of sequencing QC tools.

| Feature | RastQC | FastQC | falco | fastp | Sequali | NanoPlot | LongQC |
|---------|--------|--------|-------|-------|---------|----------|--------|
| Language | Rust | Java | C++ | C++ | C/Python | Python | Python/C |
| Short-read QC | Yes | Yes | Yes | Yes | Yes | No | No |
| Long-read QC | Yes | No | No | No | Yes | Yes | Yes |
| FastQC-compatible output | Yes | -- | Yes | No | No | No | No |
| Multi-sample summary | Built-in | No | Minimal | No | No | No | No |
| Native MultiQC JSON | Yes | No | No | fastp JSON | Yes | No | No |
| Web report viewer | Built-in | No | No | No | No | No | No |
| Single binary, no deps | Yes | No (Java) | Yes | Yes | No (Python) | No | No |
| Per-file timing (`--time`) | Yes | No | No | No | No | No | No |
| Exit codes for pipelines | Yes | No | No | No | No | No | No |

**FastQC** (Andrews, 2010) remains the standard for short-read QC. RastQC is a drop-in replacement that adds speed, lower memory, long-read support, and built-in aggregation while maintaining full output compatibility.

**falco** (de Sena Brandine & Smith, 2021) is the closest short-read alternative, reimplementing FastQC in C++ with ~3x speedup. However, falco lacks long-read modules, a web GUI, native MultiQC JSON, and multi-sample summary. Its documented `--nano` flag is not implemented.

**fastp** (Chen et al., 2018) combines QC reporting with read preprocessing (trimming, filtering, adapter removal). While highly popular, fastp does not produce FastQC-compatible output, has no long-read support (a separate tool "fastplong" exists), and lacks multi-sample aggregation.

**Sequali** (Vorderman, 2025) is the most recent competitor and the closest in scope, offering both short- and long-read QC with MultiQC integration. However, Sequali does not produce FastQC-compatible output (preventing drop-in replacement in existing pipelines), lacks built-in multi-sample aggregation, has no web GUI, and requires Python installation rather than distributing as a single binary.

**NanoPlot** (De Coster et al., 2018), **LongQC** (Fukasawa et al., 2020), **PycoQC** (Leger & Leonardi, 2019), and **MinIONQC** (Lanfear et al., 2019) provide platform-specific long-read QC but do not address short-read data, requiring users to maintain separate tool chains.

RastQC uniquely combines FastQC-compatible output, long-read metrics, built-in multi-sample summary, native MultiQC JSON, a web report viewer, and single-binary deployment in one tool, eliminating the need for separate short-read QC, long-read QC, and result aggregation tools.

### Performance Characteristics

RastQC's performance advantage stems from three architectural choices:

1. **Native binary**: Eliminates 2.5 s JVM startup per invocation. For batch processing of 1,000 amplicon panel samples (~50K reads each), this alone saves ~42 minutes.

2. **Streaming parallel pipeline**: Adaptive batch sizing targets ~4 MB per batch regardless of read length, preventing memory blowup on long reads while maintaining high throughput. On PacBio data (mean 18.8 kb reads), this yields 6.5x speedup over FastQC.

3. **Memory-bounded design**: Per-position arrays are capped at 1,000 bases, hash-map-based modules (Duplication, Overrepresented) freeze after 100K unique sequences, and channel buffers are bounded. This keeps memory usage predictable across diverse inputs.

### Long-Read QC Module Design

The three long-read modules were designed with calibrated thresholds based on platform-specific error profiles. Homopolymer Content uses 5%/10% (warn/fail) thresholds for the fraction of bases in runs of 3+ bases. In random DNA, ~16% of bases occur in such runs, and real short-read Illumina data shows 18--20%. These modules are therefore disabled by default and enabled only via `--long-read` or when processing native long-read formats (Fast5/POD5), preventing false positives on short-read data.

### Limitations

RastQC's memory usage on large files processed in parallel is comparable to FastQC rather than substantially lower, due to per-worker module state duplication. Future work could explore shared-memory module designs or more aggressive state compression to maintain the small-file memory advantage at scale. Additionally, while the 1,000-position cap on per-base tracking is sufficient for most short-read and the diagnostically relevant 5' portion of long reads, users analyzing position-specific quality across entire ultra-long reads (>100 kb) would need to adjust this parameter.

## Availability

- **Source code**: https://github.com/Huang-lab/RastQC
- **License**: MIT
- **Language**: Rust (2021 edition)
- **Installation**: `cargo install --path .` or pre-compiled binaries
- **Test suite**: 36 tests (25 unit, 11 integration)
- **System requirements**: Any platform supported by Rust (Linux, macOS, Windows)
- **Optional features**: `--features nanopore` for Fast5/POD5 support (requires HDF5 system library)

## References

Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884--i890.

De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15), 2666--2669.

de Sena Brandine, G., & Smith, A. D. (2021). Falco: high-speed FastQC emulation for quality control of sequencing data. *F1000Research*, 8, 1874.

Ewels, P., Magnusson, M., Lundin, S., & Kaller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047--3048.

Fukasawa, Y., Ermini, L., Wang, H., Carty, K., & Cheung, M. S. (2020). LongQC: A Quality Control Tool for Third Generation Sequencing Long Read Data. *G3: Genes, Genomes, Genetics*, 10(4), 1193--1196.

Lanfear, R., Schalamun, M., Kainer, D., Wang, W., & Schwessinger, B. (2019). MinIONQC: fast and simple quality control for MinION sequencing data. *Bioinformatics*, 35(3), 523--525.

Leger, A., & Leonardi, T. (2019). pycoQC, interactive quality control for Oxford Nanopore Sequencing. *Journal of Open Source Software*, 4(34), 1236.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094--3100.

Patel, R. K., & Jain, M. (2012). NGS QC Toolkit: a toolkit for quality control of next generation sequencing data. *PLoS ONE*, 7(2), e30619.

Vorderman, R. H. P. (2025). Sequali: sequence quality metrics for short and long reads. *Bioinformatics Advances*, 5(1), vbaf010.
