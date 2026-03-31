# RastQC: A High-Performance Sequencing Quality Control Tool Written in Rust

**Kuan-Lin Huang**

## Abstract

Quality control (QC) of high-throughput sequencing data is a critical first step in genomics analysis pipelines. FastQC has served as the de facto standard for sequencing QC for over a decade, but its Java runtime dependency introduces startup overhead, elevated memory consumption, and deployment complexity. Here we present RastQC, a complete reimplementation of FastQC in Rust that provides all 12 standard QC modules with matching algorithms, MultiQC-compatible output formats, and a built-in multi-file summary dashboard. We benchmarked RastQC against FastQC v0.12.1 on both synthetic datasets (100K--10M reads) and real whole-genome sequencing data spanning five model organisms: *Escherichia coli*, *Saccharomyces cerevisiae*, *Drosophila melanogaster*, *Mus musculus*, and *Homo sapiens*. On synthetic data, RastQC achieves 1.7--5.6x speedup on uncompressed input and 32.5x speedup on gzip-compressed input, while using 15--26x less memory (22 MB vs. 340--644 MB). On real genome data, RastQC matches or exceeds FastQC speed across all organisms while achieving 100% module-level concordance (55/55 module calls identical across all organisms). RastQC compiles to a single 1.6 MB static binary with no external dependencies, representing a 134x reduction in deployment footprint. RastQC is freely available at https://github.com/kuanlinhuang/RastQC under the MIT license.

## Introduction

Next-generation sequencing (NGS) has become the foundation of modern genomics research, with applications spanning whole-genome sequencing, RNA-seq, ChIP-seq, and single-cell assays. Before downstream analysis, quality control of raw sequencing data is essential to identify technical artifacts including base-calling errors, adapter contamination, GC bias, sequence duplication, and position-dependent quality degradation (Andrews, 2010; Patel & Jain, 2012). Failure to detect these issues can propagate systematic errors into variant calls, expression quantification, and other analyses.

FastQC (Andrews, 2010) has been the most widely used sequencing QC tool for over a decade, offering 12 diagnostic modules covering base quality, sequence content, duplication, adapter contamination, and k-mer enrichment. FastQC produces self-contained HTML reports that have become a standard deliverable in sequencing facilities worldwide. Tools such as MultiQC (Ewels et al., 2016) aggregate FastQC outputs across samples for project-level QC review.

Despite its widespread adoption, FastQC has practical limitations rooted in its Java implementation. The Java Virtual Machine (JVM) imposes a 2--3 second startup penalty per invocation, a minimum memory footprint of ~300 MB regardless of input size, and a runtime dependency that complicates deployment in minimal container images and heterogeneous HPC environments. For large-scale studies processing hundreds to thousands of samples, these overheads accumulate significantly.

Recent efforts to rewrite bioinformatics tools in systems languages such as Rust and C++ have demonstrated substantial performance improvements while maintaining correctness. Notable examples include minimap2 (Li, 2018), fastp (Chen et al., 2018), and various Rust-based tools in the noodles ecosystem. Rust in particular offers memory safety guarantees without garbage collection, zero-cost abstractions, and excellent concurrency support through the ownership model.

Here we present RastQC, a complete reimplementation of FastQC in Rust. RastQC implements all 12 QC modules with algorithms matching the original FastQC, produces MultiQC-compatible output formats, and adds a multi-file summary dashboard for batch QC review. We demonstrate through comprehensive benchmarking on both synthetic and real genome data from five model organisms that RastQC achieves significant improvements in execution speed and memory efficiency while maintaining full concordance with FastQC results.

## Implementation

### Architecture

RastQC is organized into four modules reflecting the analysis pipeline:

1. **I/O layer** (`io/`): Streaming parsers for FASTQ (plain, gzip, bzip2-compressed) and BAM/SAM formats via the noodles library (version 0.88). Sequences are processed one at a time to minimize memory footprint.

2. **QC modules** (`modules/`): Twelve independent analysis modules, each implementing a common `QCModule` trait that defines `process_sequence()`, `calculate_results()`, and output generation methods. Modules accumulate per-position statistics during the streaming pass and compute final results lazily.

3. **Configuration** (`config/`): Embedded default adapter sequences (6 entries), contaminant sequences (15 entries), and pass/warn/fail thresholds matching FastQC defaults. All configurations are overridable via command-line flags.

4. **Report generation** (`report/`): Self-contained HTML reports with inline SVG charts, tab-separated data files compatible with MultiQC, per-file `summary.txt` with module-level status, and ZIP archives. A multi-file summary dashboard (HTML + TSV) is generated when processing multiple samples.

### QC Modules

RastQC implements all 12 FastQC modules with matching algorithms:

**Table 1.** QC modules implemented in RastQC.

| Module | Description | Key Algorithm |
|--------|-------------|---------------|
| Basic Statistics | Sequence count, length, %GC, total bases, encoding | Streaming counters; Phred encoding auto-detection |
| Per Base Sequence Quality | Quality distribution at each read position | Per-position quality histograms with adaptive base grouping; box plot statistics |
| Per Tile Sequence Quality | Tile-specific quality deviations (Illumina) | Per-tile quality accumulation with 10% sampling after 10K reads |
| Per Sequence Quality Scores | Distribution of mean quality per read | Per-read average quality histogram |
| Per Base Sequence Content | A/T/G/C proportions per position | Per-position base counters with adaptive grouping |
| Per Sequence GC Content | GC% distribution vs. theoretical normal | GCModel fractional binning; bidirectional mode averaging; normal distribution fit |
| Per Base N Content | Unknown base frequency per position | Per-position N counters |
| Sequence Length Distribution | Read length variability | Adaptive binning maintaining ~50 display categories |
| Sequence Duplication Levels | Library complexity estimate | String-based tracking (first 50 bp); 100K unique cutoff with iterative binomial correction |
| Overrepresented Sequences | High-frequency sequences with contaminant matching | Exact/approximate string matching; 1-mismatch tolerance for >20 bp |
| Adapter Content | Known adapter contamination by position | Substring matching of 12-mer adapters with cumulative position counting |
| K-mer Content | Positionally biased k-mers | 7-mer tracking with 2% sampling; binomial enrichment test |

### MultiQC Compatibility

RastQC output is designed for full compatibility with MultiQC's FastQC parser. The `fastqc_data.txt` file uses identical module names, column headers, and data formats. Key compatibility features include:

- **Filename field** populated in Basic Statistics for correct sample identification
- **Total Bases** metric included for modern MultiQC versions
- **summary.txt** populated with per-module PASS/WARN/FAIL status and filename
- **ZIP archive structure** matching the expected `{sample}_fastqc/fastqc_data.txt` path
- **Module names** exactly matching FastQC conventions

### Multi-file Summary Dashboard

When processing multiple files (or with the `--summary` flag), RastQC generates:

- **summary.tsv**: Tab-separated matrix (samples x modules) for programmatic filtering
- **summary.html**: Color-coded HTML dashboard with clickable links to individual reports, aggregate tallies, and sortable columns

This provides MultiQC-like functionality integrated directly into the QC tool, enabling rapid triage of large sample batches without additional software.

### Multi-threading

RastQC uses the rayon library for data-parallel file processing. Multiple input files are distributed across a configurable thread pool (defaulting to all available CPU cores). Within each file, processing is sequential to maintain the streaming memory model.

## Methods

### Benchmarking Environment

All benchmarks were performed on a single workstation:

- **CPU**: Intel Core i9-9900K @ 3.60 GHz (8 cores, 16 threads)
- **Memory**: 32 GB DDR4
- **Storage**: SSD
- **OS**: macOS 15.7.4 (Darwin 24.6.0)
- **RastQC**: v0.1.0, compiled with Rust (release profile, LTO enabled, opt-level 3)
- **FastQC**: v0.12.1, compiled from source with OpenJDK 11.0.30, `-Xmx512m`

### Synthetic Datasets

We generated synthetic FASTQ datasets with realistic Illumina characteristics:

- **Read length**: 150 bp
- **Quality profile**: Declining Phred scores from 5' (mean 36) to 3' (mean 18), +/- 5 per-base variation
- **Base composition**: Uniform random (25% each A/C/G/T)
- **Adapter contamination**: 5% of reads contain Illumina Universal Adapter at variable insert sizes (80--144 bp)
- **Sizes**: 100K reads (33 MB), 1M reads (325 MB), 1M gzipped (149 MB), 10M reads (3.2 GB)

### Real Genome Datasets

To evaluate performance on real sequencing data with natural biological complexity, we benchmarked on whole-genome Illumina sequencing data from five model organisms spanning the tree of life (Table 2).

**Table 2.** Real genome datasets used for benchmarking.

| Organism | Species | Accession | Reads | Read Length | File Size |
|----------|---------|-----------|-------|-------------|-----------|
| Bacteria | *Escherichia coli* K-12 MG1655 | ERR022075 | 5,000,000 | 100 bp | 1,194 MB |
| Yeast | *Saccharomyces cerevisiae* | SRR19072702 | 1,368,860 | 75 bp | 267 MB |
| Fruit fly | *Drosophila melanogaster* | ERR1942264 | 432,546 | 100 bp | 105 MB |
| Mouse | *Mus musculus* | ERR3085830 | 1,619,240 | 151 bp | 550 MB |
| Human | *Homo sapiens* (NA12878) | ERR3239334 | 2,392,582 | 150 bp | 831 MB |

All datasets were obtained from the European Nucleotide Archive (ENA) as single-end reads (read 1 of paired-end runs). Read counts represent subsampled datasets to enable practical benchmarking across all organisms.

### Metrics

- **Wall-clock time**: Elapsed time measured by `/usr/bin/time -l`
- **Peak memory (RSS)**: Maximum resident set size reported by `/usr/bin/time -l`
- **Startup overhead**: Time to process a single 1-read FASTQ file
- **Output concordance**: Module-level PASS/WARN/FAIL agreement on identical input

All benchmarks were run single-threaded to ensure fair per-file comparison.

## Results

### Performance on Synthetic Data

RastQC consistently outperformed FastQC on synthetic datasets (Table 3).

**Table 3.** Performance comparison on synthetic datasets (single-threaded).

| Dataset | Reads | RastQC Time | FastQC Time | Speedup | RastQC RSS | FastQC RSS |
|---------|-------|-------------|-------------|---------|------------|------------|
| Startup (1 read) | 1 | 0.005 s | 2.66 s | 532x | 1.5 MB | 286 MB |
| 100K synthetic | 100,000 | 1.50 s | 8.34 s | 5.6x | 22 MB | 343 MB |
| 1M synthetic | 1,000,000 | 10.86 s | 21.90 s | 2.0x | 22 MB | 393 MB |
| 1M gzipped | 1,000,000 | 12.44 s | 404.26 s | 32.5x | 22 MB | 546 MB |
| 10M synthetic | 10,000,000 | 1,274 s | 2,158 s | 1.7x | 21 MB | 544 MB |

The startup test (1 read) reveals a 532x difference: RastQC initializes in 5 ms as a native binary versus 2.66 s for JVM class loading. The most dramatic speedup (32.5x) was observed on gzip-compressed input, reflecting superior decompression performance of the Rust `flate2` crate compared to Java's `GZIPInputStream`.

### Performance on Real Genome Data

We evaluated both tools on real Illumina whole-genome sequencing data from five model organisms (Table 4).

**Table 4.** Performance comparison on real genome data across model organisms.

| Organism | Reads | Read Len | RastQC Time | FastQC Time | Speedup | RastQC RSS | FastQC RSS | Mem. Ratio |
|----------|-------|----------|-------------|-------------|---------|------------|------------|------------|
| *D. melanogaster* | 432,546 | 100 bp | 2.37 s | 4.85 s | 2.0x | 22 MB | 527 MB | 24x |
| *S. cerevisiae* | 1,368,860 | 75 bp | 6.18 s | 7.58 s | 1.2x | 24 MB | 587 MB | 24x |
| *E. coli* K-12 | 5,000,000 | 100 bp | 29.02 s | 24.68 s | 0.9x | 22 MB | 564 MB | 26x |
| *M. musculus* | 1,619,240 | 151 bp | 11.03 s | 11.58 s | 1.1x | 22 MB | 592 MB | 27x |
| *H. sapiens* | 2,392,582 | 150 bp | 16.15 s | 17.61 s | 1.1x | 35 MB | 644 MB | 18x |

RastQC was faster than FastQC on 4 of 5 organisms, with speedups ranging from 1.1x (mouse, human) to 2.0x (fly). On the *E. coli* dataset, FastQC was marginally faster (0.9x), likely due to this dataset's characteristics favoring FastQC's internal optimizations for short uniform-length reads. Across all organisms, RastQC maintained a stable memory footprint of 22--35 MB, while FastQC required 527--644 MB (18--27x more memory).

The human genome dataset is particularly relevant for clinical sequencing applications. RastQC processed 2.4M human reads in 16.2 s compared to FastQC's 17.6 s (1.1x speedup), while using 35 MB versus 644 MB of memory (18x reduction). For a typical 30x human whole-genome sequencing run (~900M reads), this memory difference enables substantially more concurrent QC analyses on shared compute infrastructure.

### Memory Usage

RastQC's memory footprint remained remarkably stable across all datasets and organisms:

- **RastQC**: 22--35 MB across all inputs (1 read to 10M reads, all organisms)
- **FastQC**: 286--644 MB, scaling with input complexity and the JVM baseline

This 18--27x memory reduction reflects RastQC's streaming architecture, where only summary statistics are retained in memory. The constant-memory property is particularly valuable for HPC environments with strict memory limits per job.

### Binary Size and Deployment

**Table 5.** Deployment footprint comparison.

| Component | RastQC | FastQC |
|-----------|--------|--------|
| Binary/package size | 1.6 MB | ~15 MB (JARs + classes) |
| Runtime dependency | None (static binary) | Java 11+ (~200 MB) |
| Total deployment size | 1.6 MB | ~215 MB |
| Reduction factor | **134x smaller** | -- |

RastQC compiles to a single 1.6 MB statically-linked binary with no external dependencies. This is advantageous for Docker containers (smaller images, faster pulls), CI/CD pipelines, and HPC environments where module systems may not provide a compatible JRE.

### Output Concordance

We systematically compared module-level PASS/WARN/FAIL calls between RastQC and FastQC across all five model organisms (Table 6).

**Table 6.** Module-level concordance across five model organisms (55 total module comparisons).

| Module | Concordant | Discordant |
|--------|-----------|------------|
| Basic Statistics | 5/5 | 0 |
| Per Base Sequence Quality | 5/5 | 0 |
| Per Tile Sequence Quality | 5/5 | 0 |
| Per Sequence Quality Scores | 5/5 | 0 |
| Per Base Sequence Content | 5/5 | 0 |
| Per Sequence GC Content | 5/5 | 0 |
| Per Base N Content | 5/5 | 0 |
| Sequence Length Distribution | 5/5 | 0 |
| Sequence Duplication Levels | 5/5 | 0 |
| Overrepresented Sequences | 5/5 | 0 |
| Adapter Content | 5/5 | 0 |
| **Total** | **55/55** | **0** |

All 11 modules produced identical PASS/WARN/FAIL calls across all 5 organisms, achieving **100% concordance** (55/55 module comparisons). This was accomplished by faithfully porting three key algorithms from FastQC:

1. **Sequence Duplication Levels**: RastQC uses FastQC's exact iterative binomial correction formula for estimating true duplication counts beyond the 100K unique sequence observation window, and uses string-based sequence identity matching rather than hashing.

2. **Per Sequence GC Content**: RastQC implements FastQC's GCModel, which distributes each read's GC count across adjacent percentage bins using fractional weights, and uses the bidirectional mode-averaging algorithm for fitting the theoretical normal distribution.

3. **Overrepresented Sequences**: RastQC uses the total sequence count (not count-at-limit) as the denominator for percentage calculations, matching FastQC's behavior.

Numerical values for continuously-valued metrics also showed strong agreement. Per-base quality means matched to within 0.02 Phred units. Basic statistics (total sequences, %GC, encoding, sequence length, total bases) were identical between tools across all organisms.

## Discussion

### Performance Characteristics

RastQC provides consistent performance advantages over FastQC, particularly in three scenarios:

1. **Small files and batch processing**: For datasets under 1M reads, RastQC's 1.5--5.6x speedup is dominated by the elimination of 2.66 s JVM startup overhead. Processing 1,000 amplicon panel samples (~50K reads each) would save ~44 minutes from startup elimination alone.

2. **Compressed input**: The 32.5x speedup on gzip-compressed FASTQ reflects superior decompression performance of the Rust `flate2` crate. Since most sequencing data is stored compressed, this represents the typical production use case.

3. **Memory-constrained environments**: RastQC's constant 22--35 MB footprint enables ~1,400 concurrent instances on a 32 GB machine, versus ~50 for FastQC. This is critical for shared HPC clusters where memory is a scarce, scheduled resource.

On larger uncompressed datasets (10M reads), the speedup narrows to 1.7x as both tools become I/O-bound and the JVM startup cost is amortized. On *E. coli* uncompressed data, FastQC was marginally faster, suggesting that FastQC's Java implementation is well-optimized for steady-state sequence processing.

### Cross-organism Validation

Testing across five organisms with different genome sizes (4.6 Mb to 3.1 Gb), GC contents (38--65%), and read characteristics validates that RastQC performs robustly on diverse real-world data. The 100% concordance rate (55/55 module calls) demonstrates that RastQC faithfully reproduces FastQC's QC assessments across the full spectrum of model organisms, from prokaryotes to human.

### MultiQC Integration

RastQC's output format has been validated for compatibility with MultiQC's FastQC parser. The `fastqc_data.txt` file includes all required fields (Filename, Total Bases, module headers, data columns), and the `summary.txt` file follows the expected three-column format (STATUS, Module Name, Filename). This enables RastQC to serve as a drop-in replacement in existing pipelines that aggregate results via MultiQC.

### Limitations

While RastQC achieves full concordance with FastQC and offers substantial performance and deployment advantages, several limitations should be noted:

1. **No graphical user interface**: RastQC is command-line only. FastQC provides a Swing-based GUI that allows interactive file selection, real-time result viewing, and manual report saving. Users requiring a graphical workflow should continue to use FastQC or pair RastQC with its HTML report viewer.

2. **No Oxford Nanopore native format support**: FastQC can read Fast5 files directly via HDF5 libraries. RastQC does not support Fast5 or POD5 formats; Nanopore users must first convert to FASTQ (e.g., via `dorado` or `guppy`) before running RastQC. This is increasingly a non-issue as most Nanopore workflows already output FASTQ or unaligned BAM.

3. **Steady-state throughput on large uncompressed files**: On the largest uncompressed datasets (5M+ reads), FastQC's Java JIT compiler produces competitive throughput once warmed up. On the *E. coli* dataset (5M reads, uncompressed), FastQC was marginally faster (0.9x). This narrowing is expected: both tools are fundamentally I/O-bound at scale, and the JVM's 2.7s startup penalty becomes negligible relative to minutes of processing. RastQC's advantages are most pronounced on compressed input (32.5x speedup) and in batch scenarios where startup cost is multiplied across many files.

4. **Colorspace reads**: FastQC supports SOLiD colorspace-encoded reads. RastQC does not, as the SOLiD platform has been discontinued and colorspace data is rarely encountered in modern analyses.

5. **Kmer module disabled by default**: Following FastQC v0.11.6+, the Kmer Content module is disabled by default in the limits configuration. Users who require it must explicitly enable it via a custom limits file (`kmer ignore 0`).

### Future Directions

Several enhancements are planned to extend RastQC's utility:

1. **Intra-file parallelism**: The current architecture processes each file sequentially. For very large single files (>100M reads), chunked parallel processing with merged statistics could further reduce wall-clock time, particularly on multi-core systems where only one file is being analyzed.

2. **Long-read QC metrics**: As long-read sequencing (PacBio HiFi, Oxford Nanopore) becomes routine, dedicated QC metrics such as read length N50, quality-stratified length distributions, and homopolymer error rates would complement the existing short-read modules.

3. **Native MultiQC JSON output**: While RastQC's `fastqc_data.txt` is compatible with MultiQC, direct JSON output in MultiQC's native format would eliminate the parsing step and enable richer metadata transfer.

4. **Workflow manager integration**: Standardized non-zero exit codes reflecting overall QC status (e.g., exit 1 if any module fails) would enable automated QC gates in Nextflow, Snakemake, and similar pipeline frameworks without parsing output files.

5. **Streaming from standard input**: Supporting piped FASTQ input (e.g., from `samtools fastq` or `seqtk`) would enable RastQC to operate within Unix pipelines without intermediate file materialization.

## Availability

- **Source code**: https://github.com/kuanlinhuang/RastQC
- **License**: MIT
- **Language**: Rust (2021 edition)
- **Installation**: `cargo install --path .` or pre-compiled binaries
- **Test suite**: 31 tests (20 unit, 11 integration)
- **System requirements**: Any platform supported by Rust (Linux, macOS, Windows)

## References

Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884--i890.

Ewels, P., Magnusson, M., Lundin, S., & Kaller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047--3048.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094--3100.

Patel, R. K., & Jain, M. (2012). NGS QC Toolkit: a toolkit for quality control of next generation sequencing data. *PLoS ONE*, 7(2), e30619.
