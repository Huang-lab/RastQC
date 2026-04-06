# RastQC Agent Skill

RastQC is a fast sequencing quality control tool. Use this skill to run QC on FASTQ, BAM, or SAM files automatically.

## Installation

```bash
# Install from GitHub (requires Rust 1.70+)
cargo install --git https://github.com/Huang-lab/RastQC.git

# With Nanopore Fast5/POD5 support
cargo install --git https://github.com/Huang-lab/RastQC.git --features nanopore
```

After installation, `rastqc` is available in `$PATH`. Verify with `rastqc --version`.

If `cargo` is not in PATH, check `~/.cargo/bin/cargo`. If Rust is not installed: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`.

## Quick Reference

### Basic QC (short-read Illumina data)
```bash
rastqc -o <output_dir> <input_files...>
```

### Long-read QC (ONT / PacBio)
```bash
rastqc --long-read -o <output_dir> <input_files...>
```
The `--long-read` flag is auto-enabled for `.fast5` and `.pod5` files. For long-read FASTQ, you must pass it explicitly.

### Common Options

| Flag | Purpose |
|------|---------|
| `-o <DIR>` | Output directory (default: current dir) |
| `-t <N>` | Number of threads (default: all CPUs) |
| `--long-read` | Enable long-read modules (N50, Quality Stratified Length, Homopolymer) |
| `--time` | Show per-file timing breakdown |
| `--nozip` | Write HTML report only, skip ZIP |
| `--summary` | Force multi-file summary (auto for 2+ files) |
| `--multiqc-json` | Output native MultiQC JSON |
| `--exit-code` | Return 0=pass, 1=warn, 2=fail for pipeline gates |
| `--serve` | Launch web browser to view reports |
| `-q` / `--quiet` | Suppress progress output |
| `--stdin` | Read FASTQ from stdin pipe |

### Input Formats
- FASTQ (plain, `.gz`, `.bz2`)
- BAM / SAM
- Fast5 / POD5 (requires `--features nanopore` build)
- Standard input (`--stdin` or pass `-` as filename)

## Output Files

For each input file `sample.fastq.gz`:

| File | Description |
|------|-------------|
| `sample_fastqc.html` | Self-contained HTML report with SVG charts |
| `sample_fastqc.zip` | ZIP archive (contains HTML + data + summary) |
| `sample_multiqc.json` | MultiQC JSON (only with `--multiqc-json`) |

For multiple files:

| File | Description |
|------|-------------|
| `summary.html` | Cross-sample dashboard with pass/warn/fail matrix |
| `summary.tsv` | Machine-readable summary for scripting |

## Agent Workflow: Automatic QC

When asked to run QC on sequencing data, follow this workflow:

### Step 1: Check rastqc is installed
```bash
which rastqc || cargo install --git https://github.com/Huang-lab/RastQC.git
```

### Step 2: Detect input files
```bash
ls <data_dir>/*.fastq.gz <data_dir>/*.fq.gz <data_dir>/*.bam 2>/dev/null
```

### Step 3: Determine read type
Check if files are long-read by examining a few reads:
```bash
# Peek at read lengths (if mean > 500bp, it's long-read)
gzip -dc <file> 2>/dev/null | head -400 | awk 'NR%4==2 {sum+=length($0); n++} END {printf "%.0f\n", sum/n}'
```
If mean read length > 500 bp, or filename contains `ont`, `nanopore`, `pacbio`, `hifi`, `revio`, or extension is `.fast5`/`.pod5`, add `--long-read`.

### Step 4: Run RastQC
```bash
# Short-read
rastqc --time -o <output_dir> <files...>

# Long-read
rastqc --long-read --time -o <output_dir> <files...>

# With MultiQC JSON
rastqc --multiqc-json --time -o <output_dir> <files...>
```

### Step 5: Interpret results
Read the stderr output for the summary, or parse:
```bash
# Quick pass/fail check
cat <output_dir>/summary.tsv

# Check specific file's module results
unzip -p <output_dir>/sample_fastqc.zip "*/summary.txt"
```

### Step 6: Report to user
Summarize:
- Total files processed and wall-clock time
- Number of pass/warn/fail modules
- Any files with failures and which modules failed
- Open the HTML report: `open <output_dir>/summary.html`

## Interpreting Module Results

| Module | WARN means | FAIL means |
|--------|-----------|------------|
| Per Base Quality | Median < Q25 at any position | Median < Q20 at any position |
| Per Sequence Quality | Mode quality <= Q27 | Mode quality <= Q20 |
| Per Base Content | \|A-T\| or \|G-C\| > 10% | > 20% |
| Per Sequence GC | >15% deviation from normal | >30% deviation |
| N Content | N% > 5% at any position | N% > 20% |
| Duplication | <70% unique reads | <50% unique |
| Overrepresented | Any sequence > 0.1% | > 1% |
| Adapter Content | >5% adapter at any position | >10% |
| Homopolymer (LR) | >5% bases in runs | >10% bases in runs |

## Pipeline Integration Example

```bash
# Nextflow-style QC gate
rastqc --exit-code -o qc_results/ sample.fastq.gz
if [ $? -eq 2 ]; then
    echo "QC FAILED - stopping pipeline"
    exit 1
fi

# Pipe from samtools
samtools fastq aligned.bam | rastqc --stdin -o qc_results/

# Batch + MultiQC
rastqc --multiqc-json -o qc_results/ data/*.fastq.gz
multiqc qc_results/
```
