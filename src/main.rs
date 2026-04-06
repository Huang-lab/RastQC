mod config;
mod gui;
mod io;
mod modules;
mod parallel;
mod report;

use anyhow::Result;
use clap::Parser;
use rayon::prelude::*;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use config::FastQCConfig;
use io::SequenceReader;
use modules::{ModuleFactory, QCResult};
use report::{
    generate_html_report, generate_multiqc_json, generate_summary_html, generate_summary_tsv,
    generate_summary_txt, generate_text_data, write_zip_archive,
};

/// Per-file result summary used for the multi-file overview.
pub struct FileSummary {
    pub filename: String,
    pub module_results: Vec<(String, QCResult)>,
    pub total_sequences: u64,
    pub report_path: String,
    pub timing: Option<FileTiming>,
}

/// Per-file timing breakdown.
pub struct FileTiming {
    pub read_and_qc_secs: f64,
    pub report_gen_secs: f64,
    pub io_write_secs: f64,
    pub total_secs: f64,
}

#[derive(Parser, Debug)]
#[command(name = "rastqc", version = "0.1.0")]
#[command(about = "RastQC - A quality control tool for high throughput sequence data")]
struct Cli {
    /// Input files (FASTQ, BAM, SAM). Use "-" to read FASTQ from stdin.
    files: Vec<PathBuf>,

    /// Read FASTQ from stdin (equivalent to passing "-" as input)
    #[arg(long)]
    stdin: bool,

    /// Output directory
    #[arg(short, long)]
    outdir: Option<PathBuf>,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = num_cpus())]
    threads: usize,

    /// Custom contaminant file
    #[arg(short = 'c', long)]
    contaminants: Option<PathBuf>,

    /// Custom adapter file
    #[arg(short = 'a', long)]
    adapters: Option<PathBuf>,

    /// Custom limits file
    #[arg(short = 'l', long)]
    limits: Option<PathBuf>,

    /// Kmer size (default 7)
    #[arg(short = 'k', long, default_value_t = 7)]
    kmer_size: usize,

    /// Don't filter low quality reads
    #[arg(long)]
    nofilter: bool,

    /// Extract ZIP archive after creation
    #[arg(long)]
    extract: bool,

    /// Don't create ZIP, only HTML
    #[arg(long)]
    nozip: bool,

    /// Write multi-file summary report (summary.html + summary.tsv)
    #[arg(long)]
    summary: bool,

    /// Quiet mode - suppress progress
    #[arg(short = 'q', long)]
    quiet: bool,

    /// Length to truncate sequences for duplication detection
    #[arg(long, default_value_t = 50)]
    dup_length: usize,

    /// Output native MultiQC JSON (multiqc_fastqc.json) alongside standard reports
    #[arg(long)]
    multiqc_json: bool,

    /// Return QC-aware exit codes: 0=all pass, 1=warnings, 2=failures.
    /// Useful for automated QC gates in Nextflow, Snakemake, etc.
    #[arg(long)]
    exit_code: bool,

    /// Start a local web server to browse reports (default port: 8080)
    #[arg(long)]
    serve: bool,

    /// Port for the web server (used with --serve)
    #[arg(long, default_value_t = 8080)]
    port: u16,

    /// Disable intra-file parallelism (streaming parallel is on by default for files >50MB).
    #[arg(long)]
    no_parallel: bool,

    /// Show per-file and per-step timing breakdown
    #[arg(long)]
    time: bool,

    /// Enable long-read QC modules (Read Length N50, Quality Stratified Length,
    /// Homopolymer Content). Auto-enabled for Fast5/POD5 inputs.
    #[arg(long)]
    long_read: bool,
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => code,
        Err(e) => {
            eprintln!("Error: {:#}", e);
            ExitCode::from(3)
        }
    }
}

fn run() -> Result<ExitCode> {
    let mut cli = Cli::parse();

    // Handle --stdin flag: add "-" to file list
    if cli.stdin {
        cli.files.push(PathBuf::from("-"));
    }

    if cli.files.is_empty() {
        anyhow::bail!("No input files specified. Pass file paths or use --stdin to read from standard input.");
    }

    // Auto-detect long-read mode from file extensions
    let long_read = cli.long_read
        || cli.files.iter().any(|f| {
            let name = f.to_string_lossy().to_lowercase();
            name.ends_with(".fast5") || name.ends_with(".pod5")
        });

    let mut config = FastQCConfig::new(
        cli.contaminants.as_deref(),
        cli.adapters.as_deref(),
        cli.limits.as_deref(),
        cli.kmer_size,
        cli.nofilter,
        cli.dup_length,
        long_read,
    )?;

    if long_read {
        config.enable_long_read_modules();
        if !cli.quiet {
            eprintln!("Long-read mode enabled: running Read Length N50, Quality Stratified Length, Homopolymer Content");
        }
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .ok();

    let outdir = cli.outdir.clone().unwrap_or_else(|| PathBuf::from("."));
    if !outdir.exists() {
        std::fs::create_dir_all(&outdir)?;
    }

    let total_start = Instant::now();
    let total_files = cli.files.len();

    // Shared counter for progress
    let completed = AtomicUsize::new(0);
    let summaries: Mutex<Vec<FileSummary>> = Mutex::new(Vec::new());

    cli.files.par_iter().for_each(|file| {
        match process_file(file, &outdir, &config, &cli) {
            Ok(summary) => {
                let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                if !cli.quiet {
                    let status = summary_status_line(&summary);
                    if let Some(ref t) = summary.timing {
                        eprintln!(
                            "[{}/{}] {} ({} seqs) {}  {:.2}s  [qc {:.2}s | report {:.3}s | io {:.3}s]",
                            done,
                            total_files,
                            summary.filename,
                            format_count(summary.total_sequences),
                            status,
                            t.total_secs,
                            t.read_and_qc_secs,
                            t.report_gen_secs,
                            t.io_write_secs,
                        );
                    } else {
                        eprintln!(
                            "[{}/{}] {} ({} seqs) {}  {:.1}s",
                            done,
                            total_files,
                            summary.filename,
                            format_count(summary.total_sequences),
                            status,
                            total_start.elapsed().as_secs_f64()
                        );
                    }
                }
                summaries.lock().unwrap().push(summary);
            }
            Err(e) => {
                completed.fetch_add(1, Ordering::Relaxed);
                eprintln!("Error processing {}: {}", file.display(), e);
            }
        }
    });

    let mut summaries = summaries.into_inner().unwrap();
    // Sort summaries by filename for consistent output
    summaries.sort_by(|a, b| a.filename.cmp(&b.filename));

    // Print final tally
    if !cli.quiet {
        print_final_tally(&summaries, total_files, total_start.elapsed().as_secs_f64(), cli.time);
    }

    // Write multi-file summary if requested (or auto-enable for 2+ files)
    if cli.summary || total_files > 1 {
        let tsv = generate_summary_tsv(&summaries);
        let tsv_path = outdir.join("summary.tsv");
        std::fs::write(&tsv_path, &tsv)?;

        let html = generate_summary_html(&summaries);
        let html_path = outdir.join("summary.html");
        std::fs::write(&html_path, &html)?;

        if !cli.quiet {
            eprintln!("Summary: {} and {}", tsv_path.display(), html_path.display());
        }
    }

    // Start web server if requested
    if cli.serve {
        gui::start_server(&outdir, cli.port)?;
        return Ok(ExitCode::SUCCESS);
    }

    // Determine exit code based on QC results
    if cli.exit_code {
        let has_fail = summaries
            .iter()
            .any(|s| s.module_results.iter().any(|(_, r)| *r == QCResult::Fail));
        let has_warn = summaries
            .iter()
            .any(|s| s.module_results.iter().any(|(_, r)| *r == QCResult::Warn));
        if has_fail {
            return Ok(ExitCode::from(2));
        }
        if has_warn {
            return Ok(ExitCode::from(1));
        }
    }

    Ok(ExitCode::SUCCESS)
}

fn process_file(
    file: &PathBuf,
    outdir: &PathBuf,
    config: &FastQCConfig,
    cli: &Cli,
) -> Result<FileSummary> {
    let file_start = Instant::now();
    let is_stdin = file.as_os_str() == "-";

    // Choose between streaming parallel and sequential processing
    let qc_start = Instant::now();
    let (qc_modules, count) = if !is_stdin && !cli.no_parallel && parallel::should_use_parallel(file) {
        parallel::process_file_parallel(file, config, cli.threads)?
    } else {
        let mut reader = if is_stdin {
            SequenceReader::from_stdin()
        } else {
            SequenceReader::open(file)?
        };
        let mut qc_modules = ModuleFactory::create_modules(config);

        let mut count: u64 = 0;
        while let Some(seq) = reader.next_sequence()? {
            for module in &mut qc_modules {
                module.process_sequence(&seq);
            }
            count += 1;
        }

        for module in &mut qc_modules {
            module.calculate_results(config);
        }
        (qc_modules, count)
    };
    let read_and_qc_secs = qc_start.elapsed().as_secs_f64();

    // Determine output base name
    let filename = if is_stdin {
        "stdin.fastq".to_string()
    } else {
        file.file_name().unwrap().to_string_lossy().to_string()
    };
    let stem = filename
        .trim_end_matches(".gz")
        .trim_end_matches(".bz2")
        .trim_end_matches(".fastq")
        .trim_end_matches(".fq")
        .trim_end_matches(".bam")
        .trim_end_matches(".sam")
        .trim_end_matches(".fast5")
        .trim_end_matches(".pod5")
        .to_string();

    // Collect module results for summary
    let module_results: Vec<(String, QCResult)> = qc_modules
        .iter()
        .map(|m| (m.name().to_string(), m.result()))
        .collect();

    // Generate outputs
    let report_start = Instant::now();
    let html = generate_html_report(&filename, &qc_modules);
    let text = generate_text_data(&filename, &qc_modules);
    let summary_txt = generate_summary_txt(&filename, &qc_modules);

    // Generate MultiQC JSON if requested
    let multiqc_json_str = if cli.multiqc_json {
        Some(generate_multiqc_json(&filename, &qc_modules))
    } else {
        None
    };
    let report_gen_secs = report_start.elapsed().as_secs_f64();

    // Write outputs
    let io_start = Instant::now();

    if let Some(json) = multiqc_json_str {
        let json_path = outdir.join(format!("{}_multiqc.json", stem));
        std::fs::write(&json_path, &json)?;
    }

    let report_path;
    if cli.nozip {
        report_path = format!("{}_fastqc.html", stem);
        let html_path = outdir.join(&report_path);
        std::fs::write(&html_path, &html)?;
    } else {
        let zip_path = outdir.join(format!("{}_fastqc.zip", stem));
        write_zip_archive(&zip_path, &stem, &html, &text, &summary_txt)?;
        report_path = format!("{}_fastqc.html", stem);

        if cli.extract {
            let extract_dir = outdir.join(format!("{}_fastqc", stem));
            std::fs::create_dir_all(&extract_dir)?;
            std::fs::write(extract_dir.join("fastqc_report.html"), &html)?;
            std::fs::write(extract_dir.join("fastqc_data.txt"), &text)?;
            std::fs::write(extract_dir.join("summary.txt"), &summary_txt)?;
        }

        // Always write standalone HTML alongside ZIP for summary linking
        let html_path = outdir.join(&report_path);
        std::fs::write(&html_path, &html)?;
    }
    let io_write_secs = io_start.elapsed().as_secs_f64();
    let total_secs = file_start.elapsed().as_secs_f64();

    let timing = if cli.time {
        Some(FileTiming {
            read_and_qc_secs,
            report_gen_secs,
            io_write_secs,
            total_secs,
        })
    } else {
        None
    };

    Ok(FileSummary {
        filename,
        module_results,
        total_sequences: count,
        report_path,
        timing,
    })
}

fn summary_status_line(summary: &FileSummary) -> String {
    let mut pass = 0;
    let mut warn = 0;
    let mut fail = 0;
    for (_, result) in &summary.module_results {
        match result {
            QCResult::Pass => pass += 1,
            QCResult::Warn => warn += 1,
            QCResult::Fail => fail += 1,
            QCResult::NotRun => {}
        }
    }
    let mut parts = Vec::new();
    if pass > 0 {
        parts.push(format!("{} pass", pass));
    }
    if warn > 0 {
        parts.push(format!("{} warn", warn));
    }
    if fail > 0 {
        parts.push(format!("{} FAIL", fail));
    }
    parts.join(", ")
}

fn format_count(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}K", n as f64 / 1_000.0)
    } else {
        format!("{}", n)
    }
}

fn print_final_tally(summaries: &[FileSummary], total_files: usize, elapsed: f64, show_time: bool) {
    eprintln!();
    eprintln!("=== Analysis complete ===");
    eprintln!(
        "{} file(s) processed in {:.1}s",
        total_files, elapsed
    );

    let mut total_pass = 0usize;
    let mut total_warn = 0usize;
    let mut total_fail = 0usize;
    let mut files_with_fail = 0usize;
    let mut files_with_warn = 0usize;

    for summary in summaries {
        let mut file_has_fail = false;
        let mut file_has_warn = false;
        for (_, result) in &summary.module_results {
            match result {
                QCResult::Pass => total_pass += 1,
                QCResult::Warn => {
                    total_warn += 1;
                    file_has_warn = true;
                }
                QCResult::Fail => {
                    total_fail += 1;
                    file_has_fail = true;
                }
                QCResult::NotRun => {}
            }
        }
        if file_has_fail {
            files_with_fail += 1;
        } else if file_has_warn {
            files_with_warn += 1;
        }
    }

    eprintln!(
        "Module checks: {} pass, {} warn, {} fail",
        total_pass, total_warn, total_fail
    );

    if files_with_fail > 0 {
        eprintln!(
            "Files with failures: {}/{}",
            files_with_fail, total_files
        );
    }
    if files_with_warn > 0 {
        eprintln!(
            "Files with warnings: {}/{}",
            files_with_warn, total_files
        );
    }

    // List files with failures
    if files_with_fail > 0 {
        eprintln!();
        eprintln!("Failed modules by file:");
        for summary in summaries {
            let fails: Vec<&str> = summary
                .module_results
                .iter()
                .filter(|(_, r)| *r == QCResult::Fail)
                .map(|(name, _)| name.as_str())
                .collect();
            if !fails.is_empty() {
                eprintln!("  {} -- {}", summary.filename, fails.join(", "));
            }
        }
    }

    // Timing breakdown
    if show_time {
        eprintln!();
        eprintln!("--- Timing breakdown ---");
        eprintln!("{:<35} {:>8} {:>8} {:>8} {:>8}", "File", "QC", "Report", "IO", "Total");
        let mut sum_qc = 0.0_f64;
        let mut sum_report = 0.0_f64;
        let mut sum_io = 0.0_f64;
        let mut sum_total = 0.0_f64;
        for summary in summaries {
            if let Some(ref t) = summary.timing {
                eprintln!(
                    "{:<35} {:>7.2}s {:>7.3}s {:>7.3}s {:>7.2}s",
                    truncate_str(&summary.filename, 35),
                    t.read_and_qc_secs,
                    t.report_gen_secs,
                    t.io_write_secs,
                    t.total_secs,
                );
                sum_qc += t.read_and_qc_secs;
                sum_report += t.report_gen_secs;
                sum_io += t.io_write_secs;
                sum_total += t.total_secs;
            }
        }
        if summaries.len() > 1 {
            eprintln!("{:<35} {:>7.2}s {:>7.3}s {:>7.3}s {:>7.2}s", "SUM (cpu-time)", sum_qc, sum_report, sum_io, sum_total);
            eprintln!("Wall-clock: {:.2}s", elapsed);
        }
    }

    eprintln!();
}

fn truncate_str(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("...{}", &s[s.len() - (max_len - 3)..])
    }
}
