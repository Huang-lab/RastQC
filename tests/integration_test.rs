use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn binary_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_BIN_EXE_rastqc"));
    // Fall back to target/debug if env not available
    if !path.exists() {
        path = PathBuf::from("target/debug/rastqc");
    }
    path
}

fn test_dir() -> PathBuf {
    let dir = PathBuf::from("tests/tmp");
    let _ = fs::create_dir_all(&dir);
    dir
}

fn write_fastq(path: &PathBuf, records: &[(&str, &str, &str)]) {
    let mut content = String::new();
    for (name, seq, qual) in records {
        content.push_str(&format!("@{}\n{}\n+\n{}\n", name, seq, qual));
    }
    fs::write(path, content).unwrap();
}

fn cleanup(dir: &PathBuf) {
    let _ = fs::remove_dir_all(dir);
}

// ---------------------------------------------------------------------------
// Integration tests
// ---------------------------------------------------------------------------

#[test]
fn test_single_file_nozip() {
    let outdir = test_dir().join("single_nozip");
    let _ = fs::create_dir_all(&outdir);
    let input = outdir.join("test.fastq");

    write_fastq(&input, &[
        ("read1", "ACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIII"),
        ("read2", "GCTAGCTAGCTAGCTAGCTA", "IIIIIIIIIIIIIIIIIIII"),
    ]);

    let output = Command::new(binary_path())
        .args(["--nozip", "--quiet", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run rastqc");

    assert!(output.status.success(), "Exit code was not 0: {:?}", output);

    let html = outdir.join("test_fastqc.html");
    assert!(html.exists(), "HTML report not created");

    let content = fs::read_to_string(&html).unwrap();
    assert!(content.contains("Basic Statistics"));
    assert!(content.contains("Per base sequence quality"));
    assert!(content.contains("RastQC"));

    cleanup(&outdir);
}

#[test]
fn test_single_file_zip() {
    let outdir = test_dir().join("single_zip");
    let _ = fs::create_dir_all(&outdir);
    let input = outdir.join("sample.fastq");

    write_fastq(&input, &[
        ("read1", "ACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIII"),
    ]);

    let output = Command::new(binary_path())
        .args(["--quiet", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run rastqc");

    assert!(output.status.success());

    let zip_path = outdir.join("sample_fastqc.zip");
    assert!(zip_path.exists(), "ZIP not created");
    assert!(fs::metadata(&zip_path).unwrap().len() > 100, "ZIP too small");

    // HTML should also be written alongside (for summary linking)
    let html_path = outdir.join("sample_fastqc.html");
    assert!(html_path.exists(), "Standalone HTML not created alongside ZIP");

    cleanup(&outdir);
}

#[test]
fn test_zip_with_extract() {
    let outdir = test_dir().join("extract");
    let _ = fs::create_dir_all(&outdir);
    let input = outdir.join("sample.fastq");

    write_fastq(&input, &[
        ("read1", "ACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIII"),
    ]);

    let output = Command::new(binary_path())
        .args(["--extract", "--quiet", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run");

    assert!(output.status.success());

    let extract_dir = outdir.join("sample_fastqc");
    assert!(extract_dir.join("fastqc_report.html").exists());
    assert!(extract_dir.join("fastqc_data.txt").exists());

    cleanup(&outdir);
}

#[test]
fn test_multiple_files_summary() {
    let outdir = test_dir().join("multi");
    let _ = fs::create_dir_all(&outdir);

    let input1 = outdir.join("a.fastq");
    let input2 = outdir.join("b.fastq");

    write_fastq(&input1, &[
        ("r1", "ACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIII"),
        ("r2", "GCTAGCTAGCTAGCTAGCTA", "IIIIIIIIIIIIIIIIIIII"),
    ]);
    write_fastq(&input2, &[
        ("r1", "AAAAAAAAAAAAAAAAAAAA", "IIIIIIIIIIIIIIIIIIII"),
        ("r2", "TTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIII"),
    ]);

    let output = Command::new(binary_path())
        .args(["--nozip", "--summary", "--quiet", "-o"])
        .arg(&outdir)
        .arg(&input1)
        .arg(&input2)
        .output()
        .expect("Failed to run");

    assert!(output.status.success());

    // Check individual reports
    assert!(outdir.join("a_fastqc.html").exists());
    assert!(outdir.join("b_fastqc.html").exists());

    // Check summary files
    let tsv_path = outdir.join("summary.tsv");
    assert!(tsv_path.exists(), "summary.tsv not created");

    let tsv = fs::read_to_string(&tsv_path).unwrap();
    assert!(tsv.contains("a.fastq"));
    assert!(tsv.contains("b.fastq"));
    assert!(tsv.contains("PASS") || tsv.contains("WARN") || tsv.contains("FAIL"));

    let summary_html = outdir.join("summary.html");
    assert!(summary_html.exists(), "summary.html not created");

    let html = fs::read_to_string(&summary_html).unwrap();
    assert!(html.contains("RastQC Summary"));
    assert!(html.contains("a.fastq"));
    assert!(html.contains("b.fastq"));

    cleanup(&outdir);
}

#[test]
fn test_text_data_modules_present() {
    let outdir = test_dir().join("text_data");
    let _ = fs::create_dir_all(&outdir);
    let input = outdir.join("test.fastq");

    write_fastq(&input, &[
        ("r1", "ACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIII"),
        ("r2", "GCTAGCTAGCTAGCTAGCTA", "BBBBBBBBBBBBBBBBBBBB"),
        ("r3", "NNNNNACTGACTGACTGACT", "IIIIIIIIIIIIIIIIIIII"),
    ]);

    let output = Command::new(binary_path())
        .args(["--extract", "--quiet", "--long-read", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run");

    assert!(output.status.success());

    let data_path = outdir.join("test_fastqc/fastqc_data.txt");
    let data = fs::read_to_string(&data_path).unwrap();

    // Verify all expected modules are present (12 core + 3 long-read with --long-read)
    assert!(data.contains(">>Basic Statistics"));
    assert!(data.contains(">>Per base sequence quality"));
    assert!(data.contains(">>Per tile sequence quality"));
    assert!(data.contains(">>Per sequence quality scores"));
    assert!(data.contains(">>Per base sequence content"));
    assert!(data.contains(">>Per sequence GC content"));
    assert!(data.contains(">>Per base N content"));
    assert!(data.contains(">>Sequence Length Distribution"));
    assert!(data.contains(">>Sequence Duplication Levels"));
    assert!(data.contains(">>Overrepresented sequences"));
    assert!(data.contains(">>Adapter Content"));
    assert!(data.contains(">>Kmer Content"));
    assert!(data.contains(">>Read Length N50"));
    assert!(data.contains(">>Quality Stratified Length"));
    assert!(data.contains(">>Homopolymer Content"));

    // Verify each module ends properly (12 original + 3 long-read = 15)
    let module_count = data.matches(">>END_MODULE").count();
    assert!(module_count >= 15, "Expected >=15 END_MODULE markers, got {}", module_count);

    cleanup(&outdir);
}

#[test]
fn test_gzip_input() {
    use std::io::Write;

    let outdir = test_dir().join("gzip");
    let _ = fs::create_dir_all(&outdir);
    let input = outdir.join("test.fastq.gz");

    // Write gzipped FASTQ
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(b"@read1\nACTGACTGACTG\n+\nIIIIIIIIIIII\n").unwrap();
    encoder.write_all(b"@read2\nGCTAGCTAGCTA\n+\nIIIIIIIIIIII\n").unwrap();
    let compressed = encoder.finish().unwrap();
    fs::write(&input, compressed).unwrap();

    let output = Command::new(binary_path())
        .args(["--nozip", "--quiet", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run");

    assert!(output.status.success());
    assert!(outdir.join("test_fastqc.html").exists());

    cleanup(&outdir);
}

#[test]
fn test_empty_fastq_error() {
    let outdir = test_dir().join("empty");
    let _ = fs::create_dir_all(&outdir);
    let input = outdir.join("empty.fastq");
    fs::write(&input, "").unwrap();

    let output = Command::new(binary_path())
        .args(["--nozip", "--quiet", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run");

    // Should still succeed (0 sequences), not crash
    assert!(output.status.success());

    cleanup(&outdir);
}

#[test]
fn test_nonexistent_file_error() {
    let output = Command::new(binary_path())
        .args(["--quiet", "/nonexistent/file.fastq"])
        .output()
        .expect("Failed to run");

    // Should print error but not crash
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Error") || stderr.contains("error") || stderr.contains("Cannot"));
}

#[test]
fn test_help_flag() {
    let output = Command::new(binary_path())
        .arg("--help")
        .output()
        .expect("Failed to run");

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("quality control tool"));
    assert!(stdout.contains("--outdir"));
    assert!(stdout.contains("--threads"));
    assert!(stdout.contains("--adapters"));
    assert!(stdout.contains("--contaminants"));
    assert!(stdout.contains("--summary"));
}

#[test]
fn test_version_flag() {
    let output = Command::new(binary_path())
        .arg("--version")
        .output()
        .expect("Failed to run");

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("0.1.0"));
}

#[test]
fn test_summary_tsv_format() {
    let outdir = test_dir().join("tsv_format");
    let _ = fs::create_dir_all(&outdir);

    let input = outdir.join("s1.fastq");
    write_fastq(&input, &[
        ("r1", "ACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIII"),
    ]);

    Command::new(binary_path())
        .args(["--nozip", "--summary", "--quiet", "-o"])
        .arg(&outdir)
        .arg(&input)
        .output()
        .expect("Failed to run");

    let tsv = fs::read_to_string(outdir.join("summary.tsv")).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();

    assert!(lines.len() >= 2, "Need header + data row");

    // Header should have Sample + modules + Total Sequences
    let header_cols: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(header_cols[0], "Sample");
    assert!(header_cols.last().unwrap().contains("Total Sequences"));

    // Data row
    let data_cols: Vec<&str> = lines[1].split('\t').collect();
    assert!(data_cols[0].contains("s1.fastq"));
    // Last col should be a number
    let seq_count: u64 = data_cols.last().unwrap().parse().unwrap();
    assert!(seq_count > 0);

    cleanup(&outdir);
}
