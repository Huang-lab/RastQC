use anyhow::{bail, Context, Result};
use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use super::colorspace;
use super::Sequence;

pub struct FastqReader {
    reader: Box<dyn BufRead>,
    line_buf: String,
    colorspace_detected: Option<bool>,
}

impl FastqReader {
    /// Open a FASTQ file (optionally compressed with gzip or bzip2).
    pub fn open(path: &Path) -> Result<Self> {
        let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_lowercase();

        let file =
            File::open(path).with_context(|| format!("Cannot open file: {}", path.display()))?;

        let reader: Box<dyn Read> = if name.ends_with(".gz") {
            Box::new(GzDecoder::new(file))
        } else if name.ends_with(".bz2") {
            Box::new(BzDecoder::new(file))
        } else {
            Box::new(file)
        };

        Ok(FastqReader {
            reader: Box::new(BufReader::with_capacity(1024 * 1024, reader)),
            line_buf: String::with_capacity(512),
            colorspace_detected: None,
        })
    }

    /// Create a reader from stdin.
    pub fn from_stdin() -> Self {
        FastqReader {
            reader: Box::new(BufReader::with_capacity(1024 * 1024, io::stdin())),
            line_buf: String::with_capacity(512),
            colorspace_detected: None,
        }
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        let mut header_line = String::with_capacity(512);
        
        // Read header line (starts with @)
        loop {
            header_line.clear();
            match self.reader.read_line(&mut header_line) {
                Ok(0) => return Ok(None),
                Ok(_) => {
                    let trimmed = header_line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    if !trimmed.starts_with('@') {
                        bail!("Expected FASTQ header starting with '@', got: {}", trimmed);
                    }
                    break;
                }
                Err(e) => {
                    if e.kind() == io::ErrorKind::UnexpectedEof || e.to_string().contains("unexpected end of file") {
                        eprintln!("Warning: Unexpected end of compressed file. Processing remaining data.");
                        return Ok(None);
                    }
                    return Err(e).context("Error reading FASTQ header");
                }
            }
        }

        let header = header_line.trim().to_string();
        let filtered = header.contains(":Y:");

        // Read sequence line
        let mut sequence_line = String::with_capacity(512);
        let res = self.reader.read_line(&mut sequence_line);
        if let Err(e) = res {
            if e.kind() == io::ErrorKind::UnexpectedEof || e.to_string().contains("unexpected end of file") {
                eprintln!("Warning: File truncated during sequence read. Skipping last partial record.");
                return Ok(None);
            }
            return Err(e).context("Error reading sequence line");
        }
        let sequence = sequence_line.trim().as_bytes().to_vec();

        // Read separator line (+)
        let mut sep_line = String::with_capacity(16);
        let res = self.reader.read_line(&mut sep_line);
        if let Err(e) = res {
            if e.kind() == io::ErrorKind::UnexpectedEof || e.to_string().contains("unexpected end of file") {
                eprintln!("Warning: File truncated during separator read. Skipping last partial record.");
                return Ok(None);
            }
            return Err(e).context("Error reading separator line");
        }
        if !sep_line.trim().starts_with('+') {
            bail!("Expected '+' separator, got: {}", sep_line.trim());
        }

        // Read quality line
        let mut quality_line = String::with_capacity(512);
        let res = self.reader.read_line(&mut quality_line);
        if let Err(e) = res {
            if e.kind() == io::ErrorKind::UnexpectedEof || e.to_string().contains("unexpected end of file") {
                eprintln!("Warning: File truncated during quality score read. Skipping last partial record.");
                return Ok(None);
            }
            return Err(e).context("Error reading quality line");
        }
        let quality = quality_line.trim().as_bytes().to_vec();

        // Auto-detect colorspace on first sequence
        if self.colorspace_detected.is_none() {
            self.colorspace_detected = Some(colorspace::is_colorspace(&sequence));
        }

        // Decode colorspace if detected
        let sequence = if self.colorspace_detected == Some(true) {
            colorspace::decode_colorspace(&sequence).unwrap_or(sequence)
        } else {
            sequence
        };

        Ok(Some(Sequence {
            header,
            sequence,
            quality,
            filtered,
        }))
    }
}
