use anyhow::{bail, Context, Result};
use bzip2::read::BzDecoder;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use super::colorspace;
use super::Sequence;

/// Compression of a stream whose name is unknown (e.g. stdin), inferred from
/// leading magic bytes.
#[derive(Clone, Copy)]
enum StreamKind {
    Gzip,
    Bzip2,
    Plain,
}

pub struct FastqReader {
    reader: Box<dyn BufRead>,
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
            // MultiGzDecoder (not GzDecoder) is required so that multi-member
            // gzip streams are fully decoded. Large FASTQs compressed with
            // pigz/bgzip (or any chunked gzip) are concatenations of many gzip
            // members; GzDecoder stops silently after the first member, which
            // truncated large inputs to ~20k reads. See issue #3.
            Box::new(MultiGzDecoder::new(file))
        } else if name.ends_with(".bz2") {
            Box::new(BzDecoder::new(file))
        } else {
            Box::new(file)
        };

        Ok(FastqReader {
            reader: Box::new(BufReader::with_capacity(1024 * 1024, reader)),
            colorspace_detected: None,
        })
    }

    /// Create a reader from stdin.
    ///
    /// stdin carries no filename to infer compression from, so we sniff the
    /// leading magic bytes and transparently decompress gzip (including
    /// multi-member pigz/bgzip streams) and bzip2 input — making
    /// `zcat foo.fastq.gz | rastqc -` and `rastqc - < foo.fastq.gz` behave like
    /// passing the file directly. Falls back to raw bytes otherwise.
    pub fn from_stdin() -> Self {
        let mut buf = BufReader::with_capacity(1024 * 1024, io::stdin());

        // Peek (don't consume): BufReader::fill_buf exposes buffered bytes that
        // the subsequent reads — including the decoder's — still see.
        let kind = match buf.fill_buf() {
            Ok([0x1f, 0x8b, ..]) => StreamKind::Gzip,
            Ok([b'B', b'Z', b'h', ..]) => StreamKind::Bzip2,
            _ => StreamKind::Plain,
        };

        let reader: Box<dyn BufRead> = match kind {
            StreamKind::Gzip => Box::new(BufReader::with_capacity(
                1024 * 1024,
                MultiGzDecoder::new(buf),
            )),
            StreamKind::Bzip2 => {
                Box::new(BufReader::with_capacity(1024 * 1024, BzDecoder::new(buf)))
            }
            StreamKind::Plain => Box::new(buf),
        };

        FastqReader {
            reader,
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
                    if e.kind() == io::ErrorKind::UnexpectedEof
                        || e.to_string().contains("unexpected end of file")
                    {
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
            if e.kind() == io::ErrorKind::UnexpectedEof
                || e.to_string().contains("unexpected end of file")
            {
                eprintln!(
                    "Warning: File truncated during sequence read. Skipping last partial record."
                );
                return Ok(None);
            }
            return Err(e).context("Error reading sequence line");
        }
        let sequence = sequence_line.trim().as_bytes().to_vec();

        // Read separator line (+)
        let mut sep_line = String::with_capacity(16);
        let res = self.reader.read_line(&mut sep_line);
        if let Err(e) = res {
            if e.kind() == io::ErrorKind::UnexpectedEof
                || e.to_string().contains("unexpected end of file")
            {
                eprintln!(
                    "Warning: File truncated during separator read. Skipping last partial record."
                );
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
            if e.kind() == io::ErrorKind::UnexpectedEof
                || e.to_string().contains("unexpected end of file")
            {
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

        if sequence.len() != quality.len() {
            bail!(
                "Malformed FASTQ record '{}': sequence length ({}) does not match quality length ({})",
                header,
                sequence.len(),
                quality.len()
            );
        }

        Ok(Some(Sequence {
            header,
            sequence,
            quality,
            filtered,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn reader_for(data: &str) -> FastqReader {
        FastqReader {
            reader: Box::new(BufReader::new(Cursor::new(data.as_bytes().to_vec()))),
            colorspace_detected: None,
        }
    }

    #[test]
    fn well_formed_record_parses() {
        let mut r = reader_for("@read1\nACGT\n+\nIIII\n");
        let seq = r.next_sequence().unwrap().unwrap();
        assert_eq!(seq.sequence, b"ACGT");
        assert_eq!(seq.quality, b"IIII");
    }

    #[test]
    fn mismatched_sequence_and_quality_length_is_rejected() {
        let mut r = reader_for("@read1\nACGT\n+\nIII\n");
        let err = r.next_sequence().unwrap_err();
        assert!(err.to_string().contains("does not match"));
    }
}
