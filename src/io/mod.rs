mod bam;
pub mod colorspace;
mod fast5;
mod fasta;
mod fastq;
mod pod5;

use anyhow::{bail, Result};
use std::io::{self, BufRead};
use std::path::Path;

/// Generous but bounded per-line cap for text-format readers (FASTQ/FASTA).
///
/// `BufRead::read_line` grows its output buffer without limit until it finds
/// a newline (or hits EOF), so a corrupted/truncated file, or a binary file
/// misidentified as text, with a very long or missing newline can allocate
/// unbounded memory before any of the readers' own record-level validation
/// gets a chance to run. Real long-read sequencing (ONT ultra-long reads) can
/// exceed a few megabases, so this is set far above any realistic read length
/// while still bounding worst-case memory use.
const MAX_LINE_LEN: usize = 256 * 1024 * 1024;

/// Like `BufRead::read_line`, but errors instead of growing `buf` without
/// bound if no newline is found within `MAX_LINE_LEN` bytes.
///
/// `scratch` is a caller-owned byte buffer reused across calls (cleared each
/// call): without it, every call would allocate a fresh `Vec` from scratch,
/// throwing away the capacity-reuse the callers already rely on for their
/// per-record `String` buffers (one call per line, so this matters — 4x per
/// FASTQ record, 1x per FASTA line).
pub(crate) fn read_line_bounded<R: BufRead + ?Sized>(
    reader: &mut R,
    buf: &mut String,
    scratch: &mut Vec<u8>,
) -> io::Result<usize> {
    read_line_bounded_with_limit(reader, buf, scratch, MAX_LINE_LEN)
}

fn read_line_bounded_with_limit<R: BufRead + ?Sized>(
    reader: &mut R,
    buf: &mut String,
    scratch: &mut Vec<u8>,
    max_len: usize,
) -> io::Result<usize> {
    scratch.clear();
    loop {
        let (found_newline, consumed) = {
            let available = reader.fill_buf()?;
            if available.is_empty() {
                (true, 0)
            } else if let Some(pos) = available.iter().position(|&b| b == b'\n') {
                scratch.extend_from_slice(&available[..=pos]);
                (true, pos + 1)
            } else {
                scratch.extend_from_slice(available);
                (false, available.len())
            }
        };
        reader.consume(consumed);
        if scratch.len() > max_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "line exceeds maximum length of {} bytes (malformed or corrupted input?)",
                    max_len
                ),
            ));
        }
        if found_newline {
            break;
        }
    }
    let s =
        std::str::from_utf8(scratch).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let n = s.len();
    buf.push_str(s);
    Ok(n)
}

#[cfg(test)]
mod read_line_bounded_tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn reads_normal_line_including_newline() {
        let mut reader = Cursor::new(b"hello\nworld\n".to_vec());
        let mut buf = String::new();
        let mut scratch = Vec::new();
        let n = read_line_bounded_with_limit(&mut reader, &mut buf, &mut scratch, 1024).unwrap();
        assert_eq!(buf, "hello\n");
        assert_eq!(n, 6);
    }

    #[test]
    fn reads_final_line_without_trailing_newline() {
        let mut reader = Cursor::new(b"hello".to_vec());
        let mut buf = String::new();
        let mut scratch = Vec::new();
        let n = read_line_bounded_with_limit(&mut reader, &mut buf, &mut scratch, 1024).unwrap();
        assert_eq!(buf, "hello");
        assert_eq!(n, 5);
    }

    #[test]
    fn returns_zero_at_eof() {
        let mut reader = Cursor::new(Vec::new());
        let mut buf = String::new();
        let mut scratch = Vec::new();
        let n = read_line_bounded_with_limit(&mut reader, &mut buf, &mut scratch, 1024).unwrap();
        assert_eq!(n, 0);
        assert_eq!(buf, "");
    }

    #[test]
    fn errors_instead_of_growing_unbounded_when_no_newline_within_limit() {
        // A line with no newline that exceeds the cap must error, not
        // silently keep allocating (the scenario this guards against: a
        // corrupted/binary file misidentified as FASTQ/FASTA with a very
        // long or missing newline).
        let data = vec![b'A'; 10_000];
        let mut reader = Cursor::new(data);
        let mut buf = String::new();
        let mut scratch = Vec::new();
        let err =
            read_line_bounded_with_limit(&mut reader, &mut buf, &mut scratch, 100).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("exceeds maximum length"));
    }

    #[test]
    fn a_line_within_the_limit_across_multiple_fill_buf_chunks_still_succeeds() {
        // Cursor's fill_buf returns everything at once, so wrap in a reader
        // that yields one byte per call to exercise the multi-iteration
        // accumulation path (real BufReaders hand back whatever happens to
        // be in their internal buffer, not necessarily the whole line).
        struct OneByteAtATime(Cursor<Vec<u8>>);
        impl BufRead for OneByteAtATime {
            fn fill_buf(&mut self) -> io::Result<&[u8]> {
                let buf = self.0.fill_buf()?;
                let n = buf.len().min(1);
                Ok(&buf[..n])
            }
            fn consume(&mut self, amt: usize) {
                self.0.consume(amt)
            }
        }
        impl std::io::Read for OneByteAtATime {
            fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
                self.0.read(out)
            }
        }

        let mut reader = OneByteAtATime(Cursor::new(b"hello world\n".to_vec()));
        let mut buf = String::new();
        let mut scratch = Vec::new();
        let n = read_line_bounded_with_limit(&mut reader, &mut buf, &mut scratch, 1024).unwrap();
        assert_eq!(buf, "hello world\n");
        assert_eq!(n, 12);
    }
}

/// A single sequence record
#[derive(Debug, Clone)]
pub struct Sequence {
    pub header: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub filtered: bool,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[allow(dead_code)]
    pub fn gc_percent(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc = self
            .sequence
            .iter()
            .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count();
        let valid = self
            .sequence
            .iter()
            .filter(|&&b| b != b'N' && b != b'n')
            .count();
        if valid == 0 {
            return 0.0;
        }
        (gc as f64 / valid as f64) * 100.0
    }
}

/// Unified reader for FASTQ, BAM, SAM, Fast5, POD5, FASTA files
pub enum SequenceReader {
    Fastq(fastq::FastqReader),
    Bam(bam::BamReader),
    Fast5(fast5::Fast5Reader),
    Pod5(pod5::Pod5Reader),
    Fasta(fasta::FastaReader),
}

impl SequenceReader {
    pub fn open(path: &Path) -> Result<Self> {
        let name = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_lowercase();

        if name.ends_with(".fast5") {
            Ok(SequenceReader::Fast5(fast5::Fast5Reader::open(path)?))
        } else if name.ends_with(".pod5") {
            Ok(SequenceReader::Pod5(pod5::Pod5Reader::open(path)?))
        } else if name.ends_with(".bam") {
            Ok(SequenceReader::Bam(bam::BamReader::open(path)?))
        } else if name.ends_with(".sam") {
            Ok(SequenceReader::Bam(bam::BamReader::open_sam(path)?))
        } else if name.ends_with(".fastq")
            || name.ends_with(".fq")
            || name.ends_with(".fastq.gz")
            || name.ends_with(".fq.gz")
            || name.ends_with(".fastq.bz2")
            || name.ends_with(".fq.bz2")
        {
            Ok(SequenceReader::Fastq(fastq::FastqReader::open(path)?))
        } else if name.ends_with(".fasta")
            || name.ends_with(".fa")
            || name.ends_with(".fasta.gz")
            || name.ends_with(".fa.gz")
            || name.ends_with(".fasta.bz2")
            || name.ends_with(".fa.bz2")
        {
            Ok(SequenceReader::Fasta(fasta::FastaReader::open(path)?))
        } else {
            // Try FASTQ as default
            match fastq::FastqReader::open(path) {
                Ok(r) => Ok(SequenceReader::Fastq(r)),
                Err(_) => match fasta::FastaReader::open(path) {
                    Ok(r) => Ok(SequenceReader::Fasta(r)),
                    Err(_) => bail!(
                        "Unrecognized file format: {}. Supported: .fastq, .fq, .fasta, .fa, .bam, .sam, .fast5, .pod5 (with optional .gz/.bz2 compression)",
                        path.display()
                    ),
                },
            }
        }
    }

    /// Create a reader that reads FASTQ from stdin.
    pub fn from_stdin() -> Self {
        SequenceReader::Fastq(fastq::FastqReader::from_stdin())
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        match self {
            SequenceReader::Fastq(r) => r.next_sequence(),
            SequenceReader::Bam(r) => r.next_sequence(),
            SequenceReader::Fast5(r) => r.next_sequence(),
            SequenceReader::Pod5(r) => r.next_sequence(),
            SequenceReader::Fasta(r) => r.next_sequence(),
        }
    }
}
