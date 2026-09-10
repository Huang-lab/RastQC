//! Record-aligned block reading and borrowed parsing for FASTQ.
//!
//! The per-record reader ([`super::fastq::FastqReader`]) materializes every
//! line into owned `String`/`Vec<u8>` buffers. That costs roughly seven heap
//! allocations, several memcpys and a full UTF-8 validation per read — and on
//! the parallel pipeline all of it lands on the *single* reader thread, which
//! then caps total throughput no matter how many workers are running.
//!
//! This module splits that work: the reader thread only decompresses and cuts
//! the stream at record boundaries, handing whole blocks of bytes to the
//! workers, which parse records as slices borrowed from the block. Parsing
//! then scales with the worker count and the hot loop allocates nothing.

use anyhow::{anyhow, bail, Result};
use memchr::{memchr, memchr_iter};
use std::io::{self, Read};

/// Target decompressed bytes per block handed to a worker.
///
/// Large enough that the per-block channel hand-off and boundary scan are
/// negligible against the parsing work inside it, small enough that a pool of
/// one block per worker stays a rounding error against the module state each
/// worker already carries.
pub const BLOCK_SIZE: usize = 1024 * 1024;

/// Bytes requested per `read` call while filling a block.
const READ_CHUNK: usize = 128 * 1024;

/// Hard cap on a single block. A block must contain at least one whole record,
/// so this doubles as the cap on the largest single record — set far above any
/// real read (including ONT ultra-long reads) while still bounding worst-case
/// memory on corrupt or misdetected-binary input, the same role
/// [`super::MAX_LINE_LEN`] plays for the per-record reader.
const MAX_BLOCK_LEN: usize = 256 * 1024 * 1024;

/// A structural parse failure on the block fast path.
///
/// Kept as its own type so the caller can distinguish "this input doesn't fit
/// the strict four-lines-per-record fast path" from a genuine I/O failure and
/// retry the file on the more tolerant per-record reader, rather than failing
/// the whole run.
#[derive(Debug)]
pub struct BlockParseError(pub String);

impl std::fmt::Display for BlockParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for BlockParseError {}

/// Reads a FASTQ stream in blocks that each hold a whole number of records.
pub struct FastqBlockReader {
    inner: Box<dyn Read + Send>,
    /// Trailing bytes of an incomplete record, carried into the next block.
    carry: Vec<u8>,
    eof: bool,
}

impl FastqBlockReader {
    pub fn new(inner: Box<dyn Read + Send>) -> Self {
        FastqBlockReader {
            inner,
            carry: Vec::new(),
            eof: false,
        }
    }

    /// Fill `buf` with a whole number of four-line FASTQ records.
    ///
    /// Returns `Ok(false)` once the stream is exhausted. A trailing partial
    /// record is carried into the next call; one left at EOF is dropped with a
    /// warning, matching the per-record reader's behaviour on truncated input.
    pub fn next_block(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        buf.clear();
        buf.append(&mut self.carry);

        let mut target = BLOCK_SIZE;
        loop {
            while buf.len() < target && !self.eof {
                let start = buf.len();
                buf.resize(start + READ_CHUNK, 0);
                match self.inner.read(&mut buf[start..]) {
                    Ok(0) => {
                        buf.truncate(start);
                        self.eof = true;
                    }
                    Ok(n) => buf.truncate(start + n),
                    Err(ref e) if e.kind() == io::ErrorKind::Interrupted => buf.truncate(start),
                    Err(e) => {
                        buf.truncate(start);
                        // A truncated compressed stream still has usable
                        // records ahead of the cut; warn and analyze them
                        // rather than discarding the whole file.
                        if e.kind() == io::ErrorKind::UnexpectedEof
                            || e.to_string().contains("unexpected end of file")
                        {
                            eprintln!("Warning: Unexpected end of compressed file. Processing remaining data.");
                            self.eof = true;
                        } else {
                            return Err(e.into());
                        }
                    }
                }
            }

            let cut = last_record_boundary(buf);
            if cut > 0 {
                self.carry.extend_from_slice(&buf[cut..]);
                buf.truncate(cut);
                return Ok(true);
            }

            // Nothing complete in hand yet.
            if self.eof {
                if !trim_ascii(buf).is_empty() {
                    eprintln!("Warning: File truncated mid-record. Skipping last partial record.");
                }
                buf.clear();
                return Ok(false);
            }
            if target >= MAX_BLOCK_LEN {
                bail!(
                    "no complete FASTQ record found within {} bytes (malformed or corrupted input?)",
                    MAX_BLOCK_LEN
                );
            }
            target = target.saturating_mul(2).min(MAX_BLOCK_LEN);
        }
    }
}

/// Byte offset just past the last complete four-line record in `buf`.
///
/// Relies on `buf` starting on a record boundary, which `next_block`
/// maintains by carrying any partial tail into the following block. Returns 0
/// when no complete record is present.
fn last_record_boundary(buf: &[u8]) -> usize {
    let mut lines = 0usize;
    let mut cut = 0usize;
    for pos in memchr_iter(b'\n', buf) {
        lines += 1;
        if lines % 4 == 0 {
            cut = pos + 1;
        }
    }
    cut
}

/// One FASTQ record, borrowed from the block it was parsed out of.
pub struct RawRecord<'a> {
    pub header: &'a [u8],
    pub sequence: &'a [u8],
    pub quality: &'a [u8],
}

/// Iterator over the records in one block.
pub struct RecordIter<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> RecordIter<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        RecordIter { buf, pos: 0 }
    }

    fn next_line(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.buf.len() {
            return None;
        }
        let rest = &self.buf[self.pos..];
        let (line, advance) = match memchr(b'\n', rest) {
            Some(i) => (&rest[..i], i + 1),
            None => (rest, rest.len()),
        };
        self.pos += advance;
        Some(line)
    }
}

impl<'a> Iterator for RecordIter<'a> {
    type Item = Result<RawRecord<'a>, BlockParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Blank lines between records are tolerated ahead of a header, as the
        // per-record reader does.
        let header = loop {
            let line = trim_ascii(self.next_line()?);
            if !line.is_empty() {
                break line;
            }
        };
        if header[0] != b'@' {
            return Some(Err(BlockParseError(format!(
                "Expected FASTQ header starting with '@', got: {}",
                String::from_utf8_lossy(header)
            ))));
        }

        let Some(sequence) = self.next_line().map(trim_ascii) else {
            return Some(Err(BlockParseError(
                "File truncated during sequence read".into(),
            )));
        };
        let Some(separator) = self.next_line().map(trim_ascii) else {
            return Some(Err(BlockParseError(
                "File truncated during separator read".into(),
            )));
        };
        if !separator.starts_with(b"+") {
            return Some(Err(BlockParseError(format!(
                "Expected '+' separator, got: {}",
                String::from_utf8_lossy(separator)
            ))));
        }
        let Some(quality) = self.next_line().map(trim_ascii) else {
            return Some(Err(BlockParseError(
                "File truncated during quality score read".into(),
            )));
        };
        if sequence.len() != quality.len() {
            return Some(Err(BlockParseError(format!(
                "Malformed FASTQ record '{}': sequence length ({}) does not match quality length ({})",
                String::from_utf8_lossy(header),
                sequence.len(),
                quality.len()
            ))));
        }

        Some(Ok(RawRecord {
            header,
            sequence,
            quality,
        }))
    }
}

/// `[u8]::trim_ascii`, which is newer than this crate's MSRV (Rust 1.70).
///
/// Matches the `str::trim` the per-record reader applies to every line, so
/// `\r\n` line endings and stray padding are handled identically here.
pub(crate) fn trim_ascii(mut s: &[u8]) -> &[u8] {
    while let [first, rest @ ..] = s {
        if first.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    while let [rest @ .., last] = s {
        if last.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    s
}

/// Non-UTF-8 bytes in a header are a hard error, as on the per-record reader.
pub(crate) fn header_str(header: &[u8]) -> Result<&str, BlockParseError> {
    std::str::from_utf8(header)
        .map_err(|e| BlockParseError(format!("FASTQ header is not valid UTF-8: {e}")))
}

/// Whether the first record of `path` looks like SOLiD colorspace.
///
/// Colorspace decoding changes a record's length relationship with its quality
/// line, so those (legacy) files stay on the per-record reader instead of
/// growing a special case in the hot loop.
pub fn first_record_is_colorspace(path: &std::path::Path) -> Result<bool> {
    let mut reader = FastqBlockReader::new(super::fastq::open_decompressed(path)?);
    let mut buf = Vec::new();
    if !reader.next_block(&mut buf)? {
        return Ok(false);
    }
    match RecordIter::new(&buf).next() {
        Some(Ok(record)) => Ok(super::colorspace::is_colorspace(record.sequence)),
        Some(Err(e)) => Err(anyhow!(e)),
        None => Ok(false),
    }
}
