use crate::io::Sequence;
use anyhow::Result;
use bzip2::read::BzDecoder;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub struct FastaReader {
    reader: Box<dyn BufRead>,
    header: Option<String>,
    /// Reusable scratch buffer for `read_line_bounded`, shared across both
    /// the header-scan and sequence-accumulation loops.
    line_scratch: Vec<u8>,
}

impl FastaReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(&path)?;
        let path_ref = path.as_ref();
        let name = path_ref
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_lowercase();

        let reader: Box<dyn Read> = if name.ends_with(".gz") {
            // Use MultiGzDecoder to fully decode multi-member gzip streams
            // (pigz/bgzip output). See issue #3.
            Box::new(MultiGzDecoder::new(file))
        } else if name.ends_with(".bz2") {
            Box::new(BzDecoder::new(file))
        } else {
            Box::new(file)
        };

        Ok(FastaReader {
            reader: Box::new(BufReader::new(reader)),
            header: None,
            line_scratch: Vec::new(),
        })
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        let mut line = String::new();

        // If we don't have a header, find the first one
        if self.header.is_none() {
            loop {
                line.clear();
                if super::read_line_bounded(&mut *self.reader, &mut line, &mut self.line_scratch)?
                    == 0
                {
                    return Ok(None);
                }
                let trimmed = line.trim();
                if let Some(h) = trimmed.strip_prefix('>') {
                    self.header = Some(h.to_string());
                    break;
                }
            }
        }

        let header = self.header.take().unwrap();
        let mut sequence = Vec::new();

        loop {
            line.clear();
            if super::read_line_bounded(&mut *self.reader, &mut line, &mut self.line_scratch)? == 0
            {
                // End of file
                break;
            }
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(h) = trimmed.strip_prefix('>') {
                // Next header found, save it and return current sequence
                self.header = Some(h.to_string());
                break;
            }
            sequence.extend_from_slice(trimmed.as_bytes());
        }

        // `header` was already taken from `self.header` above, so by this
        // point we always have a real record to return — even a trailing
        // `>lastread` with no sequence lines before EOF is a valid
        // (zero-length) record and must not be silently dropped.
        let len = sequence.len();
        Ok(Some(Sequence {
            header,
            sequence,
            quality: vec![b'I'; len], // Default high quality 'I' (40) for FASTA
            filtered: false,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reader_for(contents: &str) -> FastaReader {
        FastaReader {
            reader: Box::new(BufReader::new(std::io::Cursor::new(contents.as_bytes().to_vec()))),
            header: None,
            line_scratch: Vec::new(),
        }
    }

    #[test]
    fn reads_a_normal_multi_record_fasta() {
        let mut r = reader_for(">seq1\nACGT\n>seq2\nTTTT\n");
        let s1 = r.next_sequence().unwrap().unwrap();
        assert_eq!(s1.header, "seq1");
        assert_eq!(s1.sequence, b"ACGT");
        let s2 = r.next_sequence().unwrap().unwrap();
        assert_eq!(s2.header, "seq2");
        assert_eq!(s2.sequence, b"TTTT");
        assert!(r.next_sequence().unwrap().is_none());
    }

    #[test]
    fn trailing_header_with_no_sequence_before_eof_is_still_returned() {
        // Regression test: a file ending in a bare `>header\n` with no
        // sequence line (or a truncated file cut off right after the
        // header) used to be silently dropped instead of yielding a
        // zero-length record.
        let mut r = reader_for(">seq1\nACGT\n>lastread\n");
        let s1 = r.next_sequence().unwrap().unwrap();
        assert_eq!(s1.header, "seq1");

        let s2 = r.next_sequence().unwrap().unwrap();
        assert_eq!(s2.header, "lastread");
        assert!(s2.sequence.is_empty());

        assert!(r.next_sequence().unwrap().is_none());
    }

    #[test]
    fn empty_file_returns_none() {
        let mut r = reader_for("");
        assert!(r.next_sequence().unwrap().is_none());
    }

    #[test]
    fn multiline_sequence_is_concatenated() {
        let mut r = reader_for(">seq1\nACGT\nTTTT\n\nGGGG\n");
        let s1 = r.next_sequence().unwrap().unwrap();
        assert_eq!(s1.sequence, b"ACGTTTTTGGGG");
    }
}
