use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use flate2::read::MultiGzDecoder;
use bzip2::read::BzDecoder;
use crate::io::Sequence;

pub struct FastaReader {
    reader: Box<dyn BufRead>,
    header: Option<String>,
    buffer: Vec<u8>,
}

impl FastaReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(&path)?;
        let path_ref = path.as_ref();
        let name = path_ref.file_name().unwrap_or_default().to_string_lossy().to_lowercase();
        
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
            buffer: Vec::new(),
        })
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        let mut line = String::new();
        
        // If we don't have a header, find the first one
        if self.header.is_none() {
            loop {
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    return Ok(None);
                }
                let trimmed = line.trim();
                if trimmed.starts_with('>') {
                    self.header = Some(trimmed[1..].to_string());
                    break;
                }
            }
        }

        let header = self.header.take().unwrap();
        let mut sequence = Vec::new();
        
        loop {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                // End of file
                break;
            }
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('>') {
                // Next header found, save it and return current sequence
                self.header = Some(trimmed[1..].to_string());
                break;
            }
            sequence.extend_from_slice(trimmed.as_bytes());
        }

        if sequence.is_empty() && self.header.is_none() {
             return Ok(None);
        }

        let len = sequence.len();
        Ok(Some(Sequence {
            header,
            sequence,
            quality: vec![b'I'; len], // Default high quality 'I' (40) for FASTA
            filtered: false,
        }))
    }
}
