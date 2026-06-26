mod bam;
pub mod colorspace;
mod fast5;
mod fasta;
mod fastq;
mod pod5;

use anyhow::{bail, Result};
use std::path::Path;

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
