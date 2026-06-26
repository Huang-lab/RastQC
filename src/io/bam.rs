use anyhow::{Context, Result};
use noodles::sam::alignment::record::QualityScores as QualityScoresTrait;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::bam;
use noodles::sam;

use super::Sequence;

pub enum BamReader {
    Bam {
        reader: bam::io::Reader<noodles::bgzf::Reader<BufReader<File>>>,
    },
    Sam {
        reader: sam::io::Reader<BufReader<File>>,
    },
}

impl BamReader {
    pub fn open(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("Cannot open BAM: {}", path.display()))?;
        let mut reader = bam::io::reader::Builder.build_from_reader(BufReader::new(file));
        let _header = reader.read_header()?;

        Ok(BamReader::Bam { reader })
    }

    pub fn open_sam(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("Cannot open SAM: {}", path.display()))?;
        let mut reader = sam::io::Reader::new(BufReader::new(file));
        let _header = reader.read_header()?;

        Ok(BamReader::Sam { reader })
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        match self {
            BamReader::Bam { reader } => {
                let mut record = bam::Record::default();
                match reader.read_record(&mut record) {
                    Ok(0) => Ok(None),
                    Ok(_) => bam_record_to_sequence(&record),
                    Err(e) => Err(e.into()),
                }
            }
            BamReader::Sam { reader } => {
                let mut record = sam::Record::default();
                match reader.read_record(&mut record) {
                    Ok(0) => Ok(None),
                    Ok(_) => sam_record_to_sequence(&record),
                    Err(e) => Err(e.into()),
                }
            }
        }
    }
}

fn bam_record_to_sequence(record: &bam::Record) -> Result<Option<Sequence>> {
    let flags = record.flags();

    if flags.is_secondary() || flags.is_supplementary() {
        return Ok(None);
    }

    let name = record.name().map(|n| n.to_string()).unwrap_or_default();

    // Extract sequence bases
    let seq_data = record.sequence();
    let seq_len = seq_data.len();
    let mut sequence = Vec::with_capacity(seq_len);
    for i in 0..seq_len {
        let base = seq_data.get(i).unwrap_or(b'N');
        sequence.push(match base.to_ascii_uppercase() {
            b'A' => b'A',
            b'C' => b'C',
            b'G' => b'G',
            b'T' => b'T',
            _ => b'N',
        });
    }

    // Extract quality scores - stored as raw phred values in BAM
    let qual_data = record.quality_scores();
    let qual_bytes: &[u8] = qual_data.as_ref();
    let quality: Vec<u8> = qual_bytes.iter().map(|&q| q + 33).collect();

    if sequence.is_empty() {
        return Ok(None);
    }

    Ok(Some(Sequence {
        header: format!("@{}", name),
        sequence,
        quality,
        filtered: flags.is_qc_fail(),
    }))
}

fn sam_record_to_sequence(record: &sam::Record) -> Result<Option<Sequence>> {
    let flags = record.flags()?;

    if flags.is_secondary() || flags.is_supplementary() {
        return Ok(None);
    }

    let name = record.name().map(|n| n.to_string()).unwrap_or_default();

    let seq_data = record.sequence();
    let seq_len = seq_data.len();
    let mut sequence = Vec::with_capacity(seq_len);
    for i in 0..seq_len {
        let base = match seq_data.get(i) {
            Some(b) => match b.to_ascii_uppercase() {
                b'A' => b'A',
                b'C' => b'C',
                b'G' => b'G',
                b'T' => b'T',
                _ => b'N',
            },
            None => b'N',
        };
        sequence.push(base);
    }

    let qual_scores = record.quality_scores();
    let mut quality = Vec::with_capacity(seq_len);
    for score in qual_scores.iter() {
        match score {
            Ok(q) => quality.push(q + 33),
            Err(_) => quality.push(33),
        }
    }

    if sequence.is_empty() {
        return Ok(None);
    }

    Ok(Some(Sequence {
        header: format!("@{}", name),
        sequence,
        quality,
        filtered: flags.is_qc_fail(),
    }))
}
