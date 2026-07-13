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

    /// Reads the next primary alignment record with a non-empty sequence,
    /// transparently skipping secondary/supplementary/empty-SEQ records.
    ///
    /// `bam_record_to_sequence`/`sam_record_to_sequence` return `Ok(None)` for
    /// records that should be skipped (not just at true end-of-stream), so a
    /// single `read_record` call cannot distinguish "skip this record" from
    /// "no more records" — this loops until a real record or EOF is found.
    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        match self {
            BamReader::Bam { reader } => {
                let mut record = bam::Record::default();
                loop {
                    match reader.read_record(&mut record) {
                        Ok(0) => return Ok(None),
                        Ok(_) => {
                            if let Some(seq) = bam_record_to_sequence(&record)? {
                                return Ok(Some(seq));
                            }
                        }
                        Err(e) => return Err(e.into()),
                    }
                }
            }
            BamReader::Sam { reader } => {
                let mut record = sam::Record::default();
                loop {
                    match reader.read_record(&mut record) {
                        Ok(0) => return Ok(None),
                        Ok(_) => {
                            if let Some(seq) = sam_record_to_sequence(&record)? {
                                return Ok(Some(seq));
                            }
                        }
                        Err(e) => return Err(e.into()),
                    }
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

    // Extract quality scores - stored as raw phred values in BAM.
    // noodles doesn't validate that QUAL bytes stay within the legal
    // 0-93 Phred range, so a malformed/corrupted BAM can carry a byte
    // >= 223 here; saturate instead of panicking (debug) or silently
    // wrapping to a bogus low value (release).
    let qual_data = record.quality_scores();
    let qual_bytes: &[u8] = qual_data.as_ref();
    let quality: Vec<u8> = qual_bytes.iter().map(|&q| q.saturating_add(33)).collect();

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
            Ok(q) => quality.push(q.saturating_add(33)),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_sam(name: &str, contents: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!("rastqc_test_{}_{}.sam", std::process::id(), name));
        let mut file = File::create(&path).unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        path
    }

    /// Regression test: a secondary/supplementary alignment in the middle of
    /// a SAM/BAM file must not be mistaken for end-of-stream. Before the fix,
    /// `bam_record_to_sequence`/`sam_record_to_sequence` returning `Ok(None)`
    /// to mean "skip this record" was indistinguishable from true EOF at the
    /// `next_sequence` call site, silently truncating the read after the
    /// first secondary alignment.
    #[test]
    fn skips_secondary_alignment_without_truncating_stream() {
        let sam = "@HD\tVN:1.6\tSO:unsorted\n\
                    @SQ\tSN:chr1\tLN:1000\n\
                    read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
                    read1\t256\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
                    read2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tGGGG\tIIII\n";
        let path = write_temp_sam("skip_secondary", sam);
        let mut reader = BamReader::open_sam(&path).unwrap();

        let mut names = Vec::new();
        while let Some(seq) = reader.next_sequence().unwrap() {
            names.push(seq.header);
        }
        let _ = std::fs::remove_file(&path);

        assert_eq!(names, vec!["@read1".to_string(), "@read2".to_string()]);
    }

    /// Same truncation hazard, but for a record with an empty `SEQ` (`*`)
    /// in the middle of the file rather than a secondary alignment flag.
    #[test]
    fn skips_empty_sequence_record_without_truncating_stream() {
        let sam = "@HD\tVN:1.6\tSO:unsorted\n\
                    @SQ\tSN:chr1\tLN:1000\n\
                    read1\t4\tchr1\t0\t0\t*\t*\t0\t0\t*\t*\n\
                    read2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tGGGG\tIIII\n";
        let path = write_temp_sam("skip_empty_seq", sam);
        let mut reader = BamReader::open_sam(&path).unwrap();

        let mut names = Vec::new();
        while let Some(seq) = reader.next_sequence().unwrap() {
            names.push(seq.header);
        }
        let _ = std::fs::remove_file(&path);

        assert_eq!(names, vec!["@read2".to_string()]);
    }
}
