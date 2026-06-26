/// Fast5 reader for Oxford Nanopore data.
///
/// Fast5 files are HDF5 containers that can hold basecalled FASTQ data
/// in groups like `/Analyses/Basecall_1D_*/BaseCalled_template/Fastq`.
///
/// This module requires the `fast5` or `nanopore` feature to be enabled
/// at compile time, which pulls in the HDF5 C library dependency.

#[cfg(feature = "hdf5")]
use anyhow::{bail, Result};
#[cfg(feature = "hdf5")]
use hdf5::File as H5File;
#[cfg(feature = "hdf5")]
use std::path::Path;

#[cfg(feature = "hdf5")]
use super::Sequence;

#[cfg(feature = "hdf5")]
pub struct Fast5Reader {
    sequences: Vec<Sequence>,
    index: usize,
}

#[cfg(feature = "hdf5")]
impl Fast5Reader {
    pub fn open(path: &Path) -> Result<Self> {
        let file = H5File::open(path)?;
        let mut sequences = Vec::new();

        // Try multi-read Fast5 format first
        if let Ok(group) = file.group("/") {
            for read_name in group.member_names()? {
                if !read_name.starts_with("read_") {
                    continue;
                }

                // Look for basecalled FASTQ data
                let analyses_path = format!("{}/Analyses", read_name);
                if let Ok(analyses) = file.group(&analyses_path) {
                    if let Some(seq) =
                        extract_basecall_from_analyses(&file, &analyses_path, &analyses)?
                    {
                        sequences.push(seq);
                    }
                }
            }
        }

        // Try single-read Fast5 format
        if sequences.is_empty() {
            if let Ok(analyses) = file.group("/Analyses") {
                if let Some(seq) = extract_basecall_from_analyses(&file, "/Analyses", &analyses)? {
                    sequences.push(seq);
                }
            }
        }

        if sequences.is_empty() {
            bail!(
                "No basecalled sequences found in Fast5 file: {}. \
                 Ensure the file contains basecalled data (not raw signal only).",
                path.display()
            );
        }

        Ok(Fast5Reader {
            sequences,
            index: 0,
        })
    }

    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        if self.index >= self.sequences.len() {
            return Ok(None);
        }
        let seq = self.sequences[self.index].clone();
        self.index += 1;
        Ok(Some(seq))
    }
}

/// Extract basecalled FASTQ from an Analyses group.
/// Tries common basecaller output paths.
#[cfg(feature = "hdf5")]
fn extract_basecall_from_analyses(
    file: &H5File,
    analyses_path: &str,
    analyses: &hdf5::Group,
) -> Result<Option<Sequence>> {
    let basecall_groups: Vec<String> = analyses
        .member_names()?
        .into_iter()
        .filter(|n| n.starts_with("Basecall_1D"))
        .collect();

    // Sort to get latest version
    let mut sorted = basecall_groups;
    sorted.sort();

    for bc_name in sorted.iter().rev() {
        let fastq_path = format!("{}/{}/BaseCalled_template/Fastq", analyses_path, bc_name);

        if let Ok(dataset) = file.dataset(&fastq_path) {
            if let Ok(fastq_str) = dataset.read_scalar::<hdf5::types::VarLenUnicode>() {
                let fastq_str = fastq_str.as_str();
                if let Some(seq) = parse_fastq_string(fastq_str) {
                    return Ok(Some(seq));
                }
            }
        }
    }

    Ok(None)
}

/// Parse a FASTQ-formatted string into a Sequence.
#[cfg(feature = "hdf5")]
fn parse_fastq_string(fastq: &str) -> Option<Sequence> {
    let lines: Vec<&str> = fastq.lines().collect();
    if lines.len() < 4 {
        return None;
    }

    let header = lines[0].to_string();
    if !header.starts_with('@') {
        return None;
    }

    let filtered = header.contains(":Y:");

    Some(Sequence {
        header,
        sequence: lines[1].trim().as_bytes().to_vec(),
        quality: lines[3].trim().as_bytes().to_vec(),
        filtered,
    })
}

// Stub when feature is not enabled
#[cfg(not(feature = "hdf5"))]
pub struct Fast5Reader;

#[cfg(not(feature = "hdf5"))]
impl Fast5Reader {
    pub fn open(_path: &std::path::Path) -> anyhow::Result<Self> {
        anyhow::bail!(
            "Fast5 support requires the 'fast5' or 'nanopore' feature. \
             Rebuild with: cargo build --features fast5"
        );
    }

    pub fn next_sequence(&mut self) -> anyhow::Result<Option<super::Sequence>> {
        Ok(None)
    }
}
