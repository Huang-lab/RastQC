/// POD5 reader for Oxford Nanopore data.
///
/// POD5 files use Apache Arrow IPC format to store sequencing data.
/// Note: POD5 files primarily contain raw signal data. Basecalled
/// sequences may be included as metadata columns in some configurations.
///
/// This module requires the `pod5` or `nanopore` feature to be enabled
/// at compile time.

#[cfg(feature = "arrow")]
use anyhow::{bail, Result};
#[cfg(feature = "arrow")]
use arrow::array::Array;
#[cfg(feature = "arrow")]
use arrow::ipc::reader::FileReader;
#[cfg(feature = "arrow")]
use std::fs::File;
#[cfg(feature = "arrow")]
use std::path::Path;

#[cfg(feature = "arrow")]
use super::Sequence;

#[cfg(feature = "arrow")]
pub struct Pod5Reader {
    sequences: Vec<Sequence>,
    index: usize,
}

#[cfg(feature = "arrow")]
impl Pod5Reader {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let reader = FileReader::try_new(file, None)?;

        let mut sequences = Vec::new();

        // POD5 files store data in Arrow record batches.
        // Look for columns containing basecalled sequence and quality data.
        for batch_result in reader {
            let batch = batch_result?;
            let schema = batch.schema();

            // Look for sequence and quality columns
            let seq_col = schema.fields().iter().position(|f| {
                let name = f.name().to_lowercase();
                name == "sequence" || name == "basecall" || name == "read_sequence"
            });

            let qual_col = schema.fields().iter().position(|f| {
                let name = f.name().to_lowercase();
                name == "quality" || name == "basecall_quality" || name == "read_quality"
            });

            let id_col = schema.fields().iter().position(|f| {
                let name = f.name().to_lowercase();
                name == "read_id" || name == "id"
            });

            if let (Some(si), Some(qi)) = (seq_col, qual_col) {
                let seq_array = batch
                    .column(si)
                    .as_any()
                    .downcast_ref::<arrow::array::StringArray>();
                let qual_array = batch
                    .column(qi)
                    .as_any()
                    .downcast_ref::<arrow::array::StringArray>();

                let id_array = id_col.and_then(|i| {
                    batch
                        .column(i)
                        .as_any()
                        .downcast_ref::<arrow::array::StringArray>()
                });

                if let (Some(seqs), Some(quals)) = (seq_array, qual_array) {
                    for row in 0..seqs.len() {
                        if seqs.is_null(row) || quals.is_null(row) {
                            continue;
                        }
                        let seq_str = seqs.value(row);
                        let qual_str = quals.value(row);
                        let read_id = id_array
                            .and_then(|a| {
                                if a.is_null(row) {
                                    None
                                } else {
                                    Some(a.value(row))
                                }
                            })
                            .unwrap_or("unknown");

                        sequences.push(Sequence {
                            header: format!("@{}", read_id),
                            sequence: seq_str.as_bytes().to_vec(),
                            quality: qual_str.as_bytes().to_vec(),
                            filtered: false,
                        });
                    }
                }
            }
        }

        if sequences.is_empty() {
            bail!(
                "No basecalled sequences found in POD5 file: {}. \
                 POD5 files often contain raw signal only. \
                 Use 'dorado basecaller' to produce basecalled FASTQ/BAM first.",
                path.display()
            );
        }

        Ok(Pod5Reader {
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

// Stub when feature is not enabled
#[cfg(not(feature = "arrow"))]
pub struct Pod5Reader;

#[cfg(not(feature = "arrow"))]
impl Pod5Reader {
    pub fn open(_path: &std::path::Path) -> anyhow::Result<Self> {
        anyhow::bail!(
            "POD5 support requires the 'pod5' or 'nanopore' feature. \
             Rebuild with: cargo build --features pod5"
        );
    }

    pub fn next_sequence(&mut self) -> anyhow::Result<Option<super::Sequence>> {
        Ok(None)
    }
}
