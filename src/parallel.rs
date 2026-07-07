use crate::config::FastQCConfig;
use crate::io::{Sequence, SequenceReader};
use crate::modules::{merge_module_sets, ModuleFactory, QCModule};
use anyhow::Result;
use crossbeam::channel;
use std::path::Path;
use std::thread;

/// Target bytes per batch (~4MB). Batch size is computed dynamically
/// from the first few reads' average length so long reads don't blow up memory.
const TARGET_BATCH_BYTES: usize = 4 * 1024 * 1024;

/// Fallback batch size when read lengths are unknown.
const DEFAULT_BATCH_SIZE: usize = 16_384;

/// Channel capacity: number of batches buffered in the channel.
const CHANNEL_CAPACITY: usize = 2;

/// Process a file using streaming parallelism.
///
/// Architecture:
///   Reader thread → bounded channel → N worker threads (each with own modules)
///   → merge all worker module states → final result
///
/// Unlike the old approach, this never buffers the entire file in memory.
/// Memory usage is bounded: O(BATCH_SIZE × CHANNEL_CAPACITY × avg_read_size).
pub fn process_file_parallel(
    path: &Path,
    config: &FastQCConfig,
    num_threads: usize,
) -> Result<(Vec<Box<dyn QCModule>>, u64)> {
    let num_workers = num_threads.max(1);

    // Bounded channel: reader sends batches, workers consume them
    let (sender, receiver) = channel::bounded::<Vec<Sequence>>(CHANNEL_CAPACITY);

    // Spawn reader thread
    let path_owned = path.to_path_buf();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut reader = SequenceReader::open(&path_owned)?;

        // Determine batch size from first reads' average length
        let mut probe_reads: Vec<Sequence> = Vec::new();
        let mut total_len: usize = 0;
        while probe_reads.len() < 100 {
            if let Some(seq) = reader.next_sequence()? {
                total_len += seq.len();
                probe_reads.push(seq);
            } else {
                break;
            }
        }

        let batch_size = if probe_reads.is_empty() {
            DEFAULT_BATCH_SIZE
        } else {
            let avg_len = (total_len / probe_reads.len()).max(1);
            // ~4MB per batch, minimum 64 reads, maximum 16K reads
            (TARGET_BATCH_BYTES / avg_len).clamp(64, DEFAULT_BATCH_SIZE)
        };

        // Send probe reads as the first batch
        let mut batch = Vec::with_capacity(batch_size);
        for seq in probe_reads {
            batch.push(seq);
            if batch.len() >= batch_size {
                if sender.send(batch).is_err() {
                    return Ok(());
                }
                batch = Vec::with_capacity(batch_size);
            }
        }

        // Continue with remaining reads
        while let Some(seq) = reader.next_sequence()? {
            batch.push(seq);
            if batch.len() >= batch_size {
                if sender.send(batch).is_err() {
                    return Ok(());
                }
                batch = Vec::with_capacity(batch_size);
            }
        }
        if !batch.is_empty() {
            let _ = sender.send(batch);
        }
        Ok(())
    });

    // Spawn worker threads, each with independent module instances
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
        let rx = receiver.clone();
        let worker_config = config.clone();
        let handle = thread::spawn(move || -> (Vec<Box<dyn QCModule>>, u64) {
            let mut modules = ModuleFactory::create_modules(&worker_config);
            let mut count: u64 = 0;

            while let Ok(batch) = rx.recv() {
                for seq in &batch {
                    for module in modules.iter_mut() {
                        module.process_sequence(seq);
                    }
                    count += 1;
                }
            }
            (modules, count)
        });
        worker_handles.push(handle);
    }

    // Drop our copy of the receiver so workers see channel close
    drop(receiver);

    // Wait for reader to finish
    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Reader thread panicked"))??;

    // Collect worker results
    let mut worker_results: Vec<(Vec<Box<dyn QCModule>>, u64)> = Vec::new();
    for handle in worker_handles {
        let result = handle
            .join()
            .map_err(|_| anyhow::anyhow!("Worker thread panicked"))?;
        worker_results.push(result);
    }

    // Merge all worker module states into the first worker's state
    let total_count: u64 = worker_results.iter().map(|(_, c)| c).sum();

    if worker_results.is_empty() {
        return Ok((ModuleFactory::create_modules(config), 0));
    }

    let (mut final_modules, _) = worker_results.remove(0);
    for (mut worker_modules, _) in worker_results {
        merge_module_sets(&mut final_modules, &mut worker_modules);
    }

    // Calculate final results on merged state
    for module in final_modules.iter_mut() {
        module.calculate_results(config);
    }

    Ok((final_modules, total_count))
}

/// Conservative estimate of gzip/bzip2 compression ratio for FASTQ base-call
/// data. Real ratios vary (~2x-10x+) with read complexity and codec; this is
/// a deliberately modest multiplier so files that don't compress as well
/// still cross the parallel threshold at a sane on-disk size.
const ASSUMED_COMPRESSION_RATIO: u64 = 4;

/// Estimate decompressed size from on-disk size and a lowercased filename.
///
/// FASTQ compresses well, so gating parallelism on the compressed on-disk
/// size keeps the common case — a `.fastq.gz` well under 50 MB on disk but
/// hundreds of MB decompressed — on the slower serial path, losing most of
/// the benefit parallel processing is meant to provide for exactly this
/// input type. Multi-member gzip (bgzip/pigz) streams make the gzip ISIZE
/// trailer unreliable as an exact size, so this uses a conservative
/// fixed-ratio estimate instead of decoding.
///
/// `lowercase_name` must already be lowercased by the caller, matching the
/// convention used for compression-format detection elsewhere (e.g.
/// `FastqReader::open`), so `sample.FASTQ.GZ` is still recognized.
fn estimate_decompressed_size(on_disk_size: u64, lowercase_name: &str) -> u64 {
    if lowercase_name.ends_with(".gz") || lowercase_name.ends_with(".bz2") {
        on_disk_size.saturating_mul(ASSUMED_COMPRESSION_RATIO)
    } else {
        on_disk_size
    }
}

/// Check if a file is large enough to benefit from parallel processing.
pub fn should_use_parallel(path: &Path) -> bool {
    let Ok(metadata) = path.metadata() else {
        return false;
    };
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_lowercase();
    estimate_decompressed_size(metadata.len(), &name) > 50 * 1024 * 1024
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::FastQCConfig;
    use std::io::Write;

    #[test]
    fn estimate_decompressed_size_scales_compressed_extensions() {
        // A 13 MB .fastq.gz is estimated at 52 MB decompressed, crossing the
        // 50 MB parallel threshold even though the on-disk file is small.
        assert_eq!(
            estimate_decompressed_size(13 * 1024 * 1024, "sample.fastq.gz"),
            52 * 1024 * 1024
        );
        assert_eq!(
            estimate_decompressed_size(13 * 1024 * 1024, "sample.fastq.bz2"),
            52 * 1024 * 1024
        );
        // Uncompressed FASTQ is taken at face value.
        assert_eq!(
            estimate_decompressed_size(13 * 1024 * 1024, "sample.fastq"),
            13 * 1024 * 1024
        );
        assert_eq!(
            estimate_decompressed_size(13 * 1024 * 1024, ""),
            13 * 1024 * 1024
        );
    }

    #[test]
    fn estimate_decompressed_size_requires_caller_to_lowercase() {
        // should_use_parallel lowercases the filename before calling this,
        // matching FastqReader::open's convention, so a real ".GZ"/".Gz"
        // file is still recognized as compressed (case handled by the
        // caller, not here — this documents/locks in that division of work).
        assert_eq!(
            estimate_decompressed_size(13 * 1024 * 1024, "sample.fastq.gz"),
            52 * 1024 * 1024
        );
        assert_eq!(
            estimate_decompressed_size(13 * 1024 * 1024, "sample.fastq.GZ"),
            13 * 1024 * 1024
        );
    }

    #[test]
    fn estimate_decompressed_size_does_not_overflow_on_huge_files() {
        assert_eq!(
            estimate_decompressed_size(u64::MAX, "sample.fastq.gz"),
            u64::MAX
        );
    }

    /// Regression test for issue #3, parallel path: large files (>50 MB) are
    /// processed by `process_file_parallel`, whose reader thread decodes the
    /// gzip stream. A multi-member gzip stream must be fully decoded here too,
    /// not just on the serial path. We call the function directly so a tiny
    /// fixture exercises the same reader/decoder a 30 GB file would.
    #[test]
    fn test_parallel_reads_all_gzip_members() {
        // Build a 2-member gzip stream (mimics pigz/bgzip chunking).
        let mut stream = Vec::new();
        for chunk in 0..2 {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            for i in 0..50 {
                let record = format!("@r_{}_{}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n", chunk, i);
                encoder.write_all(record.as_bytes()).unwrap();
            }
            stream.extend_from_slice(&encoder.finish().unwrap());
        }

        let path = std::env::temp_dir().join("rastqc_parallel_multimember_test.fastq.gz");
        std::fs::write(&path, &stream).unwrap();

        let config = FastQCConfig::new(None, None, None, 7, false, 50).unwrap();
        let (_modules, count) = process_file_parallel(&path, &config, 4).unwrap();

        let _ = std::fs::remove_file(&path);

        assert_eq!(
            count, 100,
            "parallel path must read all 100 records across both gzip members, got {count}"
        );
    }

    #[test]
    fn should_use_parallel_recognizes_uppercase_gz_extension() {
        let dir = std::env::temp_dir().join(format!("rastqc_case_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // 13 MB on disk, named with an uppercase extension.
        let path = dir.join("sample.fastq.GZ");
        std::fs::write(&path, vec![0u8; 13 * 1024 * 1024]).unwrap();

        let result = should_use_parallel(&path);

        let _ = std::fs::remove_dir_all(&dir);
        assert!(
            result,
            "a 13 MB .GZ file should be estimated at ~52 MB decompressed and cross the parallel threshold"
        );
    }
}
