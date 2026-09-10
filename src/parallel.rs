use crate::config::FastQCConfig;
use crate::io::block::{self, BlockParseError, FastqBlockReader, RecordIter, BLOCK_SIZE};
use crate::io::fastq;
use crate::io::{Sequence, SequenceReader};
use crate::modules::{merge_module_sets, should_process, ModuleFactory, QCModule};
use anyhow::Result;
use crossbeam::channel;
use std::path::Path;
use std::sync::Arc;
use std::thread;

/// Target bytes per batch (~4MB) on the generic record-at-a-time path. Batch
/// size is computed dynamically from the first few reads' average length so
/// long reads don't blow up memory.
const TARGET_BATCH_BYTES: usize = 4 * 1024 * 1024;

/// Fallback batch size when read lengths are unknown.
const DEFAULT_BATCH_SIZE: usize = 16_384;

/// Channel capacity: number of batches buffered in the channel.
const CHANNEL_CAPACITY: usize = 2;

/// Process a file using streaming parallelism.
///
/// Plain FASTQ (optionally gzip/bzip2 compressed) takes the block fast path,
/// where the reader thread only decompresses and cuts the stream at record
/// boundaries and the workers do the parsing. Every other input format — BAM,
/// SAM, FASTA, Fast5, POD5, and SOLiD colorspace FASTQ — goes through the
/// generic path, whose reader thread produces owned `Sequence` records one at
/// a time.
///
/// Neither path ever buffers the whole file: memory stays bounded by the
/// in-flight blocks (or batches) plus each worker's own module state.
pub fn process_file_parallel(
    path: &Path,
    config: &FastQCConfig,
    num_threads: usize,
) -> Result<ModuleState> {
    let num_workers = num_threads.clamp(1, max_useful_workers(path));

    if fastq_block_path_applicable(path) {
        if let Some(result) = process_fastq_blocks(path, config, num_workers)? {
            return Ok(result);
        }
    }
    process_records_parallel(path, config, num_workers)
}

/// Worker threads past which a single file stops benefiting from more of them.
///
/// One reader thread feeds every worker, so its throughput is the ceiling —
/// and each extra worker also adds a whole set of module state that has to be
/// allocated, then merged at the end. Past the ceiling, wall time gets
/// *worse* while resident memory keeps climbing, so a run given a large `-t`
/// on a single file would otherwise pay several hundred extra MB for a
/// slowdown.
///
/// The ceiling depends on what the reader has to do. On a 4M-read NextSeq
/// file, gzipped input peaked at 4 workers (inflate-bound) and uncompressed
/// input at 8. Any leftover budget is better spent on other files, which
/// `split_thread_budget` in `main.rs` already prefers.
fn max_useful_workers(path: &Path) -> usize {
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_lowercase();
    if name.ends_with(".gz") || name.ends_with(".bz2") {
        4
    } else {
        8
    }
}

/// Whether `path` can take the FASTQ block fast path.
fn fastq_block_path_applicable(path: &Path) -> bool {
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_lowercase();

    let is_fastq = ["fastq", "fq"].iter().any(|ext| {
        name.ends_with(&format!(".{ext}"))
            || name.ends_with(&format!(".{ext}.gz"))
            || name.ends_with(&format!(".{ext}.bz2"))
    });
    if !is_fastq {
        return false;
    }

    // A colorspace record decodes to one more base than its quality line has
    // scores, which the fast path deliberately doesn't model. A probe that
    // errors also defers to the generic path, so the failure gets reported
    // there with its usual context rather than from inside this check.
    !matches!(block::first_record_is_colorspace(path), Ok(true) | Err(_))
}

/// The FASTQ fast path: reader thread decompresses and cuts at record
/// boundaries, workers parse and analyze.
///
/// Returns `Ok(None)` when the input doesn't fit the strict
/// four-lines-per-record structure the block parser requires, so the caller
/// can retry on the tolerant per-record reader instead of failing the run.
fn process_fastq_blocks(
    path: &Path,
    config: &FastQCConfig,
    num_workers: usize,
) -> Result<Option<ModuleState>> {
    // Each block goes to two consumers: whichever pool worker is free, and
    // the single in-order instance of the `wants_all_reads` modules. Sharing
    // one `Arc` rather than copying keeps that free; both channels are
    // bounded, so the slower consumer applies backpressure to the reader and
    // in-flight block memory stays at ~(num_workers + 2) * BLOCK_SIZE.
    let capacity = num_workers + 2;
    let (block_tx, block_rx) = channel::bounded::<Arc<Vec<u8>>>(capacity);
    let (serial_tx, serial_rx) = channel::bounded::<Arc<Vec<u8>>>(capacity);

    let path_owned = path.to_path_buf();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut reader = FastqBlockReader::new(fastq::open_decompressed(&path_owned)?);
        loop {
            let mut buf = Vec::with_capacity(BLOCK_SIZE);
            if !reader.next_block(&mut buf)? {
                return Ok(());
            }
            let block = Arc::new(buf);
            // Send to the in-order consumer first so it never falls behind
            // the pool by more than the channel depth.
            if serial_tx.send(Arc::clone(&block)).is_err() {
                return Ok(());
            }
            if block_tx.send(block).is_err() {
                return Ok(());
            }
        }
    });

    // The in-order consumer: one module set, every block, file order.
    let serial_config = config.clone();
    let serial_handle = thread::spawn(move || -> Result<ModuleState, BlockParseError> {
        let mut modules = ModuleFactory::create_modules(&serial_config);
        let mut seq = Sequence::empty();
        while let Ok(block) = serial_rx.recv() {
            for record in RecordIter::new(&block) {
                let record = record?;
                seq.refill(
                    block::header_str(record.header)?,
                    record.sequence,
                    record.quality,
                );
                for module in modules.iter_mut() {
                    if !module.wants_all_reads()
                        || !should_process(&seq, serial_config.nofilter, module.as_ref())
                    {
                        continue;
                    }
                    module.process_sequence(&seq);
                }
            }
        }
        // Read counts come from the pool workers; this instance would
        // double them.
        Ok((modules, 0))
    });

    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
        let rx = block_rx.clone();
        let worker_config = config.clone();
        let handle = thread::spawn(move || -> Result<ModuleState, BlockParseError> {
            let mut modules = ModuleFactory::create_modules(&worker_config);
            let mut count: u64 = 0;
            // One record reused for the whole block: the modules copy out
            // whatever they need, so nothing outlives an iteration.
            let mut seq = Sequence::empty();

            while let Ok(block) = rx.recv() {
                for record in RecordIter::new(&block) {
                    let record = record?;
                    seq.refill(
                        block::header_str(record.header)?,
                        record.sequence,
                        record.quality,
                    );
                    for module in modules.iter_mut() {
                        // Left to the in-order consumer, and left empty
                        // here so it costs no memory per worker.
                        if module.wants_all_reads()
                            || !should_process(&seq, worker_config.nofilter, module.as_ref())
                        {
                            continue;
                        }
                        module.process_sequence(&seq);
                    }
                    count += 1;
                }
            }
            Ok((modules, count))
        });
        worker_handles.push(handle);
    }

    // Drop our copy so the workers see the block channel close.
    drop(block_rx);

    let reader_result = reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Reader thread panicked"))?;

    let mut worker_results: Vec<ModuleState> = Vec::new();
    let mut parse_error: Option<BlockParseError> = None;
    for handle in worker_handles {
        match handle
            .join()
            .map_err(|_| anyhow::anyhow!("Worker thread panicked"))?
        {
            Ok(result) => worker_results.push(result),
            Err(e) => {
                if parse_error.is_none() {
                    parse_error = Some(e);
                }
            }
        }
    }
    match serial_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Serial module thread panicked"))?
    {
        Ok(result) => worker_results.push(result),
        Err(e) => {
            if parse_error.is_none() {
                parse_error = Some(e);
            }
        }
    }

    if let Some(e) = parse_error {
        if !config.quiet {
            eprintln!("Note: retrying with the per-record FASTQ reader ({e})");
        }
        return Ok(None);
    }
    // Only surface reader I/O errors once parsing is known to be clean, so a
    // structural mismatch is retried rather than reported as an I/O failure.
    reader_result?;

    finish_modules(worker_results, config).map(Some)
}

/// The generic path: reader thread produces owned `Sequence` records, workers
/// consume batches of them. Used for every non-FASTQ format.
fn process_records_parallel(
    path: &Path,
    config: &FastQCConfig,
    num_workers: usize,
) -> Result<ModuleState> {
    // Bounded channels: reader sends each batch to the worker pool and to the
    // single in-order consumer of the `wants_all_reads` modules, sharing one
    // `Arc` between them.
    let (sender, receiver) = channel::bounded::<Arc<Vec<Sequence>>>(CHANNEL_CAPACITY);
    let (serial_tx, serial_rx) = channel::bounded::<Arc<Vec<Sequence>>>(CHANNEL_CAPACITY);

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

        // Hand one batch to both consumers, in-order consumer first.
        let dispatch = |batch: Vec<Sequence>| -> bool {
            let batch = Arc::new(batch);
            serial_tx.send(Arc::clone(&batch)).is_ok() && sender.send(batch).is_ok()
        };

        // Send probe reads as the first batch
        let mut batch = Vec::with_capacity(batch_size);
        for seq in probe_reads {
            batch.push(seq);
            if batch.len() >= batch_size
                && !dispatch(std::mem::replace(
                    &mut batch,
                    Vec::with_capacity(batch_size),
                ))
            {
                return Ok(());
            }
        }

        // Continue with remaining reads
        while let Some(seq) = reader.next_sequence()? {
            batch.push(seq);
            if batch.len() >= batch_size
                && !dispatch(std::mem::replace(
                    &mut batch,
                    Vec::with_capacity(batch_size),
                ))
            {
                return Ok(());
            }
        }
        if !batch.is_empty() {
            dispatch(batch);
        }
        Ok(())
    });

    // The in-order consumer: one module set, every batch, file order.
    let serial_config = config.clone();
    let serial_handle = thread::spawn(move || -> ModuleState {
        let mut modules = ModuleFactory::create_modules(&serial_config);
        while let Ok(batch) = serial_rx.recv() {
            for seq in batch.iter() {
                for module in modules.iter_mut() {
                    if !module.wants_all_reads()
                        || !should_process(seq, serial_config.nofilter, module.as_ref())
                    {
                        continue;
                    }
                    module.process_sequence(seq);
                }
            }
        }
        // Read counts come from the pool workers; this instance would double
        // them.
        (modules, 0)
    });

    // Spawn worker threads, each with independent module instances
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
        let rx = receiver.clone();
        let worker_config = config.clone();
        let handle = thread::spawn(move || -> ModuleState {
            let mut modules = ModuleFactory::create_modules(&worker_config);
            let mut count: u64 = 0;

            while let Ok(batch) = rx.recv() {
                for seq in batch.iter() {
                    for module in modules.iter_mut() {
                        // Left to the in-order consumer, and left empty here
                        // so it costs no memory per worker.
                        if module.wants_all_reads()
                            || !should_process(seq, worker_config.nofilter, module.as_ref())
                        {
                            continue;
                        }
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
    let mut worker_results: Vec<ModuleState> = Vec::new();
    for handle in worker_handles {
        let result = handle
            .join()
            .map_err(|_| anyhow::anyhow!("Worker thread panicked"))?;
        worker_results.push(result);
    }
    worker_results.push(
        serial_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Serial module thread panicked"))?,
    );

    finish_modules(worker_results, config)
}

/// One consumer's finished module state plus the reads it counted.
type ModuleState = (Vec<Box<dyn QCModule>>, u64);

/// Merge every worker's module state into one set and compute final results.
fn finish_modules(
    mut worker_results: Vec<ModuleState>,
    config: &FastQCConfig,
) -> Result<ModuleState> {
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
    fn max_useful_workers_is_lower_for_compressed_input() {
        // Compressed input is inflate-bound on the single reader thread, so
        // it saturates with fewer workers than plain input does.
        assert_eq!(max_useful_workers(Path::new("s.fastq.gz")), 4);
        assert_eq!(max_useful_workers(Path::new("s.fastq.bz2")), 4);
        assert_eq!(max_useful_workers(Path::new("s.fastq")), 8);
        // Detection is case-insensitive, matching the compression sniffing
        // the readers themselves do.
        assert_eq!(max_useful_workers(Path::new("s.FASTQ.GZ")), 4);
    }

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

    /// A file whose last record has no trailing newline must not lose that
    /// record on the parallel path. `should_use_parallel` gates on a 50 MB
    /// estimate, so this calls the pipeline directly to exercise it with a
    /// tiny fixture.
    #[test]
    fn parallel_keeps_a_final_record_with_no_trailing_newline() {
        let dir = std::env::temp_dir().join(format!("rastqc_notrail_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("notrail.fastq");

        let mut data = String::new();
        for i in 0..500 {
            data.push_str(&format!("@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n"));
        }
        let trimmed = data.trim_end_matches('\n');
        std::fs::write(&path, trimmed).unwrap();

        let config = FastQCConfig::new(None, None, None, 7, false, 50).unwrap();
        let (_modules, count) = process_file_parallel(&path, &config, 4).unwrap();

        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(
            count, 500,
            "the final record has no trailing newline but is complete; got {count}"
        );
    }

    /// The same file read with and without intra-file parallelism must give
    /// the same read count — the property that the whole two-consumer design
    /// exists to preserve.
    #[test]
    fn parallel_and_sequential_agree_on_read_count() {
        let dir = std::env::temp_dir().join(format!("rastqc_agree_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("agree.fastq");

        let mut data = String::new();
        for i in 0..2000 {
            // Repeat a small pool of sequences so the duplication and
            // overrepresented tables actually see duplicates.
            let seq = ["ACGTACGTAC", "TTTTAAAACC", "GGCCGGCCGG"][i % 3];
            data.push_str(&format!("@r{i}\n{seq}\n+\nIIIIIIIIII\n"));
        }
        std::fs::write(&path, &data).unwrap();

        let config = FastQCConfig::new(None, None, None, 7, false, 50).unwrap();
        let one = process_file_parallel(&path, &config, 1).unwrap().1;
        let many = process_file_parallel(&path, &config, 8).unwrap().1;

        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(one, 2000);
        assert_eq!(many, 2000);
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
