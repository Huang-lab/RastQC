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

/// Check if a file is large enough to benefit from parallel processing.
pub fn should_use_parallel(path: &Path) -> bool {
    path.metadata()
        .map(|m| m.len() > 50 * 1024 * 1024) // > 50 MB
        .unwrap_or(false)
}
