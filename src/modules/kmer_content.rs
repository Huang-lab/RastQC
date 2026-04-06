use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{BaseGroup, QCModule, QCResult};
use std::any::Any;
use std::collections::HashMap;

const MAX_SEQ_LENGTH: usize = 500;
const SAMPLING_RATE: u64 = 50; // process 1 in 50

struct KmerInfo {
    count: u64,
    positions: Vec<u64>,
}

pub struct KmerContent {
    kmer_size: usize,
    kmers: HashMap<Vec<u8>, KmerInfo>,
    total_kmer_counts: Vec<u64>,
    total_sequences: u64,
    max_length: usize,
    /// Reusable buffer for uppercase sequence
    upper_buf: Vec<u8>,
    // Results
    groups: Vec<BaseGroup>,
    result_entries: Vec<KmerResult>,
    min_pvalue_log10: f64,
    qc_result: QCResult,
}

pub struct KmerResult {
    pub sequence: String,
    pub count: u64,
    #[allow(dead_code)]
    pub obs_exp_overall: f64,
    pub obs_exp_max: f64,
    pub max_position: String,
    pub pvalue_log10: f64,
}

impl KmerContent {
    pub fn new(kmer_size: usize) -> Self {
        KmerContent {
            kmer_size,
            kmers: HashMap::new(),
            total_kmer_counts: Vec::new(),
            total_sequences: 0,
            max_length: 0,
            upper_buf: Vec::with_capacity(MAX_SEQ_LENGTH),
            groups: Vec::new(),
            result_entries: Vec::new(),
            min_pvalue_log10: 0.0,
            qc_result: QCResult::NotRun,
        }
    }
}

impl QCModule for KmerContent {
    fn name(&self) -> &str {
        "Kmer Content"
    }

    fn key(&self) -> &str {
        "kmer"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_sequences += 1;

        // Sample 1 in SAMPLING_RATE sequences
        if self.total_sequences % SAMPLING_RATE != 0 {
            return;
        }

        let len = seq.sequence.len().min(MAX_SEQ_LENGTH);
        if len < self.kmer_size {
            return;
        }
        // Cap tracked positions to avoid memory blowup on long reads
        let len = len.min(1000);
        if len > self.max_length {
            self.max_length = len;
        }

        // Ensure total_kmer_counts is long enough
        while self.total_kmer_counts.len() < len {
            self.total_kmer_counts.push(0);
        }

        // Reuse uppercase buffer (zero allocation in steady state)
        self.upper_buf.clear();
        self.upper_buf.extend(seq.sequence[..len].iter().map(|b| b.to_ascii_uppercase()));

        for i in 0..=(len - self.kmer_size) {
            let kmer = &self.upper_buf[i..i + self.kmer_size];

            // Skip kmers with N
            if kmer.contains(&b'N') {
                continue;
            }

            self.total_kmer_counts[i] += 1;

            // Use Vec<u8> key directly — avoids String allocation per kmer
            let entry = self.kmers.entry(kmer.to_vec()).or_insert_with(|| KmerInfo {
                count: 0,
                positions: Vec::new(),
            });
            entry.count += 1;

            while entry.positions.len() <= i {
                entry.positions.push(0);
            }
            entry.positions[i] += 1;
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.kmers.is_empty() || self.max_length == 0 {
            self.qc_result = QCResult::Pass;
            return;
        }

        self.groups = BaseGroup::make_groups(self.max_length);
        let _ng = self.groups.len();

        let total_possible: f64 = 4.0_f64.powi(self.kmer_size as i32);
        let total_kmers: u64 = self.total_kmer_counts.iter().sum();

        let warn = config.get_limit("kmer").map(|l| l.warn).unwrap_or(2.0);
        let error = config.get_limit("kmer").map(|l| l.error).unwrap_or(5.0);

        let mut results: Vec<KmerResult> = Vec::new();

        for (kmer_seq, info) in &self.kmers {
            if total_kmers == 0 {
                continue;
            }

            let expected_prop = info.count as f64 / total_kmers as f64;
            let obs_exp_overall = expected_prop * total_possible;

            // Find max obs/exp per group
            let mut max_obs_exp = 0.0_f64;
            let mut max_group_label = String::new();
            let mut min_pvalue = 1.0_f64;

            for (_gi, group) in self.groups.iter().enumerate() {
                let mut obs = 0u64;
                let mut total_in_group = 0u64;

                for pos in group.start..=group.end {
                    if pos < info.positions.len() {
                        obs += info.positions[pos];
                    }
                    if pos < self.total_kmer_counts.len() {
                        total_in_group += self.total_kmer_counts[pos];
                    }
                }

                if total_in_group > 0 && expected_prop > 0.0 {
                    let observed_prop = obs as f64 / total_in_group as f64;
                    let obs_exp = observed_prop / expected_prop;

                    if obs_exp > max_obs_exp {
                        max_obs_exp = obs_exp;
                        max_group_label = group.label();
                    }

                    // Approximate p-value using normal approximation to binomial
                    let n = total_in_group as f64;
                    let p = expected_prop;
                    let mean = n * p;
                    let std = (n * p * (1.0 - p)).sqrt();
                    if std > 0.0 {
                        let z = (obs as f64 - mean) / std;
                        // Approximate one-sided p-value
                        let pvalue = 0.5 * (-z * 0.7071).exp(); // rough approximation
                        let corrected = pvalue * total_possible;
                        if corrected < min_pvalue {
                            min_pvalue = corrected;
                        }
                    }
                }
            }

            let pvalue_log10 = if min_pvalue > 0.0 && min_pvalue < 1.0 {
                -min_pvalue.log10()
            } else {
                0.0
            };

            // Only keep significant & enriched kmers
            if pvalue_log10 > 1.0 && max_obs_exp > 5.0 {
                results.push(KmerResult {
                    sequence: String::from_utf8_lossy(kmer_seq).to_string(),
                    count: info.count,
                    obs_exp_overall,
                    obs_exp_max: max_obs_exp,
                    max_position: max_group_label,
                    pvalue_log10,
                });
            }
        }

        // Sort by p-value (most significant first)
        results.sort_by(|a, b| b.pvalue_log10.partial_cmp(&a.pvalue_log10).unwrap());
        results.truncate(20);

        self.min_pvalue_log10 = results.first().map(|r| r.pvalue_log10).unwrap_or(0.0);
        self.result_entries = results;

        self.qc_result = if self.min_pvalue_log10 > error {
            QCResult::Fail
        } else if self.min_pvalue_log10 > warn {
            QCResult::Warn
        } else {
            QCResult::Pass
        };
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn has_chart(&self) -> bool {
        false
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Kmer Content\t{}\n\
             #Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position\n",
            self.result().label()
        );

        for entry in &self.result_entries {
            out.push_str(&format!(
                "{}\t{}\t{:.6}\t{:.2}\t{}\n",
                entry.sequence,
                entry.count,
                10.0_f64.powf(-entry.pvalue_log10),
                entry.obs_exp_max,
                entry.max_position
            ));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        String::new()
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn merge_from(&mut self, other: &mut dyn QCModule) {
        if let Some(other) = other.as_any_mut().downcast_mut::<Self>() {
            self.total_sequences += other.total_sequences;
            if other.max_length > self.max_length {
                self.max_length = other.max_length;
            }
            // Merge total_kmer_counts
            while self.total_kmer_counts.len() < other.total_kmer_counts.len() {
                self.total_kmer_counts.push(0);
            }
            for (i, &val) in other.total_kmer_counts.iter().enumerate() {
                self.total_kmer_counts[i] += val;
            }
            // Merge per-kmer counts and positions
            for (kmer_str, other_info) in other.kmers.drain() {
                let entry = self.kmers.entry(kmer_str).or_insert_with(|| KmerInfo {
                    count: 0,
                    positions: Vec::new(),
                });
                entry.count += other_info.count;
                while entry.positions.len() < other_info.positions.len() {
                    entry.positions.push(0);
                }
                for (i, &val) in other_info.positions.iter().enumerate() {
                    entry.positions[i] += val;
                }
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
