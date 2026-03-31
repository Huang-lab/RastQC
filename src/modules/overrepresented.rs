use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{QCModule, QCResult};
use std::collections::HashMap;

const OBSERVATION_CUTOFF: usize = 100_000;

pub struct OverrepresentedSeqs {
    sequences: HashMap<String, u64>,
    total_count: u64,
    dup_length: usize,
    reached_limit: bool,
    count_at_limit: u64,
    // Results
    overrep_entries: Vec<OverrepEntry>,
    qc_result: QCResult,
}

pub struct OverrepEntry {
    pub sequence: String,
    pub count: u64,
    pub percentage: f64,
    pub possible_source: String,
}

impl OverrepresentedSeqs {
    pub fn new(dup_length: usize) -> Self {
        OverrepresentedSeqs {
            sequences: HashMap::new(),
            total_count: 0,
            dup_length,
            reached_limit: false,
            count_at_limit: 0,
            overrep_entries: Vec::new(),
            qc_result: QCResult::NotRun,
        }
    }

    /// Expose sequence counts for the Duplication module (mirrors FastQC architecture)
    pub fn sequence_counts(&self) -> &HashMap<String, u64> {
        &self.sequences
    }

    pub fn count_at_unique_limit(&self) -> u64 {
        self.count_at_limit
    }

    pub fn total_sequence_count(&self) -> u64 {
        self.total_count
    }

    fn find_contaminant(seq: &str, config: &FastQCConfig) -> String {
        let seq_upper = seq.to_uppercase();
        let mut best_name = String::from("No Hit");
        let mut best_len = 0usize;

        for contaminant in &config.contaminants {
            // Check forward
            if let Some(match_len) = find_match(&seq_upper, &contaminant.sequence) {
                if match_len > best_len {
                    best_len = match_len;
                    best_name = format!("{} ({}% over {}bp)", contaminant.name, 100, match_len);
                }
            }
            // Check reverse complement
            if let Some(match_len) = find_match(&seq_upper, &contaminant.reverse_complement) {
                if match_len > best_len {
                    best_len = match_len;
                    best_name = format!(
                        "{} ({}% over {}bp) - revcomp",
                        contaminant.name, 100, match_len
                    );
                }
            }
        }

        best_name
    }
}

fn find_match(query: &str, target: &str) -> Option<usize> {
    let qb = query.as_bytes();
    let tb = target.as_bytes();

    // For short queries (<= 20bp), require exact substring match
    if qb.len() <= 20 {
        if target.contains(query) || query.contains(target) {
            return Some(qb.len().min(tb.len()));
        }
        return None;
    }

    // For longer queries: sliding window with 1 mismatch tolerance
    let min_match = 20;
    let mut best_match = 0usize;

    for offset in -(tb.len() as i64 - min_match as i64)..=(qb.len() as i64 - min_match as i64) {
        let mut matches = 0usize;
        let mut mismatches = 0usize;
        let mut match_len = 0usize;
        let mut best_in_window = 0usize;

        let q_start = offset.max(0) as usize;
        let t_start = (-offset).max(0) as usize;

        let mut qi = q_start;
        let mut ti = t_start;

        while qi < qb.len() && ti < tb.len() {
            if qb[qi] == tb[ti] {
                matches += 1;
                match_len += 1;
            } else {
                mismatches += 1;
                match_len += 1;
                if mismatches > 1 {
                    if match_len - mismatches >= min_match {
                        best_in_window = best_in_window.max(match_len - mismatches);
                    }
                    matches = 0;
                    mismatches = 0;
                    match_len = 0;
                }
            }
            qi += 1;
            ti += 1;
        }

        if match_len >= min_match {
            best_in_window = best_in_window.max(matches);
        }

        best_match = best_match.max(best_in_window);
    }

    if best_match >= min_match {
        Some(best_match)
    } else {
        None
    }
}

impl QCModule for OverrepresentedSeqs {
    fn name(&self) -> &str {
        "Overrepresented sequences"
    }

    fn key(&self) -> &str {
        "overrepresented"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_count += 1;

        let end = seq.sequence.len().min(self.dup_length);
        let trunc = String::from_utf8_lossy(&seq.sequence[..end]).to_uppercase();

        if self.reached_limit {
            if let Some(count) = self.sequences.get_mut(&trunc) {
                *count += 1;
            }
        } else {
            *self.sequences.entry(trunc).or_insert(0) += 1;

            if self.sequences.len() >= OBSERVATION_CUTOFF {
                self.reached_limit = true;
                self.count_at_limit = self.total_count;
            }
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.total_count == 0 {
            self.qc_result = QCResult::Pass;
            return;
        }

        let warn_thresh = config
            .get_limit("overrepresented")
            .map(|l| l.warn)
            .unwrap_or(0.1);
        let error_thresh = config
            .get_limit("overrepresented")
            .map(|l| l.error)
            .unwrap_or(1.0);

        // Always use total count as denominator (matching FastQC behavior)
        let effective_total = self.total_count as f64;

        // Find sequences above warning threshold
        let mut entries: Vec<(&String, &u64)> = self
            .sequences
            .iter()
            .filter(|(_, &count)| {
                let pct = count as f64 / effective_total * 100.0;
                pct >= warn_thresh
            })
            .collect();

        entries.sort_by(|a, b| b.1.cmp(a.1));

        self.qc_result = QCResult::Pass;

        for (seq, &count) in entries {
            let percentage = count as f64 / effective_total * 100.0;
            let source = Self::find_contaminant(seq, config);

            if percentage > error_thresh {
                self.qc_result = QCResult::Fail;
            } else if percentage > warn_thresh && self.qc_result != QCResult::Fail {
                self.qc_result = QCResult::Warn;
            }

            self.overrep_entries.push(OverrepEntry {
                sequence: seq.clone(),
                count,
                percentage,
                possible_source: source,
            });
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn has_chart(&self) -> bool {
        false
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Overrepresented sequences\t{}\n\
             #Sequence\tCount\tPercentage\tPossible Source\n",
            self.result().label()
        );

        for entry in &self.overrep_entries {
            out.push_str(&format!(
                "{}\t{}\t{:.4}\t{}\n",
                entry.sequence, entry.count, entry.percentage, entry.possible_source
            ));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        // This module uses a table, not a chart
        String::new()
    }
}
