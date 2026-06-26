use crate::config::FastQCConfig;
use crate::io::Sequence;
use crate::modules::{QCModule, QCResult};
use std::any::Any;

// Long-read QC metrics: N50, quality-stratified length distribution,
// and homopolymer error rates. Designed for PacBio HiFi and Oxford Nanopore data.

// ─── Read Length N50 & Statistics ────────────────────────────────────────────

pub struct ReadLengthN50 {
    lengths: Vec<u64>,
    total_bases: u64,
    n50: u64,
    n90: u64,
    mean_length: f64,
    median_length: u64,
    max_length: u64,
    min_length: u64,
    qc_result: QCResult,
}

impl ReadLengthN50 {
    pub fn new() -> Self {
        ReadLengthN50 {
            lengths: Vec::new(),
            total_bases: 0,
            n50: 0,
            n90: 0,
            mean_length: 0.0,
            median_length: 0,
            max_length: 0,
            min_length: u64::MAX,
            qc_result: QCResult::Pass,
        }
    }
}

impl QCModule for ReadLengthN50 {
    fn name(&self) -> &str {
        "Read Length N50"
    }

    fn key(&self) -> &str {
        "read_length_n50"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        let len = seq.len() as u64;
        self.lengths.push(len);
        self.total_bases += len;
        if len > self.max_length {
            self.max_length = len;
        }
        if len < self.min_length {
            self.min_length = len;
        }
    }

    fn calculate_results(&mut self, _config: &FastQCConfig) {
        if self.lengths.is_empty() {
            self.qc_result = QCResult::NotRun;
            return;
        }

        self.lengths.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending
        let total = self.total_bases;

        // Calculate N50 and N90
        let mut cumulative = 0u64;
        let mut n50_found = false;
        for &len in &self.lengths {
            cumulative += len;
            if !n50_found && cumulative >= total / 2 {
                self.n50 = len;
                n50_found = true;
            }
            if cumulative >= (total * 9) / 10 {
                self.n90 = len;
                break;
            }
        }

        let count = self.lengths.len() as f64;
        self.mean_length = total as f64 / count;
        self.median_length = self.lengths[self.lengths.len() / 2];

        if self.min_length == u64::MAX {
            self.min_length = 0;
        }

        self.qc_result = QCResult::Pass;
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn has_chart(&self) -> bool {
        false
    }

    fn text_data(&self) -> String {
        let mut text = String::new();
        text.push_str(&format!(">>Read Length N50\t{}\n", self.qc_result.label()));
        text.push_str("#Metric\tValue\n");
        text.push_str(&format!("Total Reads\t{}\n", self.lengths.len()));
        text.push_str(&format!("Total Bases\t{}\n", self.total_bases));
        text.push_str(&format!("N50\t{}\n", self.n50));
        text.push_str(&format!("N90\t{}\n", self.n90));
        text.push_str(&format!("Mean Length\t{:.1}\n", self.mean_length));
        text.push_str(&format!("Median Length\t{}\n", self.median_length));
        text.push_str(&format!("Min Length\t{}\n", self.min_length));
        text.push_str(&format!("Max Length\t{}\n", self.max_length));
        text.push_str(">>END_MODULE\n");
        text
    }

    fn svg_chart(&self) -> String {
        String::new()
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn merge_from(&mut self, other: &mut dyn QCModule) {
        if let Some(other) = other.as_any_mut().downcast_mut::<Self>() {
            self.lengths.append(&mut other.lengths);
            self.total_bases += other.total_bases;
            self.max_length = self.max_length.max(other.max_length);
            self.min_length = self.min_length.min(other.min_length);
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}

// ─── Quality-Stratified Length Distribution ──────────────────────────────────

/// Bins reads by quality and tracks length distribution within each quality tier.
pub struct QualityStratifiedLength {
    /// Quality tiers: Q<10, Q10-19, Q20-29, Q30-39, Q40+
    tier_counts: [u64; 5],
    tier_total_bases: [u64; 5],
    tier_mean_length: [f64; 5],
    qc_result: QCResult,
}

const TIER_LABELS: [&str; 5] = ["Q<10", "Q10-19", "Q20-29", "Q30-39", "Q40+"];

impl QualityStratifiedLength {
    pub fn new() -> Self {
        QualityStratifiedLength {
            tier_counts: [0; 5],
            tier_total_bases: [0; 5],
            tier_mean_length: [0.0; 5],
            qc_result: QCResult::Pass,
        }
    }

    fn quality_tier(mean_q: f64) -> usize {
        if mean_q < 10.0 {
            0
        } else if mean_q < 20.0 {
            1
        } else if mean_q < 30.0 {
            2
        } else if mean_q < 40.0 {
            3
        } else {
            4
        }
    }
}

impl QCModule for QualityStratifiedLength {
    fn name(&self) -> &str {
        "Quality Stratified Length"
    }

    fn key(&self) -> &str {
        "quality_stratified_length"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if seq.quality.is_empty() {
            return;
        }
        let mean_q: f64 = seq
            .quality
            .iter()
            .map(|&q| (q.saturating_sub(33)) as f64)
            .sum::<f64>()
            / seq.quality.len() as f64;

        let len = seq.len() as u64;
        let tier = Self::quality_tier(mean_q);
        self.tier_counts[tier] += 1;
        self.tier_total_bases[tier] += len;
    }

    fn calculate_results(&mut self, _config: &FastQCConfig) {
        for i in 0..5 {
            if self.tier_counts[i] > 0 {
                self.tier_mean_length[i] =
                    self.tier_total_bases[i] as f64 / self.tier_counts[i] as f64;
            }
        }

        let total: u64 = self.tier_counts.iter().sum();
        if total == 0 {
            self.qc_result = QCResult::NotRun;
            return;
        }

        // Warn if >50% of reads are below Q20
        let low_q: u64 = self.tier_counts[0] + self.tier_counts[1];
        if low_q as f64 / total as f64 > 0.5 {
            self.qc_result = QCResult::Warn;
        } else {
            self.qc_result = QCResult::Pass;
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn has_chart(&self) -> bool {
        false
    }

    fn text_data(&self) -> String {
        let mut text = String::new();
        text.push_str(&format!(
            ">>Quality Stratified Length\t{}\n",
            self.qc_result.label()
        ));
        text.push_str("#Quality Tier\tRead Count\tTotal Bases\tMean Length\n");
        for (i, label) in TIER_LABELS.iter().enumerate() {
            text.push_str(&format!(
                "{}\t{}\t{}\t{:.1}\n",
                label, self.tier_counts[i], self.tier_total_bases[i], self.tier_mean_length[i]
            ));
        }
        text.push_str(">>END_MODULE\n");
        text
    }

    fn svg_chart(&self) -> String {
        String::new()
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn merge_from(&mut self, other: &mut dyn QCModule) {
        if let Some(other) = other.as_any_mut().downcast_mut::<Self>() {
            for i in 0..5 {
                self.tier_counts[i] += other.tier_counts[i];
                self.tier_total_bases[i] += other.tier_total_bases[i];
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}

// ─── Homopolymer Error Rates ────────────────────────────────────────────────

/// Tracks homopolymer runs in sequences and detects potential systematic errors.
pub struct HomopolymerErrors {
    /// Count of homopolymer runs by (base, length)
    /// Tracks runs of length 3+
    run_counts: [[u64; 20]; 4], // A=0, C=1, G=2, T=3; lengths 3-22 mapped to indices 0-19
    total_runs: u64,
    total_bases_in_runs: u64,
    total_sequences: u64,
    total_bases: u64,
    qc_result: QCResult,
}

impl HomopolymerErrors {
    pub fn new() -> Self {
        HomopolymerErrors {
            run_counts: [[0; 20]; 4],
            total_runs: 0,
            total_bases_in_runs: 0,
            total_sequences: 0,
            total_bases: 0,
            qc_result: QCResult::Pass,
        }
    }

    fn base_index(b: u8) -> Option<usize> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }
}

impl QCModule for HomopolymerErrors {
    fn name(&self) -> &str {
        "Homopolymer Content"
    }

    fn key(&self) -> &str {
        "homopolymer"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_sequences += 1;
        self.total_bases += seq.len() as u64;

        if seq.sequence.is_empty() {
            return;
        }

        let mut run_base = seq.sequence[0];
        let mut run_len: usize = 1;

        for &base in &seq.sequence[1..] {
            if base == run_base {
                run_len += 1;
            } else {
                // Record completed run if length >= 3
                if run_len >= 3 {
                    if let Some(bi) = Self::base_index(run_base) {
                        let idx = (run_len - 3).min(19);
                        self.run_counts[bi][idx] += 1;
                        self.total_runs += 1;
                        self.total_bases_in_runs += run_len as u64;
                    }
                }
                run_base = base;
                run_len = 1;
            }
        }
        // Handle last run
        if run_len >= 3 {
            if let Some(bi) = Self::base_index(run_base) {
                let idx = (run_len - 3).min(19);
                self.run_counts[bi][idx] += 1;
                self.total_runs += 1;
                self.total_bases_in_runs += run_len as u64;
            }
        }
    }

    fn calculate_results(&mut self, _config: &FastQCConfig) {
        if self.total_sequences == 0 {
            self.qc_result = QCResult::NotRun;
            return;
        }

        // Warn if homopolymer runs constitute >5% of total bases
        let hp_fraction = self.total_bases_in_runs as f64 / self.total_bases.max(1) as f64;
        if hp_fraction > 0.10 {
            self.qc_result = QCResult::Fail;
        } else if hp_fraction > 0.05 {
            self.qc_result = QCResult::Warn;
        } else {
            self.qc_result = QCResult::Pass;
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn has_chart(&self) -> bool {
        false
    }

    fn text_data(&self) -> String {
        let base_labels = ['A', 'C', 'G', 'T'];
        let mut text = String::new();
        text.push_str(&format!(
            ">>Homopolymer Content\t{}\n",
            self.qc_result.label()
        ));
        text.push_str("#Base\tRun Length\tCount\n");

        for (bi, &label) in base_labels.iter().enumerate() {
            for li in 0..20 {
                let count = self.run_counts[bi][li];
                if count > 0 {
                    let run_len = li + 3;
                    text.push_str(&format!("{}\t{}\t{}\n", label, run_len, count));
                }
            }
        }

        text.push_str(&format!(
            "#Total homopolymer runs (3+bp): {}\n",
            self.total_runs
        ));
        text.push_str(&format!(
            "#Total bases in homopolymer runs: {} ({:.2}%)\n",
            self.total_bases_in_runs,
            if self.total_bases > 0 {
                self.total_bases_in_runs as f64 / self.total_bases as f64 * 100.0
            } else {
                0.0
            }
        ));
        text.push_str(">>END_MODULE\n");
        text
    }

    fn svg_chart(&self) -> String {
        String::new()
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn merge_from(&mut self, other: &mut dyn QCModule) {
        if let Some(other) = other.as_any_mut().downcast_mut::<Self>() {
            for bi in 0..4 {
                for li in 0..20 {
                    self.run_counts[bi][li] += other.run_counts[bi][li];
                }
            }
            self.total_runs += other.total_runs;
            self.total_bases_in_runs += other.total_bases_in_runs;
            self.total_sequences += other.total_sequences;
            self.total_bases += other.total_bases;
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
