use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{QCModule, QCResult};
use std::collections::HashMap;

const OBSERVATION_CUTOFF: usize = 100_000;

/// Duplication level bin labels matching FastQC (16 bins, 0-indexed)
const BIN_LABELS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k",
];

pub struct DuplicationLevel {
    /// Sequence string -> count (same tracking as OverRepresentedSeqs, matching FastQC)
    sequences: HashMap<String, u64>,
    total_sequences: u64,
    count_at_unique_limit: u64,
    dup_length: usize,
    frozen: bool,
    unique_count: usize,
    // Results
    total_percentages: [f64; 16],
    percent_different: f64,
    qc_result: QCResult,
}

impl DuplicationLevel {
    pub fn new(dup_length: usize) -> Self {
        DuplicationLevel {
            sequences: HashMap::new(),
            total_sequences: 0,
            count_at_unique_limit: 0,
            dup_length,
            frozen: false,
            unique_count: 0,
            total_percentages: [0.0; 16],
            percent_different: 100.0,
            qc_result: QCResult::NotRun,
        }
    }

    /// Map a duplication count to a bin index (0-15), matching FastQC's binning.
    fn dup_slot(count: u64) -> usize {
        let c = count.saturating_sub(1); // FastQC uses tempDupSlot = dupLevel - 1
        if c > 9999 || count == 0 { 15 }
        else if c > 4999 { 14 }
        else if c > 999 { 13 }
        else if c > 499 { 12 }
        else if c > 99 { 11 }
        else if c > 49 { 10 }
        else if c > 9 { 9 }
        else { c as usize }
    }

    /// FastQC's exact getCorrectedCount: iterative binomial correction.
    /// Given that we observed `num_observations` unique sequences at duplication
    /// level `dup_level`, and we only tracked the first `count_at_limit` of
    /// `total_count` total sequences, estimate how many unique sequences we
    /// would have seen if we tracked all of them.
    fn get_corrected_count(
        count_at_limit: u64,
        total_count: u64,
        dup_level: u64,
        num_observations: u64,
    ) -> f64 {
        // Early bailouts matching FastQC
        if count_at_limit == total_count {
            return num_observations as f64;
        }
        if total_count.saturating_sub(num_observations) < count_at_limit {
            return num_observations as f64;
        }

        // Probability of NOT seeing a sequence with this duplication level
        // within the first countAtLimit sequences
        let mut p_not_seeing: f64 = 1.0;

        // Limit below which correction is negligible (<0.01 of an observation)
        let limit_of_caring =
            1.0 - (num_observations as f64 / (num_observations as f64 + 0.01));

        for i in 0..count_at_limit {
            let numerator = (total_count - i).saturating_sub(dup_level) as f64;
            let denominator = (total_count - i) as f64;
            if denominator <= 0.0 {
                break;
            }
            p_not_seeing *= numerator / denominator;

            if p_not_seeing < limit_of_caring {
                p_not_seeing = 0.0;
                break;
            }
        }

        let p_seeing = 1.0 - p_not_seeing;
        if p_seeing <= 0.0 {
            return num_observations as f64;
        }
        num_observations as f64 / p_seeing
    }
}

impl QCModule for DuplicationLevel {
    fn name(&self) -> &str {
        "Sequence Duplication Levels"
    }

    fn key(&self) -> &str {
        "duplication"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_sequences += 1;

        // Truncate for duplication detection (matching FastQC/OverRepresentedSeqs)
        let end = seq.sequence.len().min(self.dup_length);
        let trunc = String::from_utf8_lossy(&seq.sequence[..end]).to_uppercase();

        if self.frozen {
            // Only increment existing entries after freeze
            if let Some(count) = self.sequences.get_mut(&trunc) {
                *count += 1;
            }
        } else {
            let entry = self.sequences.entry(trunc).or_insert(0);
            *entry += 1;
            if *entry == 1 {
                self.unique_count += 1;
            }
            self.count_at_unique_limit = self.total_sequences;

            if self.unique_count >= OBSERVATION_CUTOFF {
                self.frozen = true;
            }
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.sequences.is_empty() {
            return;
        }

        // Collate: count how many unique sequences have each duplication level
        let mut collated: HashMap<u64, u64> = HashMap::new();
        for (_, &count) in &self.sequences {
            *collated.entry(count).or_insert(0) += 1;
        }

        // Apply correction and accumulate into bins (matching FastQC exactly)
        let mut dedup_total: f64 = 0.0;
        let mut raw_total: f64 = 0.0;

        self.total_percentages = [0.0; 16];

        for (&dup_level, &num_observations) in &collated {
            let corrected = Self::get_corrected_count(
                self.count_at_unique_limit,
                self.total_sequences,
                dup_level,
                num_observations,
            );

            dedup_total += corrected;
            raw_total += corrected * dup_level as f64;

            let slot = Self::dup_slot(dup_level);
            self.total_percentages[slot] += corrected * dup_level as f64;
        }

        // Convert to percentages
        if raw_total > 0.0 {
            for i in 0..16 {
                self.total_percentages[i] = self.total_percentages[i] / raw_total * 100.0;
            }
        }

        self.percent_different = if raw_total > 0.0 {
            (dedup_total / raw_total) * 100.0
        } else {
            100.0
        };

        let warn = config.get_limit("duplication").map(|l| l.warn).unwrap_or(70.0);
        let error = config.get_limit("duplication").map(|l| l.error).unwrap_or(50.0);

        self.qc_result = if self.percent_different < error {
            QCResult::Fail
        } else if self.percent_different < warn {
            QCResult::Warn
        } else {
            QCResult::Pass
        };
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Sequence Duplication Levels\t{}\n\
             #Total Deduplicated Percentage\t{:.2}\n\
             #Duplication Level\tPercentage of deduplicated\tPercentage of total\n",
            self.result().label(),
            self.percent_different
        );

        for i in 0..16 {
            // FastQC only outputs totalPercentages (not separate dedup line)
            out.push_str(&format!(
                "{}\t{:.2}\t{:.2}\n",
                BIN_LABELS[i], 0.0, self.total_percentages[i]
            ));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        let width = 800.0_f64;
        let height = 400.0_f64;
        let ml = 60.0;
        let mr = 20.0;
        let mt = 30.0;
        let mb = 80.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;
        let n = 16;

        let max_pct = self
            .total_percentages
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
            .max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Total line (blue)
        svg.push_str(r##"<polyline points=""##);
        for i in 0..n {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph * (1.0 - self.total_percentages[i] / max_pct);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#0000ff" stroke-width="2" />"##);

        // Axes
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{}" stroke="black" />"##,
            mt + ph
        ));
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{}" x2="{}" y2="{}" stroke="black" />"##,
            mt + ph, ml + pw, mt + ph
        ));

        // X labels
        for i in 0..n {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{}</text>"##,
                BIN_LABELS[i]
            ));
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Sequence Duplication Levels [{:.2}% remaining after dedup]</text>"##,
            width / 2.0, self.percent_different
        ));

        svg.push_str("</svg>");
        svg
    }
}
