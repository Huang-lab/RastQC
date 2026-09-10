use super::fasthash::FxBuildHasher;
use super::{format_pct_label, QCModule, QCResult};
use crate::config::FastQCConfig;
use crate::io::Sequence;
use std::any::Any;
use std::collections::HashMap;

/// Cap on distinct sequences tracked before "freezing" (matching FastQC's own
/// bounded-memory approach). Note for the streaming-parallel path: this cap
/// is per-worker, not global — each worker thread freezes independently
/// based only on the distinct sequences *it* has seen, so the specific point
/// at which tracking freezes (and thus `count_at_unique_limit`, which feeds
/// `get_corrected_count`'s binomial correction) can differ between a
/// sequential run and a parallel run of the same file, and between the
/// individual workers of a single parallel run. In practice this cap is
/// rarely reached (100k distinct sequences is a lot), so it mostly matters
/// for extremely diverse/large inputs; see the identical caveat on
/// `KmerContent`'s sampling for the same class of intra-file-parallelism
/// tradeoff.
use crate::config::OBSERVATION_BUDGET;

/// Duplication level bin labels matching FastQC (16 bins, 0-indexed).
/// Slot 9 covers counts 10+ (up to 50), so its label is ">10" not "10".
const BIN_LABELS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k",
];

pub struct DuplicationLevel {
    /// Uppercase sequence bytes -> count (avoids String allocation per read)
    sequences: HashMap<Vec<u8>, u64, FxBuildHasher>,
    total_sequences: u64,
    count_at_unique_limit: u64,
    dup_length: usize,
    /// This instance's share of [`OBSERVATION_BUDGET`]; see
    /// `FastQCConfig::observation_cutoff`.
    observation_cutoff: usize,
    frozen: bool,
    unique_count: usize,
    /// Reusable buffer for uppercase conversion
    upper_buf: Vec<u8>,
    // Results
    total_percentages: [f64; 16],
    dedup_percentages: [f64; 16],
    percent_different: f64,
    qc_result: QCResult,
}

impl DuplicationLevel {
    pub fn new(dup_length: usize, observation_cutoff: usize) -> Self {
        DuplicationLevel {
            sequences: HashMap::default(),
            total_sequences: 0,
            count_at_unique_limit: 0,
            dup_length,
            observation_cutoff,
            frozen: false,
            unique_count: 0,
            upper_buf: Vec::with_capacity(dup_length),
            total_percentages: [0.0; 16],
            dedup_percentages: [0.0; 16],
            percent_different: 100.0,
            qc_result: QCResult::NotRun,
        }
    }

    /// Map a duplication count to a bin index (0-15), matching FastQC's binning.
    fn dup_slot(count: u64) -> usize {
        let c = count.saturating_sub(1); // FastQC uses tempDupSlot = dupLevel - 1
        if c > 9999 || count == 0 {
            15
        } else if c > 4999 {
            14
        } else if c > 999 {
            13
        } else if c > 499 {
            12
        } else if c > 99 {
            11
        } else if c > 49 {
            10
        } else if c > 9 {
            9
        } else {
            c as usize
        }
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
        let limit_of_caring = 1.0 - (num_observations as f64 / (num_observations as f64 + 0.01));

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

    fn wants_all_reads(&self) -> bool {
        true
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_sequences += 1;

        // Truncate and uppercase into reusable buffer (zero allocation in steady state)
        let end = seq.sequence.len().min(self.dup_length);
        self.upper_buf.clear();
        self.upper_buf
            .extend(seq.sequence[..end].iter().map(|b| b.to_ascii_uppercase()));

        if self.frozen {
            // Only increment existing entries after freeze
            if let Some(count) = self.sequences.get_mut(&self.upper_buf) {
                *count += 1;
            }
            return;
        }

        // `HashMap::entry` needs an owned key up front, so
        // `entry(self.upper_buf.clone()).or_insert(0)` used to clone on
        // *every* call here, including the common case where this exact
        // sequence (post-truncation) has already been seen many times —
        // e.g. any read that's part of a duplicated/overrepresented cluster.
        // Look up by reference first and only pay for the clone on the
        // genuinely-new-sequence path.
        match self.sequences.get_mut(&self.upper_buf) {
            Some(count) => *count += 1,
            None => {
                self.sequences.insert(self.upper_buf.clone(), 1);
                self.unique_count += 1;
            }
        }
        self.count_at_unique_limit = self.total_sequences;

        if self.unique_count >= self.observation_cutoff {
            self.frozen = true;
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.sequences.is_empty() {
            return;
        }

        // Collate: count how many unique sequences have each duplication level
        let mut collated: HashMap<u64, u64> = HashMap::new();
        for &count in self.sequences.values() {
            *collated.entry(count).or_insert(0) += 1;
        }

        // Apply correction and accumulate into bins (matching FastQC exactly).
        // - dedup_percentages: percentage of *deduplicated* library per bin (# unique seqs)
        // - total_percentages: percentage of *total* library per bin (# reads, weighted by count)
        let mut dedup_total: f64 = 0.0;
        let mut raw_total: f64 = 0.0;

        self.total_percentages = [0.0; 16];
        self.dedup_percentages = [0.0; 16];

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
            // Total = weighted by how many reads this level contributes
            self.total_percentages[slot] += corrected * dup_level as f64;
            // Dedup = weighted only by how many unique sequences are at this level
            self.dedup_percentages[slot] += corrected;
        }

        // Convert to percentages
        if raw_total > 0.0 {
            for i in 0..16 {
                self.total_percentages[i] = self.total_percentages[i] / raw_total * 100.0;
            }
        }
        if dedup_total > 0.0 {
            for i in 0..16 {
                self.dedup_percentages[i] = self.dedup_percentages[i] / dedup_total * 100.0;
            }
        }

        self.percent_different = if raw_total > 0.0 {
            (dedup_total / raw_total) * 100.0
        } else {
            100.0
        };

        let warn = config
            .get_limit("duplication")
            .map(|l| l.warn)
            .unwrap_or(70.0);
        let error = config
            .get_limit("duplication")
            .map(|l| l.error)
            .unwrap_or(50.0);

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

        for (i, label) in BIN_LABELS.iter().enumerate() {
            out.push_str(&format!(
                "{}\t{:.2}\t{:.2}\n",
                label, self.dedup_percentages[i], self.total_percentages[i]
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
            if i > 0 {
                svg.push(' ');
            }
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
            mt + ph,
            ml + pw,
            mt + ph
        ));

        // Y-axis ticks
        let y_steps = 5;
        for i in 0..=y_steps {
            let frac = i as f64 / y_steps as f64;
            let val = max_pct * frac;
            let y = mt + ph * (1.0 - frac);
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" font-size="10">{}</text>"##,
                ml - 5.0, y, format_pct_label(val)
            ));
        }

        // X labels
        for (i, label) in BIN_LABELS.iter().enumerate().take(n) {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{label}</text>"##
            ));
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Sequence Duplication Levels [{:.2}% remaining after dedup]</text>"##,
            width / 2.0, self.percent_different
        ));

        // Y axis label
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">% of sequences</text>"##,
            mt + ph / 2.0,
            mt + ph / 2.0
        ));

        // X axis label
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" text-anchor="middle" font-size="11">Sequence Duplication Level</text>"##,
            ml + pw / 2.0,
            height - 5.0
        ));

        svg.push_str("</svg>");
        svg
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn merge_from(&mut self, other: &mut dyn QCModule) {
        if let Some(other) = other.as_any_mut().downcast_mut::<Self>() {
            for (seq, count) in other.sequences.drain() {
                *self.sequences.entry(seq).or_insert(0) += count;
            }
            self.total_sequences += other.total_sequences;
            self.count_at_unique_limit += other.count_at_unique_limit;
            self.unique_count = self.sequences.len();
            // The merged table stands for the whole file, so it is the
            // whole-file budget that decides whether tracking is frozen —
            // not the per-worker share each instance was built with.
            self.frozen = self.unique_count >= OBSERVATION_BUDGET;
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(bases: &[u8]) -> Sequence {
        Sequence {
            header: "@t".to_string(),
            sequence: bases.to_vec(),
            quality: vec![b'I'; bases.len()],
            filtered: false,
        }
    }

    #[test]
    fn repeated_sequences_are_counted_correctly_after_the_first_sighting() {
        // Regression test for the clone-avoidance refactor in
        // process_sequence: a sequence seen many times must still increment
        // the same HashMap entry every time, not just on first sight.
        let mut m = DuplicationLevel::new(50, OBSERVATION_BUDGET);
        for _ in 0..5 {
            m.process_sequence(&seq(b"ACGTACGTACGT"));
        }
        m.process_sequence(&seq(b"TTTTGGGGCCCC"));

        assert_eq!(m.total_sequences, 6);
        assert_eq!(m.unique_count, 2);
        assert_eq!(*m.sequences.get(b"ACGTACGTACGT".as_slice()).unwrap(), 5);
        assert_eq!(*m.sequences.get(b"TTTTGGGGCCCC".as_slice()).unwrap(), 1);
        assert_eq!(m.count_at_unique_limit, 6);
    }

    #[test]
    fn sequence_case_is_normalized_before_counting() {
        // Lowercase and uppercase forms of the same sequence must land in
        // the same bucket (the get_mut-then-insert refactor must not bypass
        // the existing to_ascii_uppercase() normalization).
        let mut m = DuplicationLevel::new(50, OBSERVATION_BUDGET);
        m.process_sequence(&seq(b"acgtACGT"));
        m.process_sequence(&seq(b"ACGTacgt"));

        assert_eq!(
            m.unique_count, 1,
            "case-insensitive dedup must merge both reads into one entry"
        );
        assert_eq!(*m.sequences.get(b"ACGTACGT".as_slice()).unwrap(), 2);
    }

    #[test]
    fn get_corrected_count_returns_observations_when_all_sequences_tracked() {
        // count_at_limit == total_count: every sequence was observed, so the
        // raw count needs no binomial correction.
        assert_eq!(DuplicationLevel::get_corrected_count(10, 10, 1, 7), 7.0);
    }

    #[test]
    fn get_corrected_count_returns_observations_when_untracked_pool_too_small() {
        // total_count - num_observations (4) < count_at_limit (5): FastQC
        // bails out to the raw count rather than attempting correction.
        assert_eq!(DuplicationLevel::get_corrected_count(5, 10, 1, 6), 6.0);
    }

    #[test]
    fn get_corrected_count_applies_binomial_correction() {
        // Hand-computable case: the per-step survival probabilities telescope
        // to exactly 9/10 * 8/9 * 7/8 * 6/7 * 5/6 = 5/10 = 0.5, so p_seeing =
        // 0.5 and the corrected count doubles the raw observation count.
        let corrected = DuplicationLevel::get_corrected_count(5, 10, 1, 5);
        assert!((corrected - 10.0).abs() < 1e-9, "got {corrected}");
    }

    #[test]
    fn dup_slot_matches_fastqc_binning() {
        assert_eq!(DuplicationLevel::dup_slot(1), 0);
        assert_eq!(DuplicationLevel::dup_slot(9), 8);
        assert_eq!(DuplicationLevel::dup_slot(10), 9);
        assert_eq!(DuplicationLevel::dup_slot(11), 9);
        assert_eq!(DuplicationLevel::dup_slot(51), 10);
        assert_eq!(DuplicationLevel::dup_slot(10_001), 15);
        assert_eq!(DuplicationLevel::dup_slot(0), 15);
    }
}
