use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{QCModule, QCResult};
use std::collections::HashMap;

/// GCModel: for a given read length, precomputes the percentage bin weights
/// for each possible GC count. This matches FastQC's GCModel exactly.
struct GCModel {
    /// models[gc_count] = vec of (percentage_bin, increment)
    models: Vec<Vec<(usize, f64)>>,
}

impl GCModel {
    fn new(read_length: usize) -> Self {
        let mut claiming_counts = [0i32; 101];

        // First pass: count how many GC positions claim each percentage bin
        for pos in 0..=read_length {
            let low_count = (pos as f64 - 0.5).max(0.0).min(read_length as f64);
            let high_count = (pos as f64 + 0.5).max(0.0).min(read_length as f64);

            let low_pct = ((low_count * 100.0) / read_length as f64).round() as usize;
            let high_pct = ((high_count * 100.0) / read_length as f64).round() as usize;
            let high_pct = high_pct.min(100);

            for p in low_pct..=high_pct {
                claiming_counts[p] += 1;
            }
        }

        // Second pass: build model values with weighted increments
        let mut models = Vec::with_capacity(read_length + 1);
        for pos in 0..=read_length {
            let low_count = (pos as f64 - 0.5).max(0.0).min(read_length as f64);
            let high_count = (pos as f64 + 0.5).max(0.0).min(read_length as f64);

            let low_pct = ((low_count * 100.0) / read_length as f64).round() as usize;
            let high_pct = ((high_count * 100.0) / read_length as f64).round() as usize;
            let high_pct = high_pct.min(100);

            let mut values = Vec::new();
            for p in low_pct..=high_pct {
                if claiming_counts[p] > 0 {
                    values.push((p, 1.0 / claiming_counts[p] as f64));
                }
            }
            models.push(values);
        }

        GCModel { models }
    }

    fn get_values(&self, gc_count: usize) -> &[(usize, f64)] {
        if gc_count < self.models.len() {
            &self.models[gc_count]
        } else {
            &[]
        }
    }
}

pub struct PerSequenceGC {
    /// Histogram of GC percentages (0..=100), using GCModel smoothing
    gc_distribution: [f64; 101],
    total_count: u64,
    /// Cache of GCModels by read length
    gc_models: HashMap<usize, GCModel>,
    // Results
    theoretical: [f64; 101],
    deviation_percent: f64,
    qc_result: QCResult,
}

impl PerSequenceGC {
    pub fn new() -> Self {
        PerSequenceGC {
            gc_distribution: [0.0; 101],
            total_count: 0,
            gc_models: HashMap::new(),
            theoretical: [0.0; 101],
            deviation_percent: 0.0,
            qc_result: QCResult::NotRun,
        }
    }
}

impl QCModule for PerSequenceGC {
    fn name(&self) -> &str {
        "Per sequence GC content"
    }

    fn key(&self) -> &str {
        "gc_sequence"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if seq.sequence.is_empty() {
            return;
        }

        // Count GC bases (matching FastQC: counts all bases including N for length)
        let mut gc_count = 0usize;
        let seq_bytes = &seq.sequence;
        let mut seq_len = seq_bytes.len();

        // Truncate long sequences (matching FastQC)
        if seq_len > 1000 {
            let truncated = (seq_len / 1000) * 1000;
            seq_len = truncated;
        } else if seq_len > 100 {
            let truncated = (seq_len / 100) * 100;
            seq_len = truncated;
        }

        for i in 0..seq_len {
            match seq_bytes[i] {
                b'G' | b'C' | b'g' | b'c' => gc_count += 1,
                _ => {}
            }
        }

        // Get or create GCModel for this read length
        if !self.gc_models.contains_key(&seq_len) {
            self.gc_models.insert(seq_len, GCModel::new(seq_len));
        }
        let model = &self.gc_models[&seq_len];

        // Distribute GC count across percentage bins using the model
        for &(pct, increment) in model.get_values(gc_count) {
            if pct <= 100 {
                self.gc_distribution[pct] += increment;
            }
        }

        self.total_count += 1;
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.total_count == 0 {
            return;
        }

        let total = self.total_count as f64;

        // Find firstMode (index of maximum count)
        let mut first_mode: usize = 0;
        let mut mode_count = 0.0_f64;
        for i in 0..=100 {
            if self.gc_distribution[i] > mode_count {
                mode_count = self.gc_distribution[i];
                first_mode = i;
            }
        }

        if mode_count == 0.0 {
            return;
        }

        // FastQC's bidirectional mode averaging:
        // Average indices where count > gcDistribution[firstMode] - gcDistribution[firstMode]/10
        // (i.e., within 90% of the peak height, walking outward from peak)
        let threshold = self.gc_distribution[first_mode]
            - (self.gc_distribution[first_mode] / 10.0);

        let mut mode_sum = first_mode as f64;
        let mut mode_duplicates = 1usize;

        // Walk right from firstMode
        let mut fell_off_top = true;
        for i in (first_mode + 1)..=100 {
            if self.gc_distribution[i] > threshold {
                mode_sum += i as f64;
                mode_duplicates += 1;
            } else {
                fell_off_top = false;
                break;
            }
        }

        // Walk left from firstMode-1
        let mut fell_off_bottom = true;
        if first_mode > 0 {
            for i in (0..first_mode).rev() {
                if self.gc_distribution[i] > threshold {
                    mode_sum += i as f64;
                    mode_duplicates += 1;
                } else {
                    fell_off_bottom = false;
                    break;
                }
            }
        }

        let mode = if fell_off_bottom || fell_off_top {
            // Distribution is too skewed — keep firstMode as center
            first_mode as f64
        } else {
            mode_sum / mode_duplicates as f64
        };

        // Calculate standard deviation (Bessel's correction: divide by n-1)
        let mut variance_sum = 0.0_f64;
        for i in 0..=100 {
            let diff = i as f64 - mode;
            variance_sum += diff * diff * self.gc_distribution[i];
        }
        let stdev = (variance_sum / (total - 1.0).max(1.0)).sqrt();

        // Generate theoretical normal distribution (matching FastQC's NormalDistribution)
        if stdev > 0.0 {
            let two_var = 2.0 * stdev * stdev;
            let norm_factor = 1.0 / (stdev * (2.0 * std::f64::consts::PI).sqrt());

            for i in 0..=100 {
                let diff = i as f64 - mode;
                self.theoretical[i] = norm_factor * (-diff * diff / two_var).exp() * total;
            }
        }

        // Calculate deviation (matching FastQC exactly)
        let mut deviation_sum = 0.0_f64;
        for i in 0..=100 {
            deviation_sum += (self.theoretical[i] - self.gc_distribution[i]).abs();
        }
        self.deviation_percent = deviation_sum / total * 100.0;

        let warn = config.get_limit("gc_sequence").map(|l| l.warn).unwrap_or(15.0);
        let error = config.get_limit("gc_sequence").map(|l| l.error).unwrap_or(30.0);

        self.qc_result = if self.deviation_percent > error {
            QCResult::Fail
        } else if self.deviation_percent > warn {
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
            ">>Per sequence GC content\t{}\n\
             #GC Content\tCount\n",
            self.result().label()
        );

        for i in 0..=100 {
            out.push_str(&format!("{}\t{:.1}\n", i, self.gc_distribution[i]));
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
        let mb = 50.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;

        let max_count = self
            .gc_distribution
            .iter()
            .chain(self.theoretical.iter())
            .copied()
            .fold(0.0_f64, f64::max)
            .max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Observed distribution
        svg.push_str(r##"<polyline points=""##);
        for i in 0..=100 {
            let x = ml + i as f64 / 100.0 * pw;
            let y = mt + ph * (1.0 - self.gc_distribution[i] / max_count);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#ff0000" stroke-width="2" />"##);

        // Theoretical distribution
        svg.push_str(r##"<polyline points=""##);
        for i in 0..=100 {
            let x = ml + i as f64 / 100.0 * pw;
            let y = mt + ph * (1.0 - self.theoretical[i] / max_count);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#0000ff" stroke-width="1.5" stroke-dasharray="5,3" />"##);

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
        for v in (0..=100).step_by(10) {
            let x = ml + v as f64 / 100.0 * pw;
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="middle" font-size="10">{}</text>"##,
                x, mt + ph + 15.0, v
            ));
        }

        // Legend
        svg.push_str(&format!(
            r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#ff0000" stroke-width="2" />"##,
            ml + pw - 150.0, mt + 10.0, ml + pw - 135.0, mt + 10.0
        ));
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" font-size="10">GC count per read</text>"##,
            ml + pw - 130.0, mt + 14.0
        ));
        svg.push_str(&format!(
            r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#0000ff" stroke-width="1.5" stroke-dasharray="5,3" />"##,
            ml + pw - 150.0, mt + 25.0, ml + pw - 135.0, mt + 25.0
        ));
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" font-size="10">Theoretical distribution</text>"##,
            ml + pw - 130.0, mt + 29.0
        ));

        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">GC distribution over all sequences</text>"##,
            width / 2.0
        ));

        svg.push_str("</svg>");
        svg
    }
}
