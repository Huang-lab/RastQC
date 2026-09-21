use super::{format_count_label, PhredEncoding, QCModule, QCResult};
use crate::config::FastQCConfig;
use crate::io::Sequence;
use std::any::Any;

pub struct PerSequenceQuality {
    /// Counts indexed by the read's mean quality *character*. That mean is the
    /// average of `u8` values, so it can never leave `0..=255` — a flat array
    /// indexes it exactly, and replaces a hash probe per read with a store.
    score_counts: [u64; 256],
    lowest_char: u8,
    // Results
    scores: Vec<u32>,
    counts: Vec<f64>,
    most_frequent_score: u32,
    qc_result: QCResult,
}

impl PerSequenceQuality {
    pub fn new() -> Self {
        PerSequenceQuality {
            score_counts: [0; 256],
            lowest_char: 255,
            scores: Vec::new(),
            counts: Vec::new(),
            most_frequent_score: 0,
            qc_result: QCResult::NotRun,
        }
    }
}

impl QCModule for PerSequenceQuality {
    fn name(&self) -> &str {
        "Per sequence quality scores"
    }

    fn key(&self) -> &str {
        "quality_sequence"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if seq.quality.is_empty() {
            return;
        }

        // One pass for both the running minimum and the sum; these were two
        // separate walks of every quality string.
        let mut lowest = self.lowest_char;
        let mut sum = 0u64;
        for &q in &seq.quality {
            if q < lowest {
                lowest = q;
            }
            sum += q as u64;
        }
        self.lowest_char = lowest;

        let avg = (sum as f64 / seq.quality.len() as f64).round() as usize;
        self.score_counts[avg.min(255)] += 1;
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        let Some(min_score) = self.score_counts.iter().position(|&c| c > 0) else {
            return;
        };
        let max_score = self.score_counts.iter().rposition(|&c| c > 0).unwrap();

        let offset = PhredEncoding::detect(self.lowest_char).offset() as u32;

        let mut max_count = 0u64;
        let mut most_frequent = 0u32;

        for score in min_score..=max_score {
            let adjusted = (score as u32).saturating_sub(offset);
            let count = self.score_counts[score];
            self.scores.push(adjusted);
            self.counts.push(count as f64);

            if count > max_count {
                max_count = count;
                most_frequent = adjusted;
            }
        }

        self.most_frequent_score = most_frequent;

        let warn = config
            .get_limit("quality_sequence")
            .map(|l| l.warn)
            .unwrap_or(27.0) as u32;
        let error = config
            .get_limit("quality_sequence")
            .map(|l| l.error)
            .unwrap_or(20.0) as u32;

        self.qc_result = if most_frequent <= error {
            QCResult::Fail
        } else if most_frequent <= warn {
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
            ">>Per sequence quality scores\t{}\n\
             #Quality\tCount\n",
            self.result().label()
        );

        for (i, &score) in self.scores.iter().enumerate() {
            out.push_str(&format!("{}\t{:.1}\n", score, self.counts[i]));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        if self.scores.is_empty() {
            return String::new();
        }

        let width = 800.0_f64;
        let height = 400.0_f64;
        let ml = 60.0;
        let mr = 20.0;
        let mt = 30.0;
        let mb = 50.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;

        let max_count = self.counts.iter().copied().fold(0.0_f64, f64::max).max(1.0);
        let min_score = *self.scores.first().unwrap() as f64;
        let max_score = *self.scores.last().unwrap() as f64;
        let score_range = (max_score - min_score).max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Line
        svg.push_str(r##"<polyline points=""##);
        for (i, &score) in self.scores.iter().enumerate() {
            let x = ml + (score as f64 - min_score) / score_range * pw;
            let y = mt + ph * (1.0 - self.counts[i] / max_count);
            if i > 0 {
                svg.push(' ');
            }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#ff0000" stroke-width="2" />"##);

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
            let val = max_count * frac;
            let y = mt + ph * (1.0 - frac);
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" font-size="10">{}</text>"##,
                ml - 5.0, y, format_count_label(val)
            ));
        }

        // X-axis ticks
        {
            let step = ((max_score - min_score) / 10.0).max(1.0).ceil() as u32;
            let mut v = (min_score as u32 / step) * step;
            while v <= max_score as u32 {
                if v >= min_score as u32 {
                    let x = ml + (v as f64 - min_score) / score_range * pw;
                    svg.push_str(&format!(
                        r##"<text x="{}" y="{}" text-anchor="middle" font-size="10">{}</text>"##,
                        x,
                        mt + ph + 15.0,
                        v
                    ));
                }
                v += step;
            }
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Quality score distribution over all sequences</text>"##,
            width / 2.0
        ));

        // Y axis label
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">Count</text>"##,
            mt + ph / 2.0,
            mt + ph / 2.0
        ));

        // X axis label
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" text-anchor="middle" font-size="11">Mean Sequence Quality (Phred Score)</text>"##,
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
            for (mine, &theirs) in self.score_counts.iter_mut().zip(other.score_counts.iter()) {
                *mine += theirs;
            }
            self.lowest_char = self.lowest_char.min(other.lowest_char);
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
