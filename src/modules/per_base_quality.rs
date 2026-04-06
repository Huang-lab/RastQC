use std::any::Any;
use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{BaseGroup, PhredEncoding, QCModule, QCResult};

/// Tracks quality value frequencies at each position
struct QualityCount {
    counts: [u64; 128],
    total: u64,
}

impl QualityCount {
    fn new() -> Self {
        QualityCount {
            counts: [0; 128],
            total: 0,
        }
    }

    fn add(&mut self, quality_char: u8) {
        let idx = (quality_char as usize).min(127);
        self.counts[idx] += 1;
        self.total += 1;
    }

    fn percentile(&self, fraction: f64, offset: u8) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let target = (self.total as f64 * fraction).ceil() as u64;
        let mut cumulative = 0u64;
        for i in 0..128 {
            cumulative += self.counts[i];
            if cumulative >= target {
                return i as f64 - offset as f64;
            }
        }
        0.0
    }

    fn mean(&self, offset: u8) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let sum: f64 = self
            .counts
            .iter()
            .enumerate()
            .map(|(i, &c)| (i as f64 - offset as f64) * c as f64)
            .sum();
        sum / self.total as f64
    }

    fn min_char(&self) -> u8 {
        for i in 0..128 {
            if self.counts[i] > 0 {
                return i as u8;
            }
        }
        255
    }
}

pub struct PerBaseQuality {
    quality_counts: Vec<QualityCount>,
    // Computed results
    groups: Vec<BaseGroup>,
    means: Vec<f64>,
    medians: Vec<f64>,
    lower_quartiles: Vec<f64>,
    upper_quartiles: Vec<f64>,
    p10: Vec<f64>,
    p90: Vec<f64>,
    qc_result: QCResult,
}

impl PerBaseQuality {
    pub fn new() -> Self {
        PerBaseQuality {
            quality_counts: Vec::new(),
            groups: Vec::new(),
            means: Vec::new(),
            medians: Vec::new(),
            lower_quartiles: Vec::new(),
            upper_quartiles: Vec::new(),
            p10: Vec::new(),
            p90: Vec::new(),
            qc_result: QCResult::NotRun,
        }
    }

    /// Maximum number of positions to track individually.
    /// Beyond this, per-position data would be grouped anyway for display.
    /// Capping here saves ~100MB/worker for long reads with 100kb+ sequences.
    const MAX_TRACKED_POSITIONS: usize = 1000;

    fn ensure_length(&mut self, len: usize) {
        let target = len.min(Self::MAX_TRACKED_POSITIONS);
        while self.quality_counts.len() < target {
            self.quality_counts.push(QualityCount::new());
        }
    }
}

impl QCModule for PerBaseQuality {
    fn name(&self) -> &str {
        "Per base sequence quality"
    }

    fn key(&self) -> &str {
        "quality_base"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        let len = seq.quality.len().min(Self::MAX_TRACKED_POSITIONS);
        self.ensure_length(len);
        for (i, &q) in seq.quality[..len].iter().enumerate() {
            self.quality_counts[i].add(q);
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.quality_counts.is_empty() {
            return;
        }

        // Detect encoding
        let min_char = self
            .quality_counts
            .iter()
            .map(|qc| qc.min_char())
            .min()
            .unwrap_or(33);
        let encoding = PhredEncoding::detect(min_char);
        let offset = encoding.offset();

        self.groups = BaseGroup::make_groups(self.quality_counts.len());

        for group in &self.groups {
            // Merge quality counts across the group
            let mut merged = QualityCount::new();
            for pos in group.start..=group.end {
                if pos < self.quality_counts.len() {
                    for i in 0..128 {
                        merged.counts[i] += self.quality_counts[pos].counts[i];
                        merged.total += self.quality_counts[pos].counts[i];
                    }
                }
            }

            self.means.push(merged.mean(offset));
            self.medians.push(merged.percentile(0.5, offset));
            self.lower_quartiles.push(merged.percentile(0.25, offset));
            self.upper_quartiles.push(merged.percentile(0.75, offset));
            self.p10.push(merged.percentile(0.1, offset));
            self.p90.push(merged.percentile(0.9, offset));
        }

        // Determine pass/warn/fail
        let warn_lower = config
            .get_limit("quality_base_lower")
            .map(|l| l.warn)
            .unwrap_or(10.0);
        let error_lower = config
            .get_limit("quality_base_lower")
            .map(|l| l.error)
            .unwrap_or(5.0);
        let warn_median = config
            .get_limit("quality_base_median")
            .map(|l| l.warn)
            .unwrap_or(25.0);
        let error_median = config
            .get_limit("quality_base_median")
            .map(|l| l.error)
            .unwrap_or(20.0);

        self.qc_result = QCResult::Pass;

        for i in 0..self.groups.len() {
            if self.lower_quartiles[i] < error_lower || self.medians[i] < error_median {
                self.qc_result = QCResult::Fail;
                break;
            }
            if self.lower_quartiles[i] < warn_lower || self.medians[i] < warn_median {
                self.qc_result = QCResult::Warn;
            }
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Per base sequence quality\t{}\n\
             #Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n",
            self.result().label()
        );

        for (i, group) in self.groups.iter().enumerate() {
            out.push_str(&format!(
                "{}\t{:.2}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\n",
                group.label(),
                self.means[i],
                self.medians[i],
                self.lower_quartiles[i],
                self.upper_quartiles[i],
                self.p10[i],
                self.p90[i]
            ));
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        if self.groups.is_empty() {
            return String::new();
        }

        let width = 800.0_f64;
        let height = 400.0_f64;
        let margin_left = 60.0;
        let margin_right = 20.0;
        let margin_top = 30.0;
        let margin_bottom = 80.0;
        let plot_w = width - margin_left - margin_right;
        let plot_h = height - margin_top - margin_bottom;

        let n = self.groups.len();
        let max_q: f64 = self
            .p90
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
            .max(40.0);

        let bar_w = plot_w / n as f64;

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Background quality zones
        let good_y = margin_top;
        let good_h = plot_h * (1.0 - 28.0 / max_q);
        let ok_y = margin_top + good_h;
        let ok_h = plot_h * ((28.0 - 20.0) / max_q);
        let bad_y = ok_y + ok_h;
        let bad_h = plot_h - good_h - ok_h;

        svg.push_str(&format!(
            r##"<rect x="{}" y="{}" width="{}" height="{}" fill="#e6ffe6" />"##,
            margin_left, good_y, plot_w, good_h
        ));
        svg.push_str(&format!(
            r##"<rect x="{}" y="{}" width="{}" height="{}" fill="#ffffcc" />"##,
            margin_left, ok_y, plot_w, ok_h
        ));
        svg.push_str(&format!(
            r##"<rect x="{}" y="{}" width="{}" height="{}" fill="#ffe6e6" />"##,
            margin_left, bad_y, plot_w, bad_h
        ));

        // Box plots
        for i in 0..n {
            let x = margin_left + i as f64 * bar_w;
            let cx = x + bar_w / 2.0;

            let y_fn = |val: f64| -> f64 { margin_top + plot_h * (1.0 - val / max_q) };

            let y_p90 = y_fn(self.p90[i]);
            let y_uq = y_fn(self.upper_quartiles[i]);
            let y_med = y_fn(self.medians[i]);
            let y_lq = y_fn(self.lower_quartiles[i]);
            let y_p10 = y_fn(self.p10[i]);

            // Whisker lines
            svg.push_str(&format!(
                r##"<line x1="{cx}" y1="{}" x2="{cx}" y2="{}" stroke="#333" stroke-width="0.5" />"##,
                y_p90, y_uq
            ));
            svg.push_str(&format!(
                r##"<line x1="{cx}" y1="{}" x2="{cx}" y2="{}" stroke="#333" stroke-width="0.5" />"##,
                y_lq, y_p10
            ));

            // Whisker caps
            let cap_w = bar_w * 0.3;
            svg.push_str(&format!(
                r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#333" stroke-width="0.5" />"##,
                cx - cap_w, y_p90, cx + cap_w, y_p90
            ));
            svg.push_str(&format!(
                r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#333" stroke-width="0.5" />"##,
                cx - cap_w, y_p10, cx + cap_w, y_p10
            ));

            // IQR box
            let box_w = bar_w * 0.6;
            let box_h = (y_lq - y_uq).abs().max(1.0);
            svg.push_str(&format!(
                r##"<rect x="{}" y="{}" width="{}" height="{}" fill="#ffd700" stroke="#333" stroke-width="0.5" />"##,
                cx - box_w / 2.0,
                y_uq.min(y_lq),
                box_w,
                box_h
            ));

            // Median line
            svg.push_str(&format!(
                r##"<line x1="{}" y1="{y_med}" x2="{}" y2="{y_med}" stroke="red" stroke-width="1.5" />"##,
                cx - box_w / 2.0,
                cx + box_w / 2.0
            ));
        }

        // Mean line
        svg.push_str("<polyline points=\"");
        for i in 0..n {
            let x = margin_left + i as f64 * bar_w + bar_w / 2.0;
            let y = margin_top + plot_h * (1.0 - self.means[i] / max_q);
            if i > 0 {
                svg.push(' ');
            }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str("\" fill=\"none\" stroke=\"#2222ff\" stroke-width=\"1\" />");

        // Axes
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mb}" stroke="black" stroke-width="1" />"##,
            ml = margin_left,
            mt = margin_top,
            mb = margin_top + plot_h
        ));
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mb}" x2="{mr}" y2="{mb}" stroke="black" stroke-width="1" />"##,
            ml = margin_left,
            mr = margin_left + plot_w,
            mb = margin_top + plot_h
        ));

        // Y-axis labels
        let step = if max_q > 30.0 { 10.0 } else { 5.0 };
        let mut v = 0.0;
        while v <= max_q {
            let y = margin_top + plot_h * (1.0 - v / max_q);
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" font-size="10">{}</text>"##,
                margin_left - 5.0, y, v as u32
            ));
            v += step;
        }

        // X-axis labels (show subset)
        let label_step = (n / 20).max(1);
        for i in (0..n).step_by(label_step) {
            let x = margin_left + i as f64 * bar_w + bar_w / 2.0;
            let y = margin_top + plot_h + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{}</text>"##,
                self.groups[i].label()
            ));
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Quality scores across all bases</text>"##,
            width / 2.0
        ));

        // Axis labels
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">Quality Score</text>"##,
            margin_top + plot_h / 2.0,
            margin_top + plot_h / 2.0
        ));
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" text-anchor="middle" font-size="11">Position in read (bp)</text>"##,
            margin_left + plot_w / 2.0,
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
            self.ensure_length(other.quality_counts.len());
            for (i, oqc) in other.quality_counts.iter().enumerate() {
                for j in 0..128 {
                    self.quality_counts[i].counts[j] += oqc.counts[j];
                }
                self.quality_counts[i].total += oqc.total;
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
