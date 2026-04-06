use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{BaseGroup, QCModule, QCResult, format_pct_label};
use std::any::Any;

pub struct NContent {
    n_counts: Vec<u64>,
    total_counts: Vec<u64>,
    // Results
    groups: Vec<BaseGroup>,
    percentages: Vec<f64>,
    qc_result: QCResult,
}

impl NContent {
    pub fn new() -> Self {
        NContent {
            n_counts: Vec::new(),
            total_counts: Vec::new(),
            groups: Vec::new(),
            percentages: Vec::new(),
            qc_result: QCResult::NotRun,
        }
    }

    const MAX_TRACKED_POSITIONS: usize = 1000;

    fn ensure_length(&mut self, len: usize) {
        let target = len.min(Self::MAX_TRACKED_POSITIONS);
        while self.n_counts.len() < target {
            self.n_counts.push(0);
            self.total_counts.push(0);
        }
    }
}

impl QCModule for NContent {
    fn name(&self) -> &str {
        "Per base N content"
    }

    fn key(&self) -> &str {
        "n_content"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        let len = seq.len().min(Self::MAX_TRACKED_POSITIONS);
        self.ensure_length(len);
        for (i, &b) in seq.sequence[..len].iter().enumerate() {
            self.total_counts[i] += 1;
            if b == b'N' || b == b'n' {
                self.n_counts[i] += 1;
            }
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.n_counts.is_empty() {
            return;
        }

        self.groups = BaseGroup::make_groups(self.n_counts.len());

        let warn = config.get_limit("n_content").map(|l| l.warn).unwrap_or(5.0);
        let error = config.get_limit("n_content").map(|l| l.error).unwrap_or(20.0);

        self.qc_result = QCResult::Pass;

        for group in &self.groups {
            let mut n_sum = 0u64;
            let mut total_sum = 0u64;
            for pos in group.start..=group.end {
                if pos < self.n_counts.len() {
                    n_sum += self.n_counts[pos];
                    total_sum += self.total_counts[pos];
                }
            }
            let pct = if total_sum > 0 {
                n_sum as f64 / total_sum as f64 * 100.0
            } else {
                0.0
            };
            self.percentages.push(pct);

            if pct > error {
                self.qc_result = QCResult::Fail;
            } else if pct > warn && self.qc_result != QCResult::Fail {
                self.qc_result = QCResult::Warn;
            }
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Per base N content\t{}\n\
             #Base\tN-Count\n",
            self.result().label()
        );

        for (i, group) in self.groups.iter().enumerate() {
            out.push_str(&format!(
                "{}\t{:.6}\n",
                group.label(),
                self.percentages.get(i).unwrap_or(&0.0)
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
        let ml = 60.0;
        let mr = 20.0;
        let mt = 30.0;
        let mb = 80.0;
        let pw = width - ml - mr;
        let ph = height - mt - mb;
        let n = self.groups.len();

        let max_pct = self.percentages.iter().copied().fold(0.0_f64, f64::max).max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        svg.push_str(r##"<polyline points=""##);
        for (i, &pct) in self.percentages.iter().enumerate() {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph * (1.0 - pct / max_pct);
            if i > 0 { svg.push(' '); }
            svg.push_str(&format!("{:.1},{:.1}", x, y));
        }
        svg.push_str(r##"" fill="none" stroke="#ff0000" stroke-width="1.5" />"##);

        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{}" stroke="black" />"##,
            mt + ph
        ));
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{}" x2="{}" y2="{}" stroke="black" />"##,
            mt + ph, ml + pw, mt + ph
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

        // X-axis tick labels
        let step = (n / 15).max(1);
        for i in (0..n).step_by(step) {
            let x = ml + (i as f64 + 0.5) / n as f64 * pw;
            let y = mt + ph + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{}</text>"##,
                self.groups[i].label()
            ));
        }

        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">N content across all bases</text>"##,
            width / 2.0
        ));

        // Y axis label
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">% N</text>"##,
            mt + ph / 2.0,
            mt + ph / 2.0
        ));

        // X axis label
        svg.push_str(&format!(
            r##"<text x="{}" y="{}" text-anchor="middle" font-size="11">Position in read (bp)</text>"##,
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
            self.ensure_length(other.n_counts.len());
            for i in 0..other.n_counts.len() {
                self.n_counts[i] += other.n_counts[i];
                self.total_counts[i] += other.total_counts[i];
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
