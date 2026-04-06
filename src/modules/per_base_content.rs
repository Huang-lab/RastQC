use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{BaseGroup, QCModule, QCResult};
use std::any::Any;

pub struct PerBaseContent {
    a_counts: Vec<u64>,
    t_counts: Vec<u64>,
    g_counts: Vec<u64>,
    c_counts: Vec<u64>,
    // Results
    groups: Vec<BaseGroup>,
    a_pct: Vec<f64>,
    t_pct: Vec<f64>,
    g_pct: Vec<f64>,
    c_pct: Vec<f64>,
    qc_result: QCResult,
}

impl PerBaseContent {
    pub fn new() -> Self {
        PerBaseContent {
            a_counts: Vec::new(),
            t_counts: Vec::new(),
            g_counts: Vec::new(),
            c_counts: Vec::new(),
            groups: Vec::new(),
            a_pct: Vec::new(),
            t_pct: Vec::new(),
            g_pct: Vec::new(),
            c_pct: Vec::new(),
            qc_result: QCResult::NotRun,
        }
    }

    const MAX_TRACKED_POSITIONS: usize = 1000;

    fn ensure_length(&mut self, len: usize) {
        let target = len.min(Self::MAX_TRACKED_POSITIONS);
        while self.a_counts.len() < target {
            self.a_counts.push(0);
            self.t_counts.push(0);
            self.g_counts.push(0);
            self.c_counts.push(0);
        }
    }
}

impl QCModule for PerBaseContent {
    fn name(&self) -> &str {
        "Per base sequence content"
    }

    fn key(&self) -> &str {
        "sequence"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        let len = seq.len().min(Self::MAX_TRACKED_POSITIONS);
        self.ensure_length(len);
        for (i, &b) in seq.sequence[..len].iter().enumerate() {
            match b {
                b'A' | b'a' => self.a_counts[i] += 1,
                b'T' | b't' => self.t_counts[i] += 1,
                b'G' | b'g' => self.g_counts[i] += 1,
                b'C' | b'c' => self.c_counts[i] += 1,
                _ => {}
            }
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.a_counts.is_empty() {
            return;
        }

        self.groups = BaseGroup::make_groups(self.a_counts.len());

        let warn = config.get_limit("sequence").map(|l| l.warn).unwrap_or(10.0);
        let error = config.get_limit("sequence").map(|l| l.error).unwrap_or(20.0);

        self.qc_result = QCResult::Pass;

        for group in &self.groups {
            let mut a = 0u64;
            let mut t = 0u64;
            let mut g = 0u64;
            let mut c = 0u64;
            for pos in group.start..=group.end {
                if pos < self.a_counts.len() {
                    a += self.a_counts[pos];
                    t += self.t_counts[pos];
                    g += self.g_counts[pos];
                    c += self.c_counts[pos];
                }
            }
            let total = (a + t + g + c) as f64;
            if total == 0.0 {
                self.a_pct.push(0.0);
                self.t_pct.push(0.0);
                self.g_pct.push(0.0);
                self.c_pct.push(0.0);
                continue;
            }
            let ap = a as f64 / total * 100.0;
            let tp = t as f64 / total * 100.0;
            let gp = g as f64 / total * 100.0;
            let cp = c as f64 / total * 100.0;

            self.a_pct.push(ap);
            self.t_pct.push(tp);
            self.g_pct.push(gp);
            self.c_pct.push(cp);

            let at_diff = (ap - tp).abs();
            let gc_diff = (gp - cp).abs();

            if at_diff > error || gc_diff > error {
                self.qc_result = QCResult::Fail;
            } else if (at_diff > warn || gc_diff > warn) && self.qc_result != QCResult::Fail {
                self.qc_result = QCResult::Warn;
            }
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Per base sequence content\t{}\n\
             #Base\tG\tA\tT\tC\n",
            self.result().label()
        );

        for (i, group) in self.groups.iter().enumerate() {
            out.push_str(&format!(
                "{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\n",
                group.label(),
                self.g_pct.get(i).unwrap_or(&0.0),
                self.a_pct.get(i).unwrap_or(&0.0),
                self.t_pct.get(i).unwrap_or(&0.0),
                self.c_pct.get(i).unwrap_or(&0.0)
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

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        // Lines for each base
        let colors = [
            ("#dc322f", &self.t_pct, "% T"),
            ("#268bd2", &self.c_pct, "% C"),
            ("#859900", &self.a_pct, "% A"),
            ("#333333", &self.g_pct, "% G"),
        ];

        for &(color, data, _label) in &colors {
            svg.push_str(r##"<polyline points=""##);
            for (i, &val) in data.iter().enumerate() {
                let x = ml + (i as f64 + 0.5) / n as f64 * pw;
                let y = mt + ph * (1.0 - val / 100.0);
                if i > 0 { svg.push(' '); }
                svg.push_str(&format!("{:.1},{:.1}", x, y));
            }
            svg.push_str(&format!(
                r##"" fill="none" stroke="{}" stroke-width="1.5" />"##,
                color
            ));
        }

        // Axes
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{}" stroke="black" />"##,
            mt + ph
        ));
        svg.push_str(&format!(
            r##"<line x1="{ml}" y1="{}" x2="{}" y2="{}" stroke="black" />"##,
            mt + ph, ml + pw, mt + ph
        ));

        // Y axis tick labels
        for v in (0..=100).step_by(20) {
            let y = mt + ph * (1.0 - v as f64 / 100.0);
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" font-size="10">{}%</text>"##,
                ml - 5.0, y, v
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

        // Legend
        let legend_x = ml + pw - 120.0;
        let legend_y = mt + 10.0;
        for (i, &(color, _, label)) in colors.iter().enumerate() {
            let y = legend_y + i as f64 * 15.0;
            svg.push_str(&format!(
                r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="2" />"##,
                legend_x, y, legend_x + 15.0, y, color
            ));
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" dominant-baseline="middle" font-size="10">{}</text>"##,
                legend_x + 20.0, y, label
            ));
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Sequence content across all bases</text>"##,
            width / 2.0
        ));

        // Y axis label
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">Percentage (%)</text>"##,
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
            while self.a_counts.len() < other.a_counts.len() {
                self.a_counts.push(0);
                self.t_counts.push(0);
                self.g_counts.push(0);
                self.c_counts.push(0);
            }
            for i in 0..other.a_counts.len() {
                self.a_counts[i] += other.a_counts[i];
                self.t_counts[i] += other.t_counts[i];
                self.g_counts[i] += other.g_counts[i];
                self.c_counts[i] += other.c_counts[i];
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
