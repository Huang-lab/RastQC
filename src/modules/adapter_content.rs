use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{BaseGroup, QCModule, QCResult, format_pct_label};
use std::any::Any;

struct AdapterTracker {
    name: String,
    sequence: Vec<u8>,
    positions: Vec<u64>,
}

pub struct AdapterContent {
    adapters: Vec<AdapterTracker>,
    total_count: u64,
    max_length: usize,
    // Results
    groups: Vec<BaseGroup>,
    /// enrichments[adapter_idx][group_idx] as percentage
    enrichments: Vec<Vec<f64>>,
    max_enrichment: f64,
    qc_result: QCResult,
}

impl AdapterContent {
    pub fn new(config: &FastQCConfig) -> Self {
        let adapters = config
            .adapters
            .iter()
            .map(|a| AdapterTracker {
                name: a.name.clone(),
                sequence: a.sequence.as_bytes().to_vec(),
                positions: Vec::new(),
            })
            .collect();

        AdapterContent {
            adapters,
            total_count: 0,
            max_length: 0,
            groups: Vec::new(),
            enrichments: Vec::new(),
            max_enrichment: 0.0,
            qc_result: QCResult::NotRun,
        }
    }

    const MAX_TRACKED_POSITIONS: usize = 1000;

    fn ensure_length(positions: &mut Vec<u64>, len: usize) {
        let target = len.min(Self::MAX_TRACKED_POSITIONS);
        while positions.len() < target {
            positions.push(0);
        }
    }
}

impl QCModule for AdapterContent {
    fn name(&self) -> &str {
        "Adapter Content"
    }

    fn key(&self) -> &str {
        "adapter"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_count += 1;
        let seq_len = seq.sequence.len().min(Self::MAX_TRACKED_POSITIONS);
        if seq_len > self.max_length {
            self.max_length = seq_len;
        }

        for adapter in &mut self.adapters {
            Self::ensure_length(&mut adapter.positions, seq_len);

            let adapter_len = adapter.sequence.len();
            if seq_len < adapter_len {
                continue;
            }

            // Case-insensitive search without allocating an uppercase copy
            let mut found_pos = None;
            for start in 0..=(seq_len - adapter_len) {
                let mut matches = true;
                for j in 0..adapter_len {
                    if seq.sequence[start + j].to_ascii_uppercase() != adapter.sequence[j] {
                        matches = false;
                        break;
                    }
                }
                if matches {
                    found_pos = Some(start);
                    break;
                }
            }

            // If found, increment all positions from match point onwards
            if let Some(pos) = found_pos {
                for i in pos..seq_len {
                    adapter.positions[i] += 1;
                }
            }
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.total_count == 0 || self.max_length == 0 {
            return;
        }

        self.groups = BaseGroup::make_groups(self.max_length);
        let ng = self.groups.len();

        let warn = config.get_limit("adapter").map(|l| l.warn).unwrap_or(5.0);
        let error = config.get_limit("adapter").map(|l| l.error).unwrap_or(10.0);

        self.qc_result = QCResult::Pass;
        self.max_enrichment = 0.0;

        for adapter in &self.adapters {
            let mut row = Vec::with_capacity(ng);

            for group in &self.groups {
                let mut sum = 0u64;
                let mut count = 0u64;
                for pos in group.start..=group.end {
                    if pos < adapter.positions.len() {
                        sum += adapter.positions[pos];
                        count += 1;
                    }
                }
                let enrichment = if count > 0 && self.total_count > 0 {
                    (sum as f64 / count as f64) / self.total_count as f64 * 100.0
                } else {
                    0.0
                };
                row.push(enrichment);

                if enrichment > self.max_enrichment {
                    self.max_enrichment = enrichment;
                }
            }
            self.enrichments.push(row);
        }

        if self.max_enrichment > error {
            self.qc_result = QCResult::Fail;
        } else if self.max_enrichment > warn {
            self.qc_result = QCResult::Warn;
        }
    }

    fn result(&self) -> QCResult {
        self.qc_result
    }

    fn text_data(&self) -> String {
        let mut out = format!(
            ">>Adapter Content\t{}\n#Position",
            self.result().label()
        );

        for adapter in &self.adapters {
            out.push_str(&format!("\t{}", adapter.name));
        }
        out.push('\n');

        for (gi, group) in self.groups.iter().enumerate() {
            out.push_str(&group.label());
            for ai in 0..self.adapters.len() {
                let val = self
                    .enrichments
                    .get(ai)
                    .and_then(|row| row.get(gi))
                    .unwrap_or(&0.0);
                out.push_str(&format!("\t{:.6}", val));
            }
            out.push('\n');
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        if self.groups.is_empty() || self.adapters.is_empty() {
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

        let max_val = self.max_enrichment.max(1.0);

        let colors = ["#ff0000", "#0000ff", "#00cc00", "#ff8800", "#8800ff", "#00cccc"];

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="11">"##,
            width, height
        );

        for (ai, _adapter) in self.adapters.iter().enumerate() {
            let color = colors[ai % colors.len()];
            let row = &self.enrichments[ai];

            svg.push_str(r##"<polyline points=""##);
            for (gi, &val) in row.iter().enumerate() {
                let x = ml + (gi as f64 + 0.5) / n as f64 * pw;
                let y = mt + ph * (1.0 - val / max_val);
                if gi > 0 { svg.push(' '); }
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

        // Y-axis ticks
        let y_steps = 5;
        for i in 0..=y_steps {
            let frac = i as f64 / y_steps as f64;
            let val = max_val * frac;
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

        // Legend
        for (ai, adapter) in self.adapters.iter().enumerate() {
            let color = colors[ai % colors.len()];
            let y = mt + 10.0 + ai as f64 * 14.0;
            svg.push_str(&format!(
                r##"<line x1="{}" y1="{y}" x2="{}" y2="{y}" stroke="{color}" stroke-width="2" />"##,
                ml + pw - 200.0, ml + pw - 185.0
            ));
            svg.push_str(&format!(
                r##"<text x="{}" y="{y}" dominant-baseline="middle" font-size="9">{}</text>"##,
                ml + pw - 180.0, adapter.name
            ));
        }

        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Adapter Content</text>"##,
            width / 2.0
        ));

        // Y axis label
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">% Adapter</text>"##,
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
            self.total_count += other.total_count;
            if other.max_length > self.max_length {
                self.max_length = other.max_length;
            }
            for (i, other_adapter) in other.adapters.iter().enumerate() {
                if i < self.adapters.len() {
                    let self_positions = &mut self.adapters[i].positions;
                    // Extend self if other is longer
                    while self_positions.len() < other_adapter.positions.len() {
                        self_positions.push(0);
                    }
                    for (j, &val) in other_adapter.positions.iter().enumerate() {
                        self_positions[j] += val;
                    }
                }
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
