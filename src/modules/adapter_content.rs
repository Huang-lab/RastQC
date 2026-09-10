use super::{format_pct_label, BaseGroup, QCModule, QCResult};
use crate::config::FastQCConfig;
use crate::io::Sequence;
use crate::report::html_escape;
use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
use std::any::Any;

struct AdapterTracker {
    name: String,
    sequence: Vec<u8>,
    positions: Vec<u64>,
}

pub struct AdapterContent {
    adapters: Vec<AdapterTracker>,
    /// All adapters in one Aho-Corasick automaton, used purely to answer
    /// "does this read contain any adapter at all?" in a single SIMD-
    /// accelerated pass. Almost every read in a typical library contains
    /// none, and that pass replaces what was a brute-force scan of every
    /// start position for every adapter — which profiling showed to be over
    /// half of all CPU time spent on QC.
    ///
    /// `None` when there are no (non-empty) adapters to look for.
    matcher: Option<AhoCorasick>,
    /// Start of each adapter's first match in the current read, `usize::MAX`
    /// for "not present". Kept as a field so the rare read that does contain
    /// an adapter doesn't allocate.
    first_match: Vec<usize>,
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
        let adapters: Vec<AdapterTracker> = config
            .adapters
            .iter()
            .map(|a| AdapterTracker {
                name: a.name.clone(),
                sequence: a.sequence.as_bytes().to_vec(),
                positions: Vec::new(),
            })
            .collect();

        // An empty adapter sequence would match at every position, which is
        // meaningless output; such an entry is a malformed adapter-list line,
        // so it simply never matches. Pattern indices still line up with
        // `adapters` because empty entries are replaced, not dropped.
        let patterns: Vec<&[u8]> = adapters
            .iter()
            .map(|a| {
                if a.sequence.is_empty() {
                    &b"\0"[..]
                } else {
                    a.sequence.as_slice()
                }
            })
            .collect();
        let matcher = if patterns.is_empty() {
            None
        } else {
            // Case-insensitive so a lowercase read (or adapter list) still
            // matches, as the byte-wise comparison it replaces did.
            match AhoCorasickBuilder::new()
                .ascii_case_insensitive(true)
                .build(&patterns)
            {
                Ok(matcher) => Some(matcher),
                Err(e) => {
                    eprintln!(
                        "Warning: could not build the adapter matcher ({e}); \
                         falling back to a direct scan, which is slower."
                    );
                    None
                }
            }
        };

        AdapterContent {
            first_match: Vec::with_capacity(adapters.len()),
            adapters,
            matcher,
            total_count: 0,
            max_length: 0,
            groups: Vec::new(),
            enrichments: Vec::new(),
            max_enrichment: 0.0,
            qc_result: QCResult::NotRun,
        }
    }

    const MAX_TRACKED_POSITIONS: usize = 1000;

    /// Direct search for each adapter's first occurrence, used only when no
    /// Aho-Corasick automaton could be built. This is what the module did for
    /// every read before the automaton existed.
    fn count_adapters_by_scanning(adapters: &mut [AdapterTracker], window: &[u8], seq_len: usize) {
        for adapter in adapters.iter_mut() {
            let adapter_len = adapter.sequence.len();
            if adapter_len == 0 || seq_len < adapter_len {
                continue;
            }
            let found = (0..=(seq_len - adapter_len)).find(|&start| {
                window[start..start + adapter_len]
                    .iter()
                    .zip(adapter.sequence.iter())
                    .all(|(&b, &a)| b.eq_ignore_ascii_case(&a))
            });
            if let Some(pos) = found {
                for count in &mut adapter.positions[pos..seq_len] {
                    *count += 1;
                }
            }
        }
    }

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

        // Split the borrows so the adapter trackers and the scratch buffer can
        // be touched alongside the shared matcher.
        let Self {
            adapters,
            matcher,
            first_match,
            ..
        } = self;

        for adapter in adapters.iter_mut() {
            Self::ensure_length(&mut adapter.positions, seq_len);
        }

        let window = &seq.sequence[..seq_len];
        let Some(matcher) = matcher.as_ref() else {
            // No automaton (an oversized custom adapter list can exceed
            // Aho-Corasick's build limits). Reporting 0% adapter at every
            // position would be indistinguishable from a genuinely clean
            // library, so fall back to scanning directly rather than
            // silently reporting nothing.
            Self::count_adapters_by_scanning(adapters, window, seq_len);
            return;
        };
        // One pass settles the common case: no adapter anywhere in this read.
        if matcher.find(window).is_none() {
            return;
        }

        // This read does contain at least one adapter, so locate where each
        // one first occurs. Overlapping iteration matters here: FastQC
        // searches for every adapter independently, so an adapter starting
        // inside another adapter's match must still be found.
        first_match.clear();
        first_match.resize(adapters.len(), usize::MAX);
        for m in matcher.find_overlapping_iter(window) {
            let slot = &mut first_match[m.pattern().as_usize()];
            if *slot == usize::MAX {
                *slot = m.start();
            }
        }

        // Once an adapter appears, every base from that point on is inside it.
        for (adapter, &pos) in adapters.iter_mut().zip(first_match.iter()) {
            if pos == usize::MAX {
                continue;
            }
            for count in &mut adapter.positions[pos..seq_len] {
                *count += 1;
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
        let mut out = format!(">>Adapter Content\t{}\n#Position", self.result().label());

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

        let colors = [
            "#ff0000", "#0000ff", "#00cc00", "#ff8800", "#8800ff", "#00cccc",
        ];

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
                if gi > 0 {
                    svg.push(' ');
                }
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
            mt + ph,
            ml + pw,
            mt + ph
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
                ml + pw - 200.0,
                ml + pw - 185.0
            ));
            svg.push_str(&format!(
                r##"<text x="{}" y="{y}" dominant-baseline="middle" font-size="9">{}</text>"##,
                ml + pw - 180.0,
                html_escape(&adapter.name)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::Sequence;

    fn config_with_adapter_name(name: &str) -> FastQCConfig {
        let dir = std::env::temp_dir().join(format!(
            "rastqc-adapter-test-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("adapters.tsv");
        std::fs::write(&path, format!("{}\tAGATCGGAAGAG\n", name)).unwrap();
        FastQCConfig::new(None, Some(&path), None, 7, false, 50).unwrap()
    }

    fn seq(bases: &[u8]) -> Sequence {
        Sequence {
            header: "r".to_string(),
            sequence: bases.to_vec(),
            quality: vec![b'I'; bases.len()],
            filtered: false,
        }
    }

    #[test]
    fn svg_chart_escapes_adapter_names_from_config() {
        // Regression test: an adapter name containing HTML/SVG metacharacters
        // (loaded verbatim from a user-supplied --adapters file, which does
        // no sanitization) used to be written raw into the SVG legend,
        // letting a crafted adapter name inject markup/script into the
        // generated HTML report.
        let config = config_with_adapter_name("Evil</text><script>alert(1)</script>");
        let mut module = AdapterContent::new(&config);
        module.process_sequence(&seq(b"AGATCGGAAGAGACGT"));
        module.calculate_results(&config);

        let svg = module.svg_chart();
        assert!(!svg.contains("<script>"));
        assert!(svg.contains("&lt;script&gt;"));
    }
}
