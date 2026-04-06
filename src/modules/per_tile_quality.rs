use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{BaseGroup, QCModule, QCResult};
use std::any::Any;
use std::collections::HashMap;

pub struct PerTileQuality {
    /// tile_id -> vec of (quality_sum, count) per position
    tile_data: HashMap<u32, Vec<(f64, u64)>>,
    total_sequences: u64,
    max_length: usize,
    split_position: Option<usize>,
    max_tiles: usize,
    gave_up: bool,
    // Results
    groups: Vec<BaseGroup>,
    tile_ids: Vec<u32>,
    /// Normalized means: deviation from average per group
    normalized: Vec<Vec<f64>>,
    max_deviation: f64,
    qc_result: QCResult,
}

impl PerTileQuality {
    pub fn new() -> Self {
        PerTileQuality {
            tile_data: HashMap::new(),
            total_sequences: 0,
            max_length: 0,
            split_position: None,
            max_tiles: 2500,
            gave_up: false,
            groups: Vec::new(),
            tile_ids: Vec::new(),
            normalized: Vec::new(),
            max_deviation: 0.0,
            qc_result: QCResult::NotRun,
        }
    }

    fn extract_tile(&mut self, header: &str) -> Option<u32> {
        let parts: Vec<&str> = header.split(':').collect();

        if let Some(&pos) = self.split_position.as_ref() {
            return parts.get(pos).and_then(|s| s.parse().ok());
        }

        // Auto-detect: Illumina 1.8+ has 7+ fields, tile at position 4
        let tile_pos = if parts.len() >= 7 {
            4
        } else if parts.len() >= 5 {
            2
        } else {
            return None;
        };

        self.split_position = Some(tile_pos);
        parts.get(tile_pos).and_then(|s| s.parse().ok())
    }
}

impl QCModule for PerTileQuality {
    fn name(&self) -> &str {
        "Per tile sequence quality"
    }

    fn key(&self) -> &str {
        "tile"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        if self.gave_up {
            return;
        }

        self.total_sequences += 1;

        // Sampling: process all first 10000, then 10%
        if self.total_sequences > 10000 && self.total_sequences % 10 != 0 {
            return;
        }

        let tile = match self.extract_tile(&seq.header) {
            Some(t) => t,
            None => return,
        };

        if self.tile_data.len() >= self.max_tiles && !self.tile_data.contains_key(&tile) {
            self.gave_up = true;
            return;
        }

        let len = seq.quality.len().min(1000);
        if len > self.max_length {
            self.max_length = len;
        }

        let entry = self
            .tile_data
            .entry(tile)
            .or_insert_with(|| Vec::new());

        while entry.len() < len {
            entry.push((0.0, 0));
        }

        for (i, &q) in seq.quality[..len].iter().enumerate() {
            entry[i].0 += q as f64;
            entry[i].1 += 1;
        }
    }

    fn calculate_results(&mut self, config: &FastQCConfig) {
        if self.tile_data.is_empty() {
            self.qc_result = QCResult::Pass;
            return;
        }

        // Detect encoding
        let offset = 33u8; // default Sanger

        self.groups = BaseGroup::make_groups(self.max_length);
        let ng = self.groups.len();

        // Compute mean per group across all tiles
        let mut global_mean = vec![0.0_f64; ng];
        let mut global_count = vec![0u64; ng];

        self.tile_ids = self.tile_data.keys().copied().collect();
        self.tile_ids.sort();

        for tile_data in self.tile_data.values() {
            for (gi, group) in self.groups.iter().enumerate() {
                for pos in group.start..=group.end {
                    if pos < tile_data.len() {
                        global_mean[gi] += tile_data[pos].0 - offset as f64 * tile_data[pos].1 as f64;
                        global_count[gi] += tile_data[pos].1;
                    }
                }
            }
        }

        for gi in 0..ng {
            if global_count[gi] > 0 {
                global_mean[gi] /= global_count[gi] as f64;
            }
        }

        // Compute per-tile normalized means
        self.max_deviation = 0.0;
        for &tile_id in &self.tile_ids {
            let tile_data = &self.tile_data[&tile_id];
            let mut row = Vec::with_capacity(ng);

            for (gi, group) in self.groups.iter().enumerate() {
                let mut sum = 0.0;
                let mut count = 0u64;
                for pos in group.start..=group.end {
                    if pos < tile_data.len() {
                        sum += tile_data[pos].0 - offset as f64 * tile_data[pos].1 as f64;
                        count += tile_data[pos].1;
                    }
                }
                let tile_mean = if count > 0 {
                    sum / count as f64
                } else {
                    0.0
                };
                let dev = tile_mean - global_mean[gi];
                if dev.abs() > self.max_deviation {
                    self.max_deviation = dev.abs();
                }
                row.push(dev);
            }
            self.normalized.push(row);
        }

        let warn_thresh = config.get_limit("tile").map(|l| l.warn).unwrap_or(5.0);
        let error_thresh = config.get_limit("tile").map(|l| l.error).unwrap_or(10.0);

        self.qc_result = if self.max_deviation > error_thresh {
            QCResult::Fail
        } else if self.max_deviation > warn_thresh {
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
            ">>Per tile sequence quality\t{}\n\
             #Tile\tBase\tMean\n",
            self.result().label()
        );

        for (ti, &tile_id) in self.tile_ids.iter().enumerate() {
            for (gi, group) in self.groups.iter().enumerate() {
                if ti < self.normalized.len() && gi < self.normalized[ti].len() {
                    out.push_str(&format!(
                        "{}\t{}\t{:.2}\n",
                        tile_id,
                        group.label(),
                        self.normalized[ti][gi]
                    ));
                }
            }
        }
        out.push_str(">>END_MODULE\n");
        out
    }

    fn svg_chart(&self) -> String {
        if self.tile_ids.is_empty() || self.groups.is_empty() {
            return String::new();
        }

        let width = 800.0_f64;
        let height = 400.0_f64;
        let margin_left = 80.0;
        let margin_right = 20.0;
        let margin_top = 30.0;
        let margin_bottom = 80.0;
        let plot_w = width - margin_left - margin_right;
        let plot_h = height - margin_top - margin_bottom;

        let n_tiles = self.tile_ids.len();
        let n_groups = self.groups.len();

        let cell_w = plot_w / n_groups as f64;
        let cell_h = plot_h / n_tiles as f64;

        let max_dev = self.max_deviation.max(1.0);

        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" font-family="Arial, sans-serif" font-size="10">"##,
            width, height
        );

        // Heatmap cells
        for (ti, row) in self.normalized.iter().enumerate() {
            for (gi, &dev) in row.iter().enumerate() {
                let x = margin_left + gi as f64 * cell_w;
                let y = margin_top + ti as f64 * cell_h;

                // Color: blue (negative = bad) to white (0) to red (positive deviation not typical)
                // Actually in FastQC: red = below average (bad), blue = above average
                let norm = (dev / max_dev).clamp(-1.0, 1.0);
                let (r, g, b) = if norm < 0.0 {
                    let t = -norm;
                    ((255.0 * (1.0 - t * 0.3)) as u8, (255.0 * (1.0 - t)) as u8, (255.0 * (1.0 - t)) as u8)
                } else {
                    let t = norm;
                    ((255.0 * (1.0 - t)) as u8, (255.0 * (1.0 - t)) as u8, (255.0 * (1.0 - t * 0.3)) as u8)
                };

                svg.push_str(&format!(
                    r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="rgb({},{},{})" />"##,
                    x, y, cell_w.max(1.0), cell_h.max(1.0), r, g, b
                ));
            }
        }

        // Y-axis tick labels (tile IDs, show subset)
        let tile_step = (n_tiles / 10).max(1);
        for ti in (0..n_tiles).step_by(tile_step) {
            let y = margin_top + ti as f64 * cell_h + cell_h / 2.0;
            svg.push_str(&format!(
                r##"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" font-size="9">{}</text>"##,
                margin_left - 5.0, y, self.tile_ids[ti]
            ));
        }

        // X-axis tick labels (base positions, show subset)
        let x_step = (n_groups / 15).max(1);
        for gi in (0..n_groups).step_by(x_step) {
            let x = margin_left + gi as f64 * cell_w + cell_w / 2.0;
            let y = margin_top + plot_h + 15.0;
            svg.push_str(&format!(
                r##"<text x="{x}" y="{y}" text-anchor="end" transform="rotate(-45 {x} {y})" font-size="9">{}</text>"##,
                self.groups[gi].label()
            ));
        }

        // Title
        svg.push_str(&format!(
            r##"<text x="{}" y="18" text-anchor="middle" font-size="13" font-weight="bold">Quality per tile</text>"##,
            width / 2.0
        ));

        // Y axis label
        svg.push_str(&format!(
            r##"<text x="15" y="{}" text-anchor="middle" transform="rotate(-90 15 {})" font-size="11">Tile</text>"##,
            margin_top + plot_h / 2.0,
            margin_top + plot_h / 2.0
        ));

        // X axis label
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
            self.total_sequences += other.total_sequences;
            self.max_length = self.max_length.max(other.max_length);
            if other.gave_up { self.gave_up = true; }
            for (tile_id, other_positions) in &other.tile_data {
                let entry = self.tile_data.entry(*tile_id).or_insert_with(Vec::new);
                while entry.len() < other_positions.len() {
                    entry.push((0.0, 0));
                }
                for (i, (sum, count)) in other_positions.iter().enumerate() {
                    entry[i].0 += sum;
                    entry[i].1 += count;
                }
            }
        }
    }

    fn supports_merge(&self) -> bool {
        true
    }
}
