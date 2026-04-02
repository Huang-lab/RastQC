pub mod basic_stats;
pub mod per_base_quality;
pub mod per_tile_quality;
pub mod per_sequence_quality;
pub mod per_base_content;
pub mod per_sequence_gc;
pub mod n_content;
pub mod sequence_length;
pub mod duplication;
pub mod overrepresented;
pub mod adapter_content;
pub mod kmer_content;
pub mod long_read_quality;

use crate::config::FastQCConfig;
use crate::io::Sequence;
use serde_json::Value as JsonValue;
use std::any::Any;

/// QC result status
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum QCResult {
    Pass,
    Warn,
    Fail,
    NotRun,
}

impl QCResult {
    pub fn icon(&self) -> &'static str {
        match self {
            QCResult::Pass => "pass",
            QCResult::Warn => "warn",
            QCResult::Fail => "fail",
            QCResult::NotRun => "na",
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            QCResult::Pass => "PASS",
            QCResult::Warn => "WARN",
            QCResult::Fail => "FAIL",
            QCResult::NotRun => "N/A",
        }
    }
}

/// Trait for all QC modules
pub trait QCModule: Send {
    /// Module name for display
    fn name(&self) -> &str;

    /// Short key for config lookup
    #[allow(dead_code)]
    fn key(&self) -> &str;

    /// Process one sequence
    fn process_sequence(&mut self, seq: &Sequence);

    /// Calculate final results after all sequences processed
    fn calculate_results(&mut self, config: &FastQCConfig);

    /// Get pass/warn/fail status
    fn result(&self) -> QCResult;

    /// Generate tab-separated text data for this module
    fn text_data(&self) -> String;

    /// Generate SVG chart content
    fn svg_chart(&self) -> String;

    /// Whether this module has a chart
    fn has_chart(&self) -> bool {
        true
    }

    /// For downcasting in merge operations.
    fn as_any_mut(&mut self) -> &mut dyn Any;

    /// Merge another module instance's accumulated state into this one.
    /// Called during intra-file parallel processing to combine chunk results.
    /// Default: no-op (modules that don't implement merge fall back to sequential).
    fn merge_from(&mut self, _other: &mut dyn QCModule) {}

    /// Whether this module supports merging parallel chunk results.
    #[allow(dead_code)]
    fn supports_merge(&self) -> bool {
        false
    }

    /// Export module data as a JSON value for MultiQC native format.
    /// Default implementation parses text_data() into a structured JSON object.
    fn json_data(&self) -> JsonValue {
        text_data_to_json(&self.text_data(), self.name(), self.result())
    }
}

/// Merge two module sets. Each module in `target` merges from the corresponding
/// module in `source`. Both sets must have the same modules in the same order
/// (as produced by ModuleFactory::create_modules).
pub fn merge_module_sets(
    target: &mut [Box<dyn QCModule>],
    source: &mut [Box<dyn QCModule>],
) {
    assert_eq!(target.len(), source.len());
    for (t, s) in target.iter_mut().zip(source.iter_mut()) {
        t.merge_from(s.as_mut());
    }
}

/// Parse tab-separated text_data into a JSON object with headers and rows.
fn text_data_to_json(text: &str, module_name: &str, result: QCResult) -> JsonValue {
    let mut headers: Vec<String> = Vec::new();
    let mut rows: Vec<JsonValue> = Vec::new();

    for line in text.lines() {
        if line.starts_with(">>END_MODULE") {
            break;
        }
        if line.starts_with(">>") {
            continue;
        }
        if line.starts_with('#') {
            headers = line
                .trim_start_matches('#')
                .split('\t')
                .map(|s| s.trim().to_string())
                .collect();
            continue;
        }
        if headers.is_empty() {
            continue;
        }
        let values: Vec<&str> = line.split('\t').collect();
        let mut row = serde_json::Map::new();
        for (i, val) in values.iter().enumerate() {
            let key = headers.get(i).cloned().unwrap_or_else(|| format!("col_{}", i));
            // Try to parse as number
            if let Ok(n) = val.parse::<f64>() {
                row.insert(key, serde_json::json!(n));
            } else {
                row.insert(key, serde_json::json!(val));
            }
        }
        rows.push(JsonValue::Object(row));
    }

    serde_json::json!({
        "module_name": module_name,
        "status": result.label(),
        "data": rows
    })
}

/// Base group for aggregating positions
#[derive(Debug, Clone)]
pub struct BaseGroup {
    pub start: usize,
    pub end: usize,
}

impl BaseGroup {
    pub fn label(&self) -> String {
        if self.start == self.end {
            format!("{}", self.start + 1)
        } else {
            format!("{}-{}", self.start + 1, self.end + 1)
        }
    }

    /// Create base groups for a given max length, matching FastQC's binning:
    /// 1-9: individual positions
    /// 10-49: every 5 bp (10-14, 15-19, etc.)
    /// 50-74: every 5 bp (50-54, 55-59, etc.)
    /// 75-99: every 5 bp (75-79, 80-84, etc.)
    /// 100+: every 5 bp until long read threshold, then larger bins
    pub fn make_groups(max_length: usize) -> Vec<BaseGroup> {
        if max_length == 0 {
            return vec![];
        }

        let mut groups = Vec::new();

        // 1-9 individually
        let individual = 9.min(max_length);
        for i in 0..individual {
            groups.push(BaseGroup { start: i, end: i });
        }

        let mut pos = individual;
        while pos < max_length {
            let group_size = if pos < 50 {
                // 10-49: 5bp
                5
            } else if pos < 100 {
                // 50-99: 5bp
                5
            } else if pos < 500 {
                // 100-499: 5bp or 10bp depending on total len
                if max_length > 200 { 10 } else { 5 }
            } else {
                // 500+: 50bp
                50
            };

            let end = (pos + group_size - 1).min(max_length - 1);
            groups.push(BaseGroup { start: pos, end });
            pos = end + 1;
        }

        groups
    }
}

/// Phred encoding detection
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhredEncoding {
    Sanger,      // Phred+33
    Illumina1_3, // Phred+64
    Illumina1_5, // Phred+64
}

impl PhredEncoding {
    pub fn detect(lowest_char: u8) -> Self {
        if lowest_char < 64 {
            PhredEncoding::Sanger
        } else if lowest_char == 65 {
            PhredEncoding::Illumina1_3
        } else {
            PhredEncoding::Illumina1_5
        }
    }

    pub fn offset(&self) -> u8 {
        match self {
            PhredEncoding::Sanger => 33,
            PhredEncoding::Illumina1_3 | PhredEncoding::Illumina1_5 => 64,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            PhredEncoding::Sanger => "Sanger / Illumina 1.9",
            PhredEncoding::Illumina1_3 => "Illumina 1.3",
            PhredEncoding::Illumina1_5 => "Illumina 1.5",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_group_single_position() {
        let groups = BaseGroup::make_groups(5);
        assert_eq!(groups.len(), 5);
        assert_eq!(groups[0].label(), "1");
        assert_eq!(groups[4].label(), "5");
    }

    #[test]
    fn test_base_group_short_reads() {
        let groups = BaseGroup::make_groups(50);
        assert_eq!(groups.len(), 50); // all individual
    }

    #[test]
    fn test_base_group_long_reads() {
        let groups = BaseGroup::make_groups(150);
        assert!(groups.len() < 150); // should bin
        assert_eq!(groups[0].start, 0);
        assert_eq!(groups[0].end, 0); // first 9 individual
        assert_eq!(groups.last().unwrap().end, 149);
    }

    #[test]
    fn test_base_group_empty() {
        let groups = BaseGroup::make_groups(0);
        assert!(groups.is_empty());
    }

    #[test]
    fn test_phred_sanger() {
        assert_eq!(PhredEncoding::detect(33), PhredEncoding::Sanger);
        assert_eq!(PhredEncoding::detect(50), PhredEncoding::Sanger);
        assert_eq!(PhredEncoding::Sanger.offset(), 33);
    }

    #[test]
    fn test_phred_illumina() {
        assert_eq!(PhredEncoding::detect(65), PhredEncoding::Illumina1_3);
        assert_eq!(PhredEncoding::detect(66), PhredEncoding::Illumina1_5);
        assert_eq!(PhredEncoding::Illumina1_3.offset(), 64);
    }

    #[test]
    fn test_qc_result_labels() {
        assert_eq!(QCResult::Pass.label(), "PASS");
        assert_eq!(QCResult::Warn.label(), "WARN");
        assert_eq!(QCResult::Fail.label(), "FAIL");
    }
}

/// Factory to create all modules
pub struct ModuleFactory;

impl ModuleFactory {
    pub fn create_modules(config: &FastQCConfig) -> Vec<Box<dyn QCModule>> {
        let mut modules: Vec<Box<dyn QCModule>> = Vec::new();

        modules.push(Box::new(basic_stats::BasicStats::new()));

        if !config.is_ignored("quality_base") {
            modules.push(Box::new(per_base_quality::PerBaseQuality::new()));
        }
        if !config.is_ignored("tile") {
            modules.push(Box::new(per_tile_quality::PerTileQuality::new()));
        }
        if !config.is_ignored("quality_sequence") {
            modules.push(Box::new(per_sequence_quality::PerSequenceQuality::new()));
        }
        if !config.is_ignored("sequence") {
            modules.push(Box::new(per_base_content::PerBaseContent::new()));
        }
        if !config.is_ignored("gc_sequence") {
            modules.push(Box::new(per_sequence_gc::PerSequenceGC::new()));
        }
        if !config.is_ignored("n_content") {
            modules.push(Box::new(n_content::NContent::new()));
        }
        if !config.is_ignored("sequence_length") {
            modules.push(Box::new(sequence_length::SequenceLengthDist::new()));
        }
        if !config.is_ignored("duplication") {
            modules.push(Box::new(duplication::DuplicationLevel::new(config.dup_length)));
        }
        if !config.is_ignored("overrepresented") {
            modules.push(Box::new(overrepresented::OverrepresentedSeqs::new(config.dup_length)));
        }
        if !config.is_ignored("adapter") {
            modules.push(Box::new(adapter_content::AdapterContent::new(config)));
        }
        if !config.is_ignored("kmer") {
            modules.push(Box::new(kmer_content::KmerContent::new(config.kmer_size)));
        }

        // Long-read QC metrics
        if !config.is_ignored("read_length_n50") {
            modules.push(Box::new(long_read_quality::ReadLengthN50::new()));
        }
        if !config.is_ignored("quality_stratified_length") {
            modules.push(Box::new(long_read_quality::QualityStratifiedLength::new()));
        }
        if !config.is_ignored("homopolymer") {
            modules.push(Box::new(long_read_quality::HomopolymerErrors::new()));
        }

        modules
    }
}
