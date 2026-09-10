use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;

/// Adapter entry: name and sequence
#[derive(Debug, Clone)]
pub struct Adapter {
    pub name: String,
    pub sequence: String,
}

/// Contaminant entry: name and sequence
#[derive(Debug, Clone)]
pub struct Contaminant {
    pub name: String,
    pub sequence: String,
    pub reverse_complement: String,
}

/// Limit thresholds for a module
#[derive(Debug, Clone)]
pub struct ModuleLimits {
    pub ignore: bool,
    pub warn: f64,
    pub error: f64,
}

/// Central configuration for all QC modules
#[derive(Debug, Clone)]
pub struct FastQCConfig {
    pub adapters: Vec<Adapter>,
    pub contaminants: Vec<Contaminant>,
    pub limits: HashMap<String, ModuleLimits>,
    pub kmer_size: usize,
    pub nofilter: bool,
    pub dup_length: usize,
    /// How many distinct sequences the observation-table modules
    /// (`DuplicationLevel`, `OverrepresentedSeqs`) may track before they stop
    /// admitting new ones. Defaults to [`OBSERVATION_BUDGET`].
    ///
    /// Exposed on the config so tests can shrink it; production runs use the
    /// full budget, because those modules run on a single instance that sees
    /// every read (see `QCModule::wants_all_reads`).
    pub observation_cutoff: usize,
}

/// Distinct sequences an observation-table module tracks before it freezes.
/// Matches FastQC's 100k limit.
pub const OBSERVATION_BUDGET: usize = 100_000;

const DEFAULT_ADAPTERS: &str = "\
Illumina Universal Adapter\tAGATCGGAAGAG
Illumina Small RNA 3' Adapter\tTGGAATTCTCGG
Illumina Small RNA 5' Adapter\tGATCGTCGGACT
Nextera Transposase Sequence\tCTGTCTCTTATA
PolyA\tAAAAAAAAAAAA
PolyG\tGGGGGGGGGGGG";

const DEFAULT_CONTAMINANTS: &str = "\
Illumina Single End Adapter 1\tGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
Illumina Single End Adapter 2\tCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
Illumina Single End PCR Primer 1\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Illumina Single End PCR Primer 2\tCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
Illumina Single End Sequencing Primer\tACACTCTTTCCCTACACGACGCTCTTCCGATCT
Illumina Paired End Adapter 1\tACACTCTTTCCCTACACGACGCTCTTCCGATCT
Illumina Paired End Adapter 2\tGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
Illumina Paried End PCR Primer 1\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Illumina Paired End PCR Primer 2\tCAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
Illumina Paried End Sequencing Primer 1\tACACTCTTTCCCTACACGACGCTCTTCCGATCT
Illumina Paired End Sequencing Primer 2\tCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
TruSeq Universal Adapter\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Nextera Transposase Sequence\tCTGTCTCTTATACACATCT
Illumina Multiplexing Adapter 1\tGATCGGAAGAGCACACGTCT
Illumina Multiplexing Adapter 2\tACACTCTTTCCCTACACGACGCTCTTCCGATCT";

impl FastQCConfig {
    pub fn new(
        contaminants_path: Option<&Path>,
        adapters_path: Option<&Path>,
        limits_path: Option<&Path>,
        kmer_size: usize,
        nofilter: bool,
        dup_length: usize,
    ) -> Result<Self> {
        let adapters = if let Some(path) = adapters_path {
            let content = std::fs::read_to_string(path)
                .with_context(|| format!("Failed to read adapters file: {}", path.display()))?;
            parse_adapter_list(&content)
        } else {
            parse_adapter_list(DEFAULT_ADAPTERS)
        };

        let contaminants = if let Some(path) = contaminants_path {
            let content = std::fs::read_to_string(path)
                .with_context(|| format!("Failed to read contaminants file: {}", path.display()))?;
            parse_contaminant_list(&content)
        } else {
            parse_contaminant_list(DEFAULT_CONTAMINANTS)
        };

        let limits = if let Some(path) = limits_path {
            let content = std::fs::read_to_string(path)
                .with_context(|| format!("Failed to read limits file: {}", path.display()))?;
            parse_limits(&content)
        } else {
            default_limits()
        };

        Ok(FastQCConfig {
            adapters,
            contaminants,
            limits,
            kmer_size,
            nofilter,
            dup_length,
            observation_cutoff: OBSERVATION_BUDGET,
        })
    }

    pub fn get_limit(&self, key: &str) -> Option<&ModuleLimits> {
        self.limits.get(key)
    }

    pub fn is_ignored(&self, module_key: &str) -> bool {
        self.limits
            .get(module_key)
            .map(|l| l.ignore)
            .unwrap_or(false)
    }

    /// Enable long-read QC modules (un-ignore them).
    pub fn enable_long_read_modules(&mut self) {
        for key in &[
            "read_length_n50",
            "quality_stratified_length",
            "homopolymer",
        ] {
            if let Some(lim) = self.limits.get_mut(*key) {
                lim.ignore = false;
            }
        }
    }
}

fn parse_adapter_list(content: &str) -> Vec<Adapter> {
    content
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                Some(Adapter {
                    name: parts[0].trim().to_string(),
                    sequence: parts.last().unwrap().trim().to_uppercase(),
                })
            } else {
                None
            }
        })
        .collect()
}

/// Complement a single IUPAC nucleotide code (uppercase input expected —
/// callers upper-case the sequence before calling this).
///
/// Contaminant/adapter lists occasionally use IUPAC ambiguity codes (e.g. `R`
/// for A-or-G, `K` for G-or-T) rather than only literal A/C/G/T/N. Before
/// this, any such code fell through to `other => other` and was passed
/// through unchanged instead of complemented, silently corrupting the
/// `reverse_complement` used for contaminant matching in
/// `OverrepresentedSeqs::find_contaminant` (e.g. `R` should complement to
/// `Y`, not stay `R`) — real sequences would then fail to match their
/// contaminant's reverse-complement orientation.
fn complement_base(c: char) -> char {
    match c {
        'A' => 'T',
        'T' | 'U' => 'A',
        'C' => 'G',
        'G' => 'C',
        'R' => 'Y', // A/G -> T/C
        'Y' => 'R', // C/T -> G/A
        'S' => 'S', // C/G -> G/C (self-complementary)
        'W' => 'W', // A/T -> T/A (self-complementary)
        'K' => 'M', // G/T -> C/A
        'M' => 'K', // A/C -> T/G
        'B' => 'V', // C/G/T -> G/C/A
        'V' => 'B', // A/C/G -> T/G/C
        'D' => 'H', // A/G/T -> T/C/A
        'H' => 'D', // A/C/T -> T/G/A
        'N' => 'N',
        other => other,
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_base).collect()
}

fn parse_contaminant_list(content: &str) -> Vec<Contaminant> {
    content
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                let seq = parts.last().unwrap().trim().to_uppercase();
                let rc = reverse_complement(&seq);
                Some(Contaminant {
                    name: parts[0].trim().to_string(),
                    sequence: seq,
                    reverse_complement: rc,
                })
            } else {
                None
            }
        })
        .collect()
}

fn parse_limits(content: &str) -> HashMap<String, ModuleLimits> {
    let mut limits: HashMap<String, ModuleLimits> = HashMap::new();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let key = parts[0].to_string();
            let level = parts[1];
            let value: f64 = match parts[2].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };

            let entry = limits.entry(key).or_insert(ModuleLimits {
                ignore: false,
                warn: f64::NAN,
                error: f64::NAN,
            });

            match level {
                "ignore" => entry.ignore = value != 0.0,
                "warn" => entry.warn = value,
                "error" => entry.error = value,
                _ => {}
            }
        }
    }

    limits
}

fn default_limits() -> HashMap<String, ModuleLimits> {
    let content = "\
duplication\tignore\t0
kmer\tignore\t0
n_content\tignore\t0
overrepresented\tignore\t0
quality_base\tignore\t0
sequence\tignore\t0
gc_sequence\tignore\t0
quality_sequence\tignore\t0
tile\tignore\t0
sequence_length\tignore\t0
adapter\tignore\t0
duplication\twarn\t70
duplication\terror\t50
kmer\twarn\t2
kmer\terror\t5
n_content\twarn\t5
n_content\terror\t20
overrepresented\twarn\t0.1
overrepresented\terror\t1
quality_base_lower\twarn\t10
quality_base_lower\terror\t5
quality_base_median\twarn\t25
quality_base_median\terror\t20
sequence\twarn\t10
sequence\terror\t20
gc_sequence\twarn\t15
gc_sequence\terror\t30
quality_sequence\twarn\t27
quality_sequence\terror\t20
tile\twarn\t5
tile\terror\t10
sequence_length\twarn\t1
sequence_length\terror\t1
adapter\twarn\t5
adapter\terror\t10
read_length_n50\tignore\t1
read_length_n50\twarn\t0
read_length_n50\terror\t0
quality_stratified_length\tignore\t1
quality_stratified_length\twarn\t0
quality_stratified_length\terror\t0
homopolymer\tignore\t1
homopolymer\twarn\t5
homopolymer\terror\t10";

    parse_limits(content)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = FastQCConfig::new(None, None, None, 7, false, 50).unwrap();
        assert_eq!(config.kmer_size, 7);
        assert_eq!(config.dup_length, 50);
        assert!(!config.adapters.is_empty());
        assert!(!config.contaminants.is_empty());
    }

    #[test]
    fn test_default_adapters_count() {
        let config = FastQCConfig::new(None, None, None, 7, false, 50).unwrap();
        assert_eq!(config.adapters.len(), 6); // 6 default adapters
    }

    #[test]
    fn test_adapter_parsing() {
        let content = "Adapter1\tACTG\nAdapter2\tGGGG\n# comment\n";
        let adapters = parse_adapter_list(content);
        assert_eq!(adapters.len(), 2);
        assert_eq!(adapters[0].name, "Adapter1");
        assert_eq!(adapters[0].sequence, "ACTG");
    }

    #[test]
    fn test_contaminant_revcomp() {
        let content = "Test\tACGT\n";
        let contaminants = parse_contaminant_list(content);
        assert_eq!(contaminants[0].reverse_complement, "ACGT"); // palindrome
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("AATTCCGG"), "CCGGAATT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
    }

    #[test]
    fn reverse_complement_handles_iupac_ambiguity_codes() {
        // Regression test: ambiguity codes used to fall through to
        // `other => other` (passed through unchanged) instead of being
        // complemented, silently corrupting reverse-complement contaminant
        // matching for any contaminant/adapter sequence using them.
        assert_eq!(reverse_complement("R"), "Y"); // A/G -> complement C/T
        assert_eq!(reverse_complement("Y"), "R");
        assert_eq!(reverse_complement("K"), "M");
        assert_eq!(reverse_complement("M"), "K");
        assert_eq!(reverse_complement("S"), "S"); // self-complementary
        assert_eq!(reverse_complement("W"), "W"); // self-complementary
        assert_eq!(reverse_complement("B"), "V");
        assert_eq!(reverse_complement("D"), "H");
        // Full round trip on a mixed IUPAC sequence, read in reverse order.
        assert_eq!(reverse_complement("ACGTRYSWKMBDHVN"), "NBDHVKMWSRYACGT");
    }

    #[test]
    fn test_limits_parsing() {
        let content = "quality_base_lower\twarn\t10\nquality_base_lower\terror\t5\n";
        let limits = parse_limits(content);
        let qbl = limits.get("quality_base_lower").unwrap();
        assert_eq!(qbl.warn, 10.0);
        assert_eq!(qbl.error, 5.0);
        assert!(!qbl.ignore);
    }

    #[test]
    fn test_ignore_module() {
        let content = "kmer\tignore\t1\n";
        let limits = parse_limits(content);
        assert!(limits.get("kmer").unwrap().ignore);
    }

    #[test]
    fn test_is_ignored() {
        let config = FastQCConfig::new(None, None, None, 7, false, 50).unwrap();
        assert!(!config.is_ignored("kmer")); // kmer enabled by default
        assert!(!config.is_ignored("quality_base"));
    }
}
