use crate::config::FastQCConfig;
use crate::io::Sequence;
use super::{PhredEncoding, QCModule, QCResult};

pub struct BasicStats {
    total_sequences: u64,
    filtered_sequences: u64,
    min_length: usize,
    max_length: usize,
    total_bases: u64,
    a_count: u64,
    t_count: u64,
    g_count: u64,
    c_count: u64,
    n_count: u64,
    lowest_char: u8,
    encoding: Option<PhredEncoding>,
    gc_percent: f64,
}

impl BasicStats {
    pub fn new() -> Self {
        BasicStats {
            total_sequences: 0,
            filtered_sequences: 0,
            min_length: usize::MAX,
            max_length: 0,
            total_bases: 0,
            a_count: 0,
            t_count: 0,
            g_count: 0,
            c_count: 0,
            n_count: 0,
            lowest_char: 255,
            encoding: None,
            gc_percent: 0.0,
        }
    }

    #[allow(dead_code)]
    pub fn total_sequences(&self) -> u64 {
        self.total_sequences
    }

    #[allow(dead_code)]
    pub fn min_length(&self) -> usize {
        if self.min_length == usize::MAX { 0 } else { self.min_length }
    }

    #[allow(dead_code)]
    pub fn max_length(&self) -> usize {
        self.max_length
    }

    #[allow(dead_code)]
    pub fn encoding(&self) -> Option<PhredEncoding> {
        self.encoding
    }
}

impl QCModule for BasicStats {
    fn name(&self) -> &str {
        "Basic Statistics"
    }

    fn key(&self) -> &str {
        "basic_statistics"
    }

    fn process_sequence(&mut self, seq: &Sequence) {
        self.total_sequences += 1;
        if seq.filtered {
            self.filtered_sequences += 1;
        }

        let len = seq.len();
        if len < self.min_length {
            self.min_length = len;
        }
        if len > self.max_length {
            self.max_length = len;
        }
        self.total_bases += len as u64;

        for &b in &seq.sequence {
            match b {
                b'A' | b'a' => self.a_count += 1,
                b'T' | b't' => self.t_count += 1,
                b'G' | b'g' => self.g_count += 1,
                b'C' | b'c' => self.c_count += 1,
                _ => self.n_count += 1,
            }
        }

        for &q in &seq.quality {
            if q < self.lowest_char {
                self.lowest_char = q;
            }
        }
    }

    fn calculate_results(&mut self, _config: &FastQCConfig) {
        if self.lowest_char < 255 {
            self.encoding = Some(PhredEncoding::detect(self.lowest_char));
        }
        let atgc = self.a_count + self.t_count + self.g_count + self.c_count;
        if atgc > 0 {
            self.gc_percent = ((self.g_count + self.c_count) as f64 / atgc as f64) * 100.0;
        }
    }

    fn result(&self) -> QCResult {
        QCResult::Pass // Basic stats is always informational
    }

    fn has_chart(&self) -> bool {
        false
    }

    fn text_data(&self) -> String {
        let length_str = if self.min_length() == self.max_length {
            format!("{}", self.max_length)
        } else {
            format!("{}-{}", self.min_length(), self.max_length)
        };

        let encoding_str = self
            .encoding
            .map(|e| e.name().to_string())
            .unwrap_or_else(|| "Unknown".to_string());

        let total_bases_str = if self.total_bases >= 1_000_000_000 {
            format!("{:.1} Gbp", self.total_bases as f64 / 1_000_000_000.0)
        } else if self.total_bases >= 1_000_000 {
            format!("{} Mbp", self.total_bases / 1_000_000)
        } else {
            format!("{}", self.total_bases)
        };

        format!(
            ">>Basic Statistics\t{}\n\
             #Measure\tValue\n\
             Filename\t{{}}\n\
             File type\tConventional base calls\n\
             Encoding\t{}\n\
             Total Sequences\t{}\n\
             Total Bases\t{}\n\
             Sequences flagged as poor quality\t{}\n\
             Sequence length\t{}\n\
             %GC\t{}\n\
             >>END_MODULE\n",
            self.result().label(),
            encoding_str,
            self.total_sequences,
            total_bases_str,
            self.filtered_sequences,
            length_str,
            self.gc_percent.round() as u64
        )
    }

    fn svg_chart(&self) -> String {
        String::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::FastQCConfig;

    fn make_seq(seq: &[u8], qual: &[u8]) -> Sequence {
        Sequence {
            header: "@test".to_string(),
            sequence: seq.to_vec(),
            quality: qual.to_vec(),
            filtered: false,
        }
    }

    fn default_config() -> FastQCConfig {
        FastQCConfig::new(None, None, None, 7, false, 50).unwrap()
    }

    #[test]
    fn test_counts_sequences() {
        let mut m = BasicStats::new();
        let config = default_config();
        m.process_sequence(&make_seq(b"ACTG", b"IIII"));
        m.process_sequence(&make_seq(b"ACTG", b"IIII"));
        m.calculate_results(&config);
        assert_eq!(m.total_sequences(), 2);
    }

    #[test]
    fn test_gc_percent() {
        let mut m = BasicStats::new();
        let config = default_config();
        // 50% GC
        m.process_sequence(&make_seq(b"AACCGGTT", b"IIIIIIII"));
        m.calculate_results(&config);
        assert!((m.gc_percent - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_length_range() {
        let mut m = BasicStats::new();
        let config = default_config();
        m.process_sequence(&make_seq(b"ACTG", b"IIII"));
        m.process_sequence(&make_seq(b"ACTGACTG", b"IIIIIIII"));
        m.calculate_results(&config);
        assert_eq!(m.min_length(), 4);
        assert_eq!(m.max_length(), 8);
    }

    #[test]
    fn test_encoding_sanger() {
        let mut m = BasicStats::new();
        let config = default_config();
        m.process_sequence(&make_seq(b"ACTG", b"!@#$")); // lowest is '!' = 33
        m.calculate_results(&config);
        assert_eq!(m.encoding(), Some(PhredEncoding::Sanger));
    }

    #[test]
    fn test_text_data_format() {
        let mut m = BasicStats::new();
        let config = default_config();
        m.process_sequence(&make_seq(b"ACTG", b"IIII"));
        m.calculate_results(&config);
        let text = m.text_data();
        assert!(text.starts_with(">>Basic Statistics\tPASS"));
        assert!(text.contains("Total Sequences\t1"));
        assert!(text.ends_with(">>END_MODULE\n"));
    }
}
