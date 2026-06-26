//! SOLiD colorspace to basespace decoder.
//!
//! SOLiD platform uses di-base encoding where each "color" (0-3) represents
//! a transition between two consecutive bases. The first character is a primer
//! base (A/C/G/T), followed by color codes.
//!
//! Transition table:
//!   Color 0: A→A, C→C, G→G, T→T (same base)
//!   Color 1: A→C, C→A, G→T, T→G
//!   Color 2: A→G, C→T, G→A, T→C
//!   Color 3: A→T, C→G, G→C, T→A

/// Decode a colorspace-encoded sequence (primer base + color digits) into basespace.
/// Returns None if the input is not valid colorspace.
pub fn decode_colorspace(csfasta: &[u8]) -> Option<Vec<u8>> {
    if csfasta.len() < 2 {
        return None;
    }

    // First character must be a base (the primer/adapter base)
    let mut prev_base = match csfasta[0] {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        _ => return None,
    };

    let mut basespace = Vec::with_capacity(csfasta.len());
    basespace.push(prev_base);

    for &color in &csfasta[1..] {
        let next = match (prev_base, color) {
            // Color 0: same base
            (b'A', b'0') => b'A',
            (b'C', b'0') => b'C',
            (b'G', b'0') => b'G',
            (b'T', b'0') => b'T',
            // Color 1
            (b'A', b'1') => b'C',
            (b'C', b'1') => b'A',
            (b'G', b'1') => b'T',
            (b'T', b'1') => b'G',
            // Color 2
            (b'A', b'2') => b'G',
            (b'C', b'2') => b'T',
            (b'G', b'2') => b'A',
            (b'T', b'2') => b'C',
            // Color 3
            (b'A', b'3') => b'T',
            (b'C', b'3') => b'G',
            (b'G', b'3') => b'C',
            (b'T', b'3') => b'A',
            // Period (.) or N means unknown color → N
            (_, b'.' | b'4' | b'5' | b'6') => b'N',
            _ => b'N',
        };
        basespace.push(next);
        if next != b'N' {
            prev_base = next;
        }
        // If we hit an N, keep prev_base unchanged so next decode is relative to last known base
    }

    Some(basespace)
}

/// Check if a sequence looks like colorspace encoding.
/// Colorspace reads start with a base letter followed by digits 0-3.
pub fn is_colorspace(seq: &[u8]) -> bool {
    if seq.len() < 2 {
        return false;
    }
    // First character must be a base
    let first_is_base = matches!(
        seq[0],
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
    );
    if !first_is_base {
        return false;
    }
    // Remaining characters should be color digits (0-3) or dots
    seq[1..]
        .iter()
        .all(|&b| matches!(b, b'0' | b'1' | b'2' | b'3' | b'.'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_simple() {
        // A followed by 0s = all A's
        let result = decode_colorspace(b"A0000").unwrap();
        assert_eq!(result, b"AAAAA");
    }

    #[test]
    fn test_decode_transitions() {
        // A→C (1), C→G (3), G→T (1), T→A (3)
        let result = decode_colorspace(b"A1313").unwrap();
        assert_eq!(result, b"ACGTA");
    }

    #[test]
    fn test_is_colorspace_true() {
        assert!(is_colorspace(b"A0123"));
        assert!(is_colorspace(b"T0000"));
        assert!(is_colorspace(b"G123."));
    }

    #[test]
    fn test_is_colorspace_false() {
        assert!(!is_colorspace(b"ACGT"));
        assert!(!is_colorspace(b"0123"));
        assert!(!is_colorspace(b""));
        assert!(!is_colorspace(b"A"));
    }

    #[test]
    fn test_decode_with_dot() {
        let result = decode_colorspace(b"A01.2").unwrap();
        // A→A(0), A→C(1), C→?(.) = N, keep prev=C, C→T(2)
        assert_eq!(result, b"AACNT");
    }
}
