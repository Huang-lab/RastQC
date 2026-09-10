//! A fast, non-cryptographic hasher for the modules' internal lookup tables.
//!
//! The observation tables (`DuplicationLevel`, `OverrepresentedSeqs`) and the
//! kmer table are all keyed by short byte strings and are probed once or more
//! per read, which makes hashing one of the hottest things RastQC does.
//! `std`'s default `SipHash-1-3` is chosen to be resistant to
//! hash-flooding by untrusted keys; these tables are keyed by base calls from
//! a local sequencing file, and RastQC exposes no service where an attacker
//! chooses those keys, so that resistance buys nothing here and costs roughly
//! an order of magnitude per probe.
//!
//! This is the same FxHash construction rustc uses for its own internal maps:
//! multiply-and-rotate per word, no finalization.

use std::hash::{BuildHasherDefault, Hasher};

/// Drop-in replacement for `RandomState` on the modules' internal maps.
pub type FxBuildHasher = BuildHasherDefault<FxHasher>;

/// Fractional part of the golden ratio scaled to 64 bits — the multiplier
/// spreads each input word's bits across the whole accumulator.
const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;
const ROTATE: u32 = 5;

#[derive(Default)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    #[inline]
    fn add_to_hash(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(ROTATE) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut rest = bytes;
        // Consume 8 bytes at a time, then whatever is left. Kmer and
        // truncated-read keys are short (7-50 bytes), so this is a handful of
        // multiplies per probe.
        while rest.len() >= 8 {
            let (word, tail) = rest.split_at(8);
            self.add_to_hash(u64::from_le_bytes(word.try_into().unwrap()));
            rest = tail;
        }
        if !rest.is_empty() {
            let mut buf = [0u8; 8];
            buf[..rest.len()].copy_from_slice(rest);
            self.add_to_hash(u64::from_le_bytes(buf));
        }
        // Fold the length in so that keys differing only by trailing zero
        // bytes (which the padding above would otherwise collapse) still
        // hash apart.
        self.add_to_hash(bytes.len() as u64);
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add_to_hash(i as u64);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add_to_hash(i as u64);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add_to_hash(i);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add_to_hash(i as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::hash::Hash;

    fn hash_of(bytes: &[u8]) -> u64 {
        let mut h = FxHasher::default();
        h.write(bytes);
        h.finish()
    }

    #[test]
    fn equal_keys_hash_equally() {
        assert_eq!(hash_of(b"ACGTACG"), hash_of(b"ACGTACG"));
    }

    #[test]
    fn keys_differing_only_in_length_hash_apart() {
        // The 8-byte padding would map "ACG" and "ACG\0" to the same word, so
        // the length must participate in the hash.
        assert_ne!(hash_of(b"ACG"), hash_of(b"ACG\0"));
        assert_ne!(hash_of(b""), hash_of(b"\0"));
    }

    #[test]
    fn distinct_kmers_mostly_land_on_distinct_hashes() {
        // Every 7-mer over ACGT: 16384 keys, which is exactly the table the
        // kmer module builds. Collisions are allowed but should be rare
        // enough that the map stays O(1) in practice.
        let bases = b"ACGT";
        let mut hashes = std::collections::HashSet::new();
        let mut count = 0;
        for i in 0..4usize.pow(7) {
            let mut kmer = [0u8; 7];
            let mut n = i;
            for slot in kmer.iter_mut() {
                *slot = bases[n % 4];
                n /= 4;
            }
            hashes.insert(hash_of(&kmer));
            count += 1;
        }
        assert_eq!(count, 16384);
        assert!(
            hashes.len() > 16_300,
            "expected near-perfect spread over 7-mers, got {} distinct hashes",
            hashes.len()
        );
    }

    #[test]
    fn works_as_a_hashmap_build_hasher() {
        let mut map: HashMap<Vec<u8>, u32, FxBuildHasher> = HashMap::default();
        for i in 0..1000u32 {
            map.insert(format!("seq{i}").into_bytes(), i);
        }
        assert_eq!(map.len(), 1000);
        for i in 0..1000u32 {
            assert_eq!(map.get(format!("seq{i}").as_bytes()), Some(&i));
        }
        // Borrowed lookups must keep working, which requires `Vec<u8>` and
        // `[u8]` to hash identically through this hasher.
        let key = b"seq500".to_vec();
        let mut a = FxHasher::default();
        key.hash(&mut a);
        let mut b = FxHasher::default();
        key[..].hash(&mut b);
        assert_eq!(a.finish(), b.finish());
    }
}
