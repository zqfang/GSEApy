/// A self-contained, portable Xoshiro256++ PRNG implementation.
///
/// This implementation uses the xoshiro256++ algorithm from David Blackman and
/// Sebastiano Vigna. The algorithm is explicitly documented and fixed, ensuring
/// reproducible results regardless of the `rand` crate version.
///
/// # Motivation
///
/// The `rand::rngs::SmallRng` type is explicitly documented as non-portable:
/// its underlying algorithm may change between `rand` major versions.
/// Specifically, `SmallRng` changed from an XorShift-based algorithm in
/// rand 0.7.x to `Xoshiro256PlusPlus` in rand 0.8.x and later, which caused
/// non-reproducible GSEA results for users who had different `rand` versions
/// compiled into their gseapy installation.
///
/// This self-contained implementation guarantees stable, reproducible output
/// for the same seed across all future versions of this library.
use rand::{RngCore, SeedableRng};

/// Xoshiro256++ random number generator.
///
/// Algorithm from: http://xoshiro.di.unimi.it/xoshiro256plusplus.c
/// by David Blackman and Sebastiano Vigna.
#[derive(Debug, Clone)]
pub struct Xoshiro256PlusPlus {
    s: [u64; 4],
}

impl SeedableRng for Xoshiro256PlusPlus {
    type Seed = [u8; 32];

    /// Create a new `Xoshiro256PlusPlus` from a 32-byte seed.
    #[inline]
    fn from_seed(seed: [u8; 32]) -> Self {
        let mut s = [0u64; 4];
        for (i, chunk) in seed.chunks_exact(8).enumerate() {
            s[i] = u64::from_le_bytes(chunk.try_into().unwrap());
        }
        // If all-zero, use seed_from_u64(0) to avoid the degenerate state.
        if s.iter().all(|&x| x == 0) {
            return Self::seed_from_u64(0);
        }
        Xoshiro256PlusPlus { s }
    }

    /// Create a new `Xoshiro256PlusPlus` from a `u64` seed using SplitMix64.
    ///
    /// This uses the same SplitMix64 expansion as the reference implementation
    /// and as `rand` 0.8/0.9's `SmallRng::seed_from_u64` on 64-bit platforms,
    /// ensuring identical output for the same seed.
    #[inline]
    fn seed_from_u64(mut state: u64) -> Self {
        const PHI: u64 = 0x9e3779b97f4a7c15;
        let mut s = [0u64; 4];
        for i in s.iter_mut() {
            state = state.wrapping_add(PHI);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^= z >> 31;
            *i = z;
        }
        Xoshiro256PlusPlus { s }
    }
}

impl RngCore for Xoshiro256PlusPlus {
    #[inline]
    fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let res = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);

        let t = self.s[1] << 17;

        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];

        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);

        res
    }

    #[inline]
    fn fill_bytes(&mut self, dst: &mut [u8]) {
        let mut left = dst;
        while left.len() >= 8 {
            let bytes = self.next_u64().to_le_bytes();
            let (block, rest) = left.split_at_mut(8);
            block.copy_from_slice(&bytes);
            left = rest;
        }
        if !left.is_empty() {
            let bytes = self.next_u64().to_le_bytes();
            left.copy_from_slice(&bytes[..left.len()]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Xoshiro256PlusPlus;
    use rand::{RngCore, SeedableRng};

    /// Reference values from http://xoshiro.di.unimi.it/xoshiro256plusplus.c
    #[test]
    fn test_reference_values() {
        let mut rng = Xoshiro256PlusPlus::from_seed([
            1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0,
            0, 0, 0, 0,
        ]);
        let expected: [u64; 10] = [
            41943041,
            58720359,
            3588806011781223,
            3591011842654386,
            9228616714210784205,
            9973669472204895162,
            14011001112246962877,
            12406186145184390807,
            15849039046786891736,
            10450023813501588000,
        ];
        for &e in &expected {
            assert_eq!(rng.next_u64(), e);
        }
    }

    /// Verify that seed_from_u64 produces the same sequence as rand 0.8/0.9 SmallRng
    /// on 64-bit platforms, ensuring backward compatibility for existing seeds.
    #[test]
    fn test_seed_from_u64_stability() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        let expected: [u64; 10] = [
            5987356902031041503,
            7051070477665621255,
            6633766593972829180,
            211316841551650330,
            9136120204379184874,
            379361710973160858,
            15813423377499357806,
            15596884590815070553,
            5439680534584881407,
            1369371744833522710,
        ];
        for &e in &expected {
            assert_eq!(rng.next_u64(), e);
        }
    }
}
