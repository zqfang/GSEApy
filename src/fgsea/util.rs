//! Translated from `fgsea/src/util.cpp` + `fgsea/src/util.h`.
//!
//! Only the active (`#ifndef USE_STD_UID`) branch of `uid_wrapper` is translated;
//! the `USE_STD_UID` branch is dead under the default build.
//!
//! Vendored from the faithful fgsea-rs translation (boost::mt19937 bit-for-bit).

use special::Gamma;

/// The Mersenne Twister random engine used throughout the multilevel code.
///
/// `rand_mt::Mt` reproduces `boost::mt19937` bit-for-bit when seeded with a single
/// `u32`, preserving reproducibility with the original C++.
///
/// C++: `using random_engine_t = boost::mt19937;`
pub type random_engine_t = rand_mt::Mt;

/// Draws a uniform integer from the closed interval `[from, from + len)` using
/// rejection sampling on the engine's raw output, avoiding modulo bias.
///
/// This mirrors the active (`#ifndef USE_STD_UID`) branch of the C++ `uid_wrapper`:
/// values at or above `completePart` are rejected so that the remaining range divides
/// evenly by `len`.
///
/// Unlike the C++ struct it does not store the RNG reference (Rust aliasing rules);
/// the RNG is threaded through [`uid_wrapper::call`] instead. `completePart` is
/// precomputed from `rng.max()` (a compile-time constant `u32::MAX` for mt19937) in
/// [`uid_wrapper::new`], exactly as the C++ constructor does.
///
/// C++: `struct uid_wrapper`
pub struct uid_wrapper {
    pub from: i32,
    pub len: i32,
    pub completePart: u32,
}

impl uid_wrapper {
    /// Builds a sampler for the closed interval `[_from, _to]`, precomputing the
    /// rejection threshold `completePart` from the engine's maximum value.
    ///
    /// C++: `uid_wrapper::uid_wrapper(int _from, int _to, random_engine_t& _rng)`
    pub fn new(_from: i32, _to: i32, _rng: &random_engine_t) -> uid_wrapper {
        let from = _from;
        let len = _to - _from + 1;
        // unsigned maxVal = rng.max();  (== u32::MAX for boost::mt19937)
        let maxVal: u32 = u32::MAX;
        let completePart = maxVal - maxVal % (len as u32);
        uid_wrapper {
            from,
            len,
            completePart,
        }
    }

    /// Draws one sample, repeatedly pulling from `rng` until the raw value falls
    /// below `completePart`, then mapping it into `[from, from + len)`.
    ///
    /// C++: `int uid_wrapper::operator()()`
    pub fn call(&self, rng: &mut random_engine_t) -> i32 {
        let mut x: u32;
        loop {
            x = rng.next_u32();
            if x < self.completePart {
                break;
            }
        }
        self.from + (x % (self.len as u32)) as i32
    }
}

/// Generates `k` distinct numbers from the closed interval `[a, b]`.
///
/// `a` should be non-negative, usually 0 or 1. When `k` is small relative to the
/// range, numbers are drawn with rejection until `k` distinct values are collected
/// (on average fewer than two draws each); otherwise a selection-plus-shuffle
/// strategy is used. The Fisher-Yates (Knuth) shuffle is implemented explicitly for
/// reproducibility between platforms.
///
/// # Arguments
/// * `a` - inclusive lower bound of the interval (non-negative).
/// * `b` - inclusive upper bound of the interval.
/// * `k` - number of distinct values to draw.
/// * `rng` - random engine.
///
/// # Returns
/// A vector of `k` distinct integers in `[a, b]`.
///
/// C++: `std::vector<int> combination(const int&, const int&, const int&, random_engine_t&)`
pub fn combination(a: i32, b: i32, k: i32, rng: &mut random_engine_t) -> Vec<i32> {
    let uni = uid_wrapper::new(a, b, rng);
    let mut v: Vec<i32> = Vec::with_capacity(k as usize);

    let n = b - a + 1;
    let mut used = vec![0u8; n as usize]; // std::vector<char> used(n);

    if (k as f64) < n as f64 * 1.0 / 2.0 {
        for _i in 0..k {
            // average < 2 iterations
            loop {
                let x = uni.call(rng);
                if used[(x - a) as usize] == 0 {
                    v.push(x);
                    used[(x - a) as usize] = 1;
                    break;
                }
            }
        }
    } else {
        for r in (n - k)..n {
            let x = uid_wrapper::new(0, r, rng).call(rng);
            if used[x as usize] == 0 {
                v.push(a + x);
                used[x as usize] = 1;
            } else {
                v.push(a + r);
                used[r as usize] = 1;
            }
        }

        // implemented explicitly for reproducibility between platforms
        // Fisher–Yates (Knuth) shuffle
        let mut i = v.len() as i32 - 1;
        while i > 0 {
            // pick j in [0..i]
            let j = uid_wrapper::new(0, i, rng).call(rng);
            v.swap(i as usize, j as usize);
            i -= 1;
        }
    }

    v
}

/// Mean of the log of a Beta(`a`, `b + 1 - a`) variate, computed as
/// `digamma(a) - digamma(b + 1)`.
///
/// C++: `double betaMeanLog(unsigned long a, unsigned long b)`
pub fn betaMeanLog(a: u64, b: u64) -> f64 {
    (a as f64).digamma() - ((b + 1) as f64).digamma()
}

/// Estimates the standard error (in log2 units) accumulated over `level` multilevel
/// steps, derived from the per-level variance for a given `sampleSize`.
///
/// C++: `double multilevelError(int level, int sampleSize)`
pub fn multilevelError(level: i32, sampleSize: i32) -> f64 {
    let singleLevelError =
        (((sampleSize + 1) / 2) as f64).trigamma() - ((sampleSize + 1) as f64).trigamma();
    ((level as f64) * singleLevelError).sqrt() / (2.0_f64).ln()
}
