//! Translated from `fgsea/src/esCalculation.cpp` + `fgsea/src/esCalculation.h`.
//!
//! `int128` (boost::multiprecision or `__int128`) maps to native `i128`.
//!
//! Vendored from the faithful fgsea-rs translation.

/// An exact (rational) enrichment-score value, where
/// `score = coef_NS / NS - coef_const / diff`.
///
/// `struct score_t`
#[derive(Clone, Copy, Debug)]
pub struct score_t {
    pub NS: i64,
    pub coef_NS: i64,
    pub diff: i64, // n - k
    pub coef_const: i64,
}

impl score_t {
    /// Returns the score as a floating-point value (`coef_NS / NS - coef_const / diff`).
    ///
    /// `double getDouble() const`
    pub fn getDouble(&self) -> f64 {
        1.0 * self.coef_NS as f64 / self.NS as f64 - 1.0 * self.coef_const as f64 / self.diff as f64
    }

    /// Returns the numerator of the score (`coef_NS * diff - coef_const * NS`).
    ///
    /// `int64_t getNumerator() const`
    pub fn getNumerator(&self) -> i64 {
        self.coef_NS * self.diff - self.coef_const * self.NS
    }

    /// Compares this score with `other` exactly, returning a 128-bit signed value
    /// whose sign gives the ordering (cross-multiplied to avoid division).
    ///
    /// `int128 Compare(score_t const& other) const`
    pub fn Compare(&self, other: &score_t) -> i128 {
        let P1: i64 = self.coef_NS * self.diff + self.NS * (other.coef_const - self.coef_const);
        let Q1: i64 = self.NS * self.diff;
        let P2: i64 = other.coef_NS;
        let Q2: i64 = other.NS;
        (P1 as i128) * (Q2 as i128) - (P2 as i128) * (Q1 as i128)
    }

    /// Returns the maximum supported `NS` value (`1 << 30`).
    ///
    /// `static int64_t getMaxNS()`
    pub fn getMaxNS() -> i64 {
        1i64 << 30
    }

    /// Returns the absolute value of the score, i.e. `std::max(*this, -(*this))`.
    ///
    /// `score_t abs() const`
    pub fn abs(&self) -> score_t {
        (*self).max(-(*self))
    }
}

/// Negates the score by negating `coef_NS` and `coef_const`.
///
/// `score_t operator-() const`
impl std::ops::Neg for score_t {
    type Output = score_t;
    fn neg(self) -> score_t {
        score_t { NS: self.NS, coef_NS: -self.coef_NS, diff: self.diff, coef_const: -self.coef_const }
    }
}

// The six comparison operators all delegate to `Compare`, exactly as in C++.
impl PartialEq for score_t {
    fn eq(&self, other: &Self) -> bool {
        self.Compare(other) == 0
    }
}
impl Eq for score_t {}
impl PartialOrd for score_t {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for score_t {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // C++ `score_t::operator<` is exactly `Compare(other) < 0` — NOT a total
        // order. When `diff` differs (e.g. the diff=1 ES probe vs a diff=(n-k)
        // sample), `Compare` is not antisymmetric: it can report NEITHER a<b nor
        // b<a. `std::pair<score_t,hash_t>` (gsea_t) then treats such scores as
        // EQUIVALENT and falls through to the hash tiebreaker. Reproduce that by
        // returning `Equal` in the "neither is less" case (returning `Greater`
        // here, as a naive `Compare.cmp(&0)` would, breaks the pair tiebreak and
        // diverges from the C++ p-values).
        if self.Compare(other) < 0 {
            std::cmp::Ordering::Less
        } else if other.Compare(self) < 0 {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

/// Computes the (two-sided) enrichment score of the gene set `p` over `ranks`,
/// with `NS` (the sum of selected ranks) supplied by the caller.
///
/// `p` must be sorted.
///
/// `score_t calcES(const vector<int64_t>&, const vector<int>&, int64_t NS)`
pub fn calcES(ranks: &[i64], p: &[i32], NS: i64) -> score_t {
    // p must be sorted
    let n: i32 = ranks.len() as i32;
    let k: i32 = p.len() as i32;
    let mut res = score_t { NS, coef_NS: 0, diff: (n - k) as i64, coef_const: 0 };
    let mut cur = score_t { NS, coef_NS: 0, diff: (n - k) as i64, coef_const: 0 };
    let mut last: i32 = -1;
    for &pos in p {
        cur.coef_const += (pos - last - 1) as i64;
        if res.abs() < cur.abs() {
            res = cur;
        }
        cur.coef_NS += ranks[pos as usize];
        if res.abs() < cur.abs() {
            res = cur;
        }
        last = pos;
    }
    res
}

/// Computes the (two-sided) enrichment score of the gene set `p` over `ranks`,
/// computing `NS` internally as the sum of the selected ranks.
///
/// `p` must be sorted.
///
/// `score_t calcES(const vector<int64_t>&, const vector<int>&)`
pub fn calcES_auto(ranks: &[i64], p: &[i32]) -> score_t {
    // p must be sorted
    let mut NS: i64 = 0;
    for &pos in p {
        NS += ranks[pos as usize];
    }
    calcES(ranks, p, NS)
}

/// Computes the positive (one-sided) enrichment score of the gene set `p` over
/// `ranks`, with `NS` supplied by the caller.
///
/// `p` must be sorted.
///
/// `score_t calcPositiveES(const vector<int64_t>&, const vector<int>&, int64_t NS)`
pub fn calcPositiveES(ranks: &[i64], p: &[i32], NS: i64) -> score_t {
    // p must be sorted
    let n: i32 = ranks.len() as i32;
    let k: i32 = p.len() as i32;
    let mut res = score_t { NS, coef_NS: 0, diff: (n - k) as i64, coef_const: 0 };
    let mut cur = score_t { NS, coef_NS: 0, diff: (n - k) as i64, coef_const: 0 };
    let mut last: i32 = -1;
    for &pos in p {
        cur.coef_NS += ranks[pos as usize];
        cur.coef_const += (pos - last - 1) as i64;
        res = res.max(cur);
        last = pos;
    }
    res
}

/// Computes the positive (one-sided) enrichment score of the gene set `p` over
/// `ranks`, computing `NS` internally as the sum of the selected ranks.
///
/// `p` must be sorted.
///
/// `score_t calcPositiveES(const vector<int64_t>&, const vector<int>&)`
pub fn calcPositiveES_auto(ranks: &[i64], p: &[i32]) -> score_t {
    // p must be sorted
    let mut NS: i64 = 0;
    for &pos in p {
        NS += ranks[pos as usize];
    }
    calcPositiveES(ranks, p, NS)
}

/// Returns `true` as soon as the running enrichment statistic of the gene set
/// `p` exceeds `bound` (a fast early-exit check), otherwise `false`.
///
/// `p` must be sorted.
///
/// `bool compareStat(const vector<double>&, const vector<int>&, double NS, double bound)`
pub fn compareStat(ranks: &[f64], p: &[i32], NS: f64, bound: f64) -> bool {
    // p must be sorted
    let n: i32 = ranks.len() as i32;
    let k: i32 = p.len() as i32;
    let mut cur: f64 = 0.0;
    let q1: f64 = 1.0 / (n - k) as f64;
    let q2: f64 = 1.0 / NS;
    let mut last: i32 = -1;
    for &pos in p {
        cur += q2 * ranks[pos as usize] - q1 * (pos - last - 1) as f64;
        if cur > bound {
            return true;
        }
        last = pos;
    }
    false
}
