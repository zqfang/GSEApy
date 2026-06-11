//! Translated from `fgsea/src/fastGSEA.cpp` + `fgsea/src/fastGSEA.h`.
//!
//! C++ templates become Rust generics:
//!   * `SegmentTree<T>` -> `SegmentTree<T>`
//!   * `order<T>` -> `order<T>` (the `IndirectCmp<T>` functor is folded into
//!     `order`'s sort comparator — it has no standalone Rust counterpart).
//!
//! The Rcpp return aggregates become plain structs:
//!   * `calcGseaStatCumulativeBatch` -> `GseaStatCumulativeBatch`
//!
//! Vendored from the faithful fgsea-rs translation.

use crate::fgsea::util::{combination, random_engine_t};

pub const eps: f64 = 1e-13;

/// "Not actually a segment tree, but a square root heuristic" — a data
/// structure supporting point increments and prefix-sum queries over `[0, n)`.
///
/// `inc` runs in `O(sqrt(n))` and `queryR` in `O(1)`.
///
/// `template <class T> class SegmentTree`
pub struct SegmentTree<T> {
    t: Vec<T>,
    b: Vec<T>,
    n: i32,
    k: i32,
    k2: i32,
    logK: i32,
    blockMask: i32,
}

impl<T> SegmentTree<T>
where
    T: Copy + Default + std::ops::AddAssign + std::ops::Add<Output = T>,
{
    /// Constructs a tree sized for `n_` positions.
    ///
    /// `SegmentTree(int n_)`
    pub fn new(n_: i32) -> SegmentTree<T> {
        let mut k = 1;
        let mut logK = 0;
        while k * k < n_ {
            k <<= 1;
            logK += 1;
        }
        let k2 = (n_ - 1) / k + 1;
        let n = k * k;

        let blockMask = k - 1;

        let t = vec![T::default(); n as usize];
        let b = vec![T::default(); k2 as usize];

        SegmentTree {
            t,
            b,
            n,
            k,
            k2,
            logK,
            blockMask,
        }
    }

    /// Increases the value at position `p` by `delta`. Works in `O(sqrt(n))`.
    ///
    /// `void inc(int p, T delta)`
    pub fn inc(&mut self, p: i32, delta: T) {
        let mut p = p;
        let blockEnd = p - (p & self.blockMask) + self.blockMask + 1;
        while p < blockEnd {
            self.t[p as usize] += delta;
            p += 1;
        }

        let mut p1 = p >> self.logK;
        while p1 < self.k2 {
            self.b[p1 as usize] += delta;
            p1 += 1;
        }
    }

    /// Returns the sum over the interval `[0, r)`. Works in `O(1)`.
    ///
    /// `T queryR(int r)`
    pub fn queryR(&self, r: i32) -> T {
        if r == 0 {
            return T::default();
        }
        let r = r - 1;
        self.t[r as usize] + self.b[(r >> self.logK) as usize]
    }
}

/// Returns the indices that sort `x` in ascending order (an "order" / argsort).
///
/// Folds in the C++ `IndirectCmp<T>` comparator, which compared elements
/// indirectly through `x`.
///
/// `template<class T> vector<int> order(T const& x)`
pub fn order<T: PartialOrd>(x: &[T]) -> Vec<i32> {
    let mut res: Vec<i32> = (0..x.len() as i32).collect();
    res.sort_unstable_by(|&i, &j| x[i as usize].partial_cmp(&x[j as usize]).unwrap());
    res
}

/// Inverts an order permutation: returns the rank of each element from its
/// `order` (sorted-position) array.
///
/// `vector<int> ranksFromOrder(vector<int> const& order)`
pub fn ranksFromOrder(order: &[i32]) -> Vec<i32> {
    let mut res = vec![0i32; order.len()];

    for i in 0..order.len() as i32 {
        res[order[i as usize] as usize] = i;
    }

    res
}

/// Computes the running GSEA statistic (maximum enrichment deviation) for all
/// prefixes of a gene set, using the square-root convex-hull heuristic.
///
/// # Arguments
/// * `stats` - gene-level statistics, sorted in decreasing order (not checked).
/// * `selectedStats` - 1-based indices of the selected genes into `stats`.
/// * `selectedOrder` - the order permutation of `selectedStats` (see [`order`]).
/// * `gseaParam` - GSEA weight parameter (0 is unweighted, suggested value 1).
/// * `rev` - when true, computes the statistic in the reversed (negative) direction.
///
/// # Returns
/// A vector of the GSEA statistic for each prefix of `selectedStats`.
///
/// `NumericVector gseaStats1(NumericVector const& stats, IntegerVector const& selectedStats, vector<int> const& selectedOrder, double gseaParam, bool rev = false)`
pub fn gseaStats1(
    stats: &[f64],
    selectedStats: &[i32],
    selectedOrder: &[i32],
    gseaParam: f64,
    rev: bool,
) -> Vec<f64> {
    let n = stats.len() as i32;
    let mut NR: f64 = 0.0;
    let k = selectedStats.len() as i32;

    let mut res = vec![0.0f64; k as usize];

    let mut xs = SegmentTree::<i32>::new(k + 1);
    let mut ys = SegmentTree::<f64>::new(k + 1);

    let mut selectedRanks = ranksFromOrder(selectedOrder);

    if !rev {
        let mut prev = -1;
        for i in 0..k {
            let j = selectedOrder[i as usize];
            let t = selectedStats[j as usize] - 1;
            debug_assert!(t - prev >= 1);
            xs.inc(i, t - prev);
            prev = t;
        }
        xs.inc(k, n - 1 - prev);
    } else {
        let mut prev = n;
        for i in 0..k {
            selectedRanks[i as usize] = k - 1 - selectedRanks[i as usize];
            let j = selectedOrder[(k - 1 - i) as usize];
            let t = selectedStats[j as usize] - 1;
            debug_assert!(prev - t >= 1);
            xs.inc(i, prev - t);
            prev = t;
        }
        xs.inc(k, prev - 0);
    }

    // previous node in stack
    let mut stPrev = vec![0i32; (k + 2) as usize];
    let mut stNext = vec![0i32; (k + 2) as usize];

    let k1 = ((k + 1) as f64).sqrt() as i32;
    let k2 = (k + 1) / k1 + 1;
    let mut blockSummit = vec![0i32; k2 as usize];
    let mut blockStart = vec![0i32; k2 as usize];
    let mut blockEnd = vec![0i32; k2 as usize];

    let mut i = 0;
    while i <= k + 1 {
        let block = i / k1;
        blockStart[block as usize] = i;
        blockEnd[block as usize] = std::cmp::min(i + k1 - 1, k + 1);
        let mut j = 1;
        while i + j <= blockEnd[block as usize] {
            stPrev[(i + j) as usize] = i + j - 1;
            stNext[(i + j - 1) as usize] = i + j;
            j += 1;
        }
        stPrev[i as usize] = i;
        stNext[blockEnd[block as usize] as usize] = blockEnd[block as usize];
        blockSummit[block as usize] = blockEnd[block as usize];
        i += k1;
    }

    let mut statEps: f64 = 1e-5;

    for i in 0..k {
        let t = selectedStats[i as usize] - 1;
        let xx = stats[t as usize].abs();
        if xx > 0.0 {
            statEps = xx.min(statEps);
        }
    }
    statEps /= 1024.0;

    for i in 0..k {
        let t = selectedStats[i as usize] - 1;
        let tRank = selectedRanks[i as usize];
        // 0 values make problems, replacing with epsilon
        let adjStat = stats[t as usize].abs().max(statEps).powf(gseaParam);

        xs.inc(tRank, -1);
        ys.inc(tRank, adjStat);
        NR += adjStat;

        let m = i + 1;

        let curBlock = (tRank + 1) / k1;
        let bS = blockStart[curBlock as usize];
        let bE = blockEnd[curBlock as usize];

        // Redoing upper convext hull for invalidated block `curBlock`

        let mut curTop = std::cmp::max(tRank, bS);

        let mut j = tRank + 1;
        while j <= bE {
            let c = j;
            let xc = xs.queryR(c) as f64;
            let yc = ys.queryR(c);

            let mut b = curTop;

            let mut xb = xs.queryR(b) as f64;
            let mut yb = ys.queryR(b);

            while stPrev[curTop as usize] != curTop {
                let a = stPrev[curTop as usize];

                let xa = xs.queryR(a) as f64;
                let ya = ys.queryR(a);

                let mut pr = (xb - xa) * (yc - yb) - (yb - ya) * (xc - xb);
                if yc - ya < eps {
                    pr = 0.0;
                }
                if pr <= 0.0 {
                    // right turn
                    break;
                }
                // left turn
                curTop = a;
                stNext[b as usize] = -1;
                b = a;
                xb = xa;
                yb = ya;
            }

            stPrev[c as usize] = curTop;
            stNext[curTop as usize] = c;
            curTop = c;
            if stNext[c as usize] != -1 {
                break;
            }
            j += 1;
        }

        let coef = (n - m) as f64 / NR;

        // Finding farthest points for `curBlock` from scratch

        let mut maxP: f64 = 0.0;

        // So that blockSummit[curBlock] is to the right of
        // actual summit as for every other block
        blockSummit[curBlock as usize] = std::cmp::max(curTop, blockSummit[curBlock as usize]);

        // Updating farthest points in valid blocks

        for block in 0..k2 {
            let mut curSummit = blockSummit[block as usize];

            let mut curDist = ys.queryR(curSummit) * coef - xs.queryR(curSummit) as f64;

            loop {
                let nextSummit = stPrev[curSummit as usize];
                let nextDist = ys.queryR(nextSummit) * coef - xs.queryR(nextSummit) as f64;

                if nextDist <= curDist {
                    break;
                }
                curDist = nextDist;
                curSummit = nextSummit;
            }

            blockSummit[block as usize] = curSummit;
            maxP = maxP.max(curDist);
        }

        maxP /= (n - m) as f64;

        res[i as usize] = maxP;
    }
    res
}

/// Gathers the elements of `from` at the given 1-based `indices`.
///
/// `NumericVector subvector(NumericVector const& from, IntegerVector const& indices)`
pub fn subvector(from: &[f64], indices: &[i32]) -> Vec<f64> {
    let mut result = vec![0.0f64; indices.len()];
    for i in 0..indices.len() {
        result[i] = from[(indices[i] - 1) as usize];
    }
    result
}

/// Calculates GSEA statistic values for all the prefixes of a gene set.
///
/// Takes `O(k^{3/2})` time, where `k` is the size of `selectedStats`.
///
/// # Arguments
/// * `stats` - gene-level statistics, sorted in decreasing order (order is not checked).
/// * `selectedStats` - indexes of selected genes in the `stats` array (1-based).
/// * `gseaParam` - GSEA weight parameter (0 is unweighted, suggested value is 1).
/// * `scoreType` - one of `"std"`, `"pos"`, `"neg"`.
///
/// # Returns
/// Numeric vector of GSEA statistics for all prefixes of `selectedStats`.
///
/// `NumericVector calcGseaStatCumulative(NumericVector const& stats, IntegerVector const& selectedStats, double gseaParam, std::string scoreType = "std")`
pub fn calcGseaStatCumulative(
    stats: &[f64],
    selectedStats: &[i32],
    gseaParam: f64,
    scoreType: &str,
) -> Vec<f64> {
    let selectedOrder = order(selectedStats);

    if scoreType != "std" && scoreType != "pos" && scoreType != "neg" {
        panic!("scoreType must take values from (\"std\", \"pos\", \"neg\")");
    }

    if scoreType == "std" {
        let mut res = gseaStats1(stats, selectedStats, &selectedOrder, gseaParam, false);
        let resDown = gseaStats1(stats, selectedStats, &selectedOrder, gseaParam, true);
        for i in 0..selectedStats.len() {
            if res[i] == resDown[i] {
                res[i] = 0.0;
            } else if res[i] < resDown[i] {
                res[i] = -resDown[i];
            }
        }
        res
    } else if scoreType == "pos" {
        let res = gseaStats1(stats, selectedStats, &selectedOrder, gseaParam, false);
        res
    } else {
        let mut res = gseaStats1(stats, selectedStats, &selectedOrder, gseaParam, true);
        for v in res.iter_mut() {
            *v = -*v;
        }
        res
    }
}

/// Draws a random size-`k` gene set from `1..=n` and computes its cumulative
/// GSEA statistics via [`calcGseaStatCumulative`].
///
/// `NumericVector calcRandomGseaStatCumulative(NumericVector const& stats, int n, int k, double gseaParam, random_engine_t& rng, std::string scoreType)`
pub fn calcRandomGseaStatCumulative(
    stats: &[f64],
    n: i32,
    k: i32,
    gseaParam: f64,
    rng: &mut random_engine_t,
    scoreType: &str,
) -> Vec<f64> {
    let selectedStats = combination(1, n, k, rng);

    calcGseaStatCumulative(stats, &selectedStats, gseaParam, scoreType)
}

/// Return value of `calcGseaStatCumulativeBatch` (the Rcpp `List` of six named vectors).
pub struct GseaStatCumulativeBatch {
    pub leEs: Vec<f64>,
    pub geEs: Vec<f64>,
    pub leZero: Vec<f64>,
    pub geZero: Vec<f64>,
    pub leZeroSum: Vec<f64>,
    pub geZeroSum: Vec<f64>,
}

/// Runs `iterations` random samples and, for each pathway size, accumulates the
/// counts and sums used to estimate enrichment-score p-values: how often a
/// random ES is `<=` / `>=` the observed `pathwayScores`, how often it is
/// `<=` / `>=` zero, and the running sums of the negative and positive parts.
///
/// `List calcGseaStatCumulativeBatch(NumericVector const& stats, double gseaParam, NumericVector const& pathwayScores, IntegerVector const& pathwaysSizes, int iterations, int seed, std::string scoreType)`
pub fn calcGseaStatCumulativeBatch(
    stats: &[f64],
    gseaParam: f64,
    pathwayScores: &[f64],
    pathwaysSizes: &[i32],
    iterations: i32,
    seed: i32,
    scoreType: &str,
) -> GseaStatCumulativeBatch {
    let n = stats.len() as i32;
    let k = *pathwaysSizes.iter().max().unwrap();
    let m = pathwaysSizes.len();

    let mut leEs = vec![0.0f64; m];
    let mut geEs = vec![0.0f64; m];
    let mut leZero = vec![0.0f64; m];
    let mut geZero = vec![0.0f64; m];
    let mut leZeroSum = vec![0.0f64; m];
    let mut geZeroSum = vec![0.0f64; m];

    let mut rng = random_engine_t::new(seed as u32);

    for _ in 0..iterations {
        let randEs = calcRandomGseaStatCumulative(stats, n, k, gseaParam, &mut rng, scoreType);
        let randEsP = subvector(&randEs, pathwaysSizes);

        for j in 0..m {
            leEs[j] += (randEsP[j] <= pathwayScores[j]) as i32 as f64;
            geEs[j] += (randEsP[j] >= pathwayScores[j]) as i32 as f64;
            leZero[j] += (randEsP[j] <= 0.0) as i32 as f64;
            geZero[j] += (randEsP[j] >= 0.0) as i32 as f64;
            leZeroSum[j] += randEsP[j].min(0.0);
            geZeroSum[j] += randEsP[j].max(0.0);
        }
    }

    GseaStatCumulativeBatch {
        leEs,
        geEs,
        leZero,
        geZero,
        leZeroSum,
        geZeroSum,
    }
}

/// Calculates GSEA statistic values for all gene sets in the `selectedGenes` list.
///
/// Takes `O(n + m*k*log(k))` time, where `n` is the number of genes, `m` is the
/// number of gene sets, and `k` is the mean gene set size.
///
/// # Arguments
/// * `stats` - gene-level statistics, sorted in decreasing order.
/// * `selectedGenes` - list of integer vectors of gene IDs (from 1 to n).
/// * `geneRanks` - gene ranks (1-based).
///
/// # Returns
/// Numeric vector of GSEA statistics, the same length as `selectedGenes`.
///
/// `NumericVector calcGseaStatBatchCpp(NumericVector const& stats, List const& selectedGenes, IntegerVector const& geneRanks)`
pub fn calcGseaStatBatchCpp(
    stats: &[f64],
    selectedGenes: &[Vec<i32>],
    geneRanks: &[i32],
) -> Vec<f64> {
    let n = stats.len() as i32;
    let m = selectedGenes.len();
    let mut gseaStats = vec![0.0f64; m];

    for i in 0..m {
        let mut S: Vec<i32> = selectedGenes[i].clone();
        for j in 0..S.len() {
            S[j] = geneRanks[(S[j] - 1) as usize];
        }
        S.sort();

        let k = S.len() as i32;
        gseaStats[i] = k as f64;

        let mut NR: f64 = 0.0;
        for j in 0..k {
            let tRank = S[j as usize]; // 1-based indexing
            let cur = stats[(tRank - 1) as usize].abs();
            NR += cur;
        }

        #[allow(unused_assignments)]
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;

        let coef = (n - k) as f64 / NR;

        let mut maxP: f64 = 0.0;
        let mut minP: f64 = 0.0;

        for j in 0..k {
            let tRank = S[j as usize]; // 1-based indexing
            let cur = stats[(tRank - 1) as usize].abs();

            x = (tRank - j - 1) as f64;

            let bottom = y * coef - x;

            y += cur;

            let top = y * coef - x;
            maxP = maxP.max(top);
            minP = minP.min(bottom);
        }

        if maxP > -minP {
            gseaStats[i] = maxP;
        } else if maxP < -minP {
            gseaStats[i] = minP;
        } else {
            gseaStats[i] = 0.0;
        }
        gseaStats[i] /= (n - k) as f64;
    }

    gseaStats
}
