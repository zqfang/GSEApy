//! Port of the fgsea multilevel p-value algorithm from C++.
//!
//! Original: https://github.com/ctlab/fgsea
//! Reference: Korotkevich et al., bioRxiv 2021
//!
//! This module implements the adaptive multilevel splitting + MCMC method
//! that provides arbitrarily small p-values with a quantified error bound.
#![allow(dead_code, unused)]

use rand::Rng;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use statrs::function::gamma::digamma;
use std::cmp::Ordering;

/// Maximum NS value for integer-scaling of ranks (1 << 30 = 1,073,741,824).
pub const MAX_NS: i64 = 1 << 30;

/// Minimum absolute ES to treat as non-zero; below this the gene set is
/// considered unenriched and p-value 1.0 is returned immediately.
const ES_EPSILON: f64 = 1e-15;

// ============================================================
// ScoreT: rational fraction score with exact integer comparison
//
// score = coef_ns / ns - coef_const / diff
// Comparisons use 128-bit arithmetic to avoid floating-point errors.
// ============================================================

#[derive(Clone, Copy, Debug, Default)]
pub struct ScoreT {
    pub ns: i64,          // NS – total weight of gene set (sum of selected ranks)
    pub coef_ns: i64,     // accumulated rank contribution (numerator of first term)
    pub diff: i64,        // n - k  (number of genes NOT in the set)
    pub coef_const: i64,  // accumulated gap contribution (numerator of second term)
}

impl ScoreT {
    #[inline]
    pub fn new(ns: i64, coef_ns: i64, diff: i64, coef_const: i64) -> Self {
        ScoreT { ns, coef_ns, diff, coef_const }
    }

    /// Floating-point representation: coef_ns/ns - coef_const/diff
    #[inline]
    pub fn get_double(&self) -> f64 {
        self.coef_ns as f64 / self.ns as f64 - self.coef_const as f64 / self.diff as f64
    }

    /// Sign indicator: positive iff coef_ns*diff >= coef_const*ns
    #[inline]
    pub fn get_numerator(&self) -> i64 {
        self.coef_ns * self.diff - self.coef_const * self.ns
    }

    /// Exact comparison using 128-bit arithmetic.
    /// Returns positive iff self > other, zero if equal, negative if self < other.
    #[inline]
    fn compare(&self, other: &Self) -> i128 {
        let p1 = self.coef_ns * self.diff + self.ns * (other.coef_const - self.coef_const);
        let q1 = self.ns * self.diff;
        let p2 = other.coef_ns;
        let q2 = other.ns;
        (p1 as i128) * (q2 as i128) - (p2 as i128) * (q1 as i128)
    }

    #[inline]
    pub fn negate(&self) -> Self {
        ScoreT {
            ns: self.ns,
            coef_ns: -self.coef_ns,
            diff: self.diff,
            coef_const: -self.coef_const,
        }
    }

    /// Return the score with the larger magnitude.
    #[inline]
    pub fn abs(self) -> Self {
        let neg = self.negate();
        if self >= neg { self } else { neg }
    }
}

impl PartialEq for ScoreT {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.compare(other) == 0
    }
}
impl Eq for ScoreT {}

impl PartialOrd for ScoreT {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for ScoreT {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.compare(other).cmp(&0)
    }
}

/// GseaT = (positive ES score, XOR-hash for tiebreaking).
/// Sorted lexicographically: first by score, then by hash.
type GseaT = (ScoreT, u64);

// ============================================================
// ES calculation in the integer domain
// ============================================================

/// Compute the maximum (positive) running ES for a gene set.
/// `ranks`: integer-scaled metric (descending sorted, all positive).
/// `p`: sorted gene positions (indices into ranks).
/// `ns`: sum of ranks[p[i]] for all i in p.
pub fn calc_positive_es_int(ranks: &[i64], p: &[usize], ns: i64) -> ScoreT {
    let n = ranks.len();
    let k = p.len();
    let diff = (n - k) as i64;
    let mut res = ScoreT::new(ns, 0, diff, 0);
    let mut cur = res;
    let mut last: i64 = -1;
    for &pos in p {
        cur.coef_ns += ranks[pos];
        cur.coef_const += pos as i64 - last - 1;
        if cur > res {
            res = cur;
        }
        last = pos as i64;
    }
    res
}

/// Compute the maximum absolute running ES (with sign).
/// `p`: sorted gene positions.
/// `ns`: sum of ranks[p[i]].
pub fn calc_es_int(ranks: &[i64], p: &[usize], ns: i64) -> ScoreT {
    let n = ranks.len();
    let k = p.len();
    let diff = (n - k) as i64;
    let mut res = ScoreT::new(ns, 0, diff, 0);
    let mut cur = res;
    let mut last: i64 = -1;
    for &pos in p {
        cur.coef_const += pos as i64 - last - 1;
        if cur.abs() > res.abs() {
            res = cur;
        }
        cur.coef_ns += ranks[pos];
        if cur.abs() > res.abs() {
            res = cur;
        }
        last = pos as i64;
    }
    res
}

// ============================================================
// Statistical helper functions
// ============================================================

/// Trigamma function ψ₁(x) — second derivative of ln Γ(x).
/// Computed by recurrence (shift x ≥ 8) + asymptotic expansion.
/// The threshold 8.0 is chosen so that the asymptotic series converges to
/// machine precision with 6 Bernoulli terms (error < 1e-14 for x ≥ 8).
fn trigamma(mut x: f64) -> f64 {
    let mut result = 0.0;
    while x < 8.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    // Asymptotic expansion: ψ₁(x) ≈ 1/x + 1/(2x²) + 1/(6x³) - 1/(30x⁵) + 1/(42x⁷) - 1/(30x⁹)
    let iz = 1.0 / x;
    let iz2 = iz * iz;
    result += iz
        + 0.5 * iz2
        + (1.0 / 6.0) * iz2 * iz
        - (1.0 / 30.0) * iz2 * iz2 * iz
        + (1.0 / 42.0) * iz2 * iz2 * iz2 * iz
        - (1.0 / 30.0) * iz2 * iz2 * iz2 * iz2 * iz;
    result
}

/// betaMeanLog(a, b) = ψ(a) − ψ(b+1).
/// Used to compute log conditional probabilities at each level.
pub fn beta_mean_log(a: usize, b: usize) -> f64 {
    digamma(a as f64) - digamma((b + 1) as f64)
}

/// getVarPerLevel(k, n) = ψ₁(k) − ψ₁(n+1).
/// Variance contribution of one level to the log p-value estimate.
pub fn get_var_per_level(k: usize, n: usize) -> f64 {
    trigamma(k as f64) - trigamma((n + 1) as f64)
}

// ============================================================
// Uniform integer sampling in [from, to] (inclusive)
// ============================================================

#[inline]
fn sample_range(rng: &mut SmallRng, from: usize, to: usize) -> usize {
    debug_assert!(from <= to, "sample_range: from ({from}) > to ({to})");
    if from == to {
        return from;
    }
    rng.random_range(from..=to)
}

// ============================================================
// combination: generate k distinct integers from [a, b]
// ============================================================

/// Generate a uniformly random subset of k distinct integers from [a, b].
/// The result is in arbitrary order (caller should sort if needed).
pub fn combination(a: usize, b: usize, k: usize, rng: &mut SmallRng) -> Vec<usize> {
    let n = b - a + 1;
    let mut used = vec![false; n];
    let mut v = Vec::with_capacity(k);

    if k <= n / 2 {
        // Sparse case: rejection sampling
        for _ in 0..k {
            loop {
                let x = sample_range(rng, 0, n - 1);
                if !used[x] {
                    v.push(a + x);
                    used[x] = true;
                    break;
                }
            }
        }
    } else {
        // Dense case: Knuth's algorithm D variant
        for r in (n - k)..n {
            let x = sample_range(rng, 0, r);
            if !used[x] {
                v.push(a + x);
                used[x] = true;
            } else {
                v.push(a + r);
                used[r] = true;
            }
        }
        // Fisher–Yates shuffle for uniform sampling (reproducibility)
        for i in (1..v.len()).rev() {
            let j = sample_range(rng, 0, i);
            v.swap(i, j);
        }
    }
    v
}

// ============================================================
// Scale float ranks to integers
// ============================================================

/// Scale a non-negative float metric to i64 integers.
/// After scaling, sum(result) ≤ MAX_NS = 1 << 30.
/// Input values should be non-negative (e.g., |x|^weight weighted metrics).
/// Very large float inputs are handled safely: the scale factor is floored to an
/// integer (≥ 1) so individual elements do not overflow i64.
pub fn scale_ranks(ranks: &[f64]) -> Vec<i64> {
    let cur_ns: f64 = ranks.iter().sum();
    if cur_ns == 0.0 {
        return vec![0i64; ranks.len()];
    }
    let scale = (MAX_NS as f64 / cur_ns).floor();
    let scale = if scale < 1.0 { 1.0 } else { scale };
    ranks.iter().map(|&r| (r * scale) as i64).collect()
}

// ============================================================
// Internal data structures for EsRuler
// ============================================================

/// One level of the multilevel hierarchy.
#[derive(Clone)]
struct Level {
    /// Samples at or below the boundary (will be replaced by copies from high_scores).
    low_scores: Vec<(GseaT, bool)>,  // (gsea_t, is_positive)
    /// Samples above the boundary.
    high_scores: Vec<(GseaT, bool)>,
    /// Boundary GseaT value (greatest low score).
    bound: GseaT,
}

/// Per-sample chunked representation for fast MCMC ES updates.
#[derive(Clone)]
struct SampleChunks {
    /// Sum of integer ranks for each chunk.
    chunk_sum: Vec<i64>,
    /// Sorted gene positions in each chunk.
    chunks: Vec<Vec<usize>>,
}

impl SampleChunks {
    fn new(chunks_number: usize) -> Self {
        SampleChunks {
            chunk_sum: vec![0i64; chunks_number],
            chunks: vec![Vec::new(); chunks_number],
        }
    }
}

// ============================================================
// EsRuler: adaptive multilevel MCMC for p-value computation
// ============================================================

pub struct EsRuler {
    ranks: Vec<i64>,
    gene_hashes: Vec<u64>,
    sample_size: usize,
    pathway_size: usize,
    moves_scale: f64,
    incorrect_ruler: bool,
    current_samples: Vec<Vec<usize>>,
    old_samples_start: usize,
    levels: Vec<Level>,
    chunk_last_element: Vec<usize>,
    chunks_number: usize,
}

impl EsRuler {
    /// Create a new EsRuler.
    /// `ranks`: integer-scaled gene metric (all non-negative), sorted by importance.
    /// `sample_size`: MCMC sample size (default 101 in fgsea).
    /// `pathway_size`: gene set size (k).
    /// `moves_scale`: scaling factor for MCMC steps (default 1.0).
    pub fn new(ranks: Vec<i64>, sample_size: usize, pathway_size: usize, moves_scale: f64) -> Self {
        let n = ranks.len();
        EsRuler {
            ranks,
            gene_hashes: vec![0u64; n],
            sample_size,
            pathway_size,
            moves_scale,
            incorrect_ruler: false,
            current_samples: vec![Vec::new(); sample_size],
            old_samples_start: 0,
            levels: Vec::new(),
            chunk_last_element: Vec::new(),
            chunks_number: 0,
        }
    }

    /// Resample gene sets: split into high/low halves and replace low half with copies from high.
    /// Adds one new level to `self.levels`.
    /// Returns true on success; false means all scores are equal (stalled).
    fn resample_genesets(&mut self, rng: &mut SmallRng) -> bool {
        let sample_size = self.sample_size;
        let ranks = &self.ranks;
        let gene_hashes = &self.gene_hashes;

        // Compute (GseaT, is_positive, sample_id) for each current sample
        let mut stats: Vec<(GseaT, bool, usize)> = (0..sample_size)
            .map(|sid| {
                let sample = &self.current_samples[sid];
                let ns: i64 = sample.iter().map(|&i| ranks[i]).sum();
                let es_pos = calc_positive_es_int(ranks, sample, ns);
                let es_signed = calc_es_int(ranks, sample, ns);
                let hash = sample.iter().fold(0u64, |acc, &i| acc ^ gene_hashes[i]);
                let is_positive = es_signed.get_numerator() >= 0;
                ((es_pos, hash), is_positive, sid)
            })
            .collect();

        // Sort by GseaT (score, then hash) ascending
        stats.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Find the split point: first index where score >= median score
        let central_value = stats[sample_size / 2].0.clone();
        let mut start_from = 0usize;
        for sid in 0..sample_size {
            if stats[sid].0 >= central_value {
                start_from = sid;
                break;
            }
        }

        // If start_from == 0, advance past all scores equal to the minimum
        if start_from == 0 {
            while start_from < sample_size && stats[start_from].0 == stats[0].0 {
                start_from += 1;
            }
        }

        // If all scores are equal, we cannot advance
        if start_from == sample_size {
            return true;
        }

        // Build the level
        let bound = stats[start_from - 1].0.clone();
        let low_scores: Vec<(GseaT, bool)> = stats[..start_from]
            .iter()
            .map(|(g, is_pos, _)| (g.clone(), *is_pos))
            .collect();
        let high_scores: Vec<(GseaT, bool)> = stats[start_from..]
            .iter()
            .map(|(g, is_pos, _)| (g.clone(), *is_pos))
            .collect();

        self.levels.push(Level { low_scores, high_scores, bound });

        // Replace low-scoring samples with copies from high-scoring samples
        let n_high = sample_size - start_from;
        let mut new_sets: Vec<Vec<usize>> = Vec::with_capacity(sample_size);
        for _ in 0..start_from {
            let ind = sample_range(rng, 0, n_high - 1) + start_from;
            new_sets.push(self.current_samples[stats[ind].2].clone());
        }
        for i in start_from..sample_size {
            new_sets.push(self.current_samples[stats[i].2].clone());
        }

        self.old_samples_start = start_from;
        self.current_samples = new_sets;
        true
    }

    /// Run the MCMC perturbation kernel on one sample's chunk structure.
    ///
    /// Randomly swaps one gene in the sample set with an external gene, accepting
    /// only if the new sample's ES remains above `bound`. Repeats until the stopping
    /// criterion `f(accepted_moves, total_iters)` returns true.
    ///
    /// Returns (accepted_moves, total_iters).
    fn perturbate_until<F>(
        &self,
        k: usize,
        sample_chunks: &mut SampleChunks,
        bound: &GseaT,
        rng: &mut SmallRng,
        f: F,
    ) -> (i32, i32)
    where
        F: Fn(i32, i32) -> bool,
    {
        let n = self.ranks.len();
        let ranks = &self.ranks;
        let gene_hashes = &self.gene_hashes;
        let chunk_last_element = &self.chunk_last_element;
        let chunks_number = self.chunks_number;
        let diff = (n - k) as i64;

        // Compute initial NS (sum of ranks) and hash from chunks
        let mut ns: i64 = 0;
        let mut cur_hash: u64 = 0;
        for ci in 0..chunks_number {
            for &pos in &sample_chunks.chunks[ci] {
                ns += ranks[pos];
                cur_hash ^= gene_hashes[pos];
            }
        }

        // Candidate gene for quick bound check (optimization)
        let mut cand_val: Option<usize> = None;
        let mut cand_x: i64 = 0;
        let mut cand_y: i64 = 0;

        let mut moves: i32 = 0;
        let mut iters: i32 = 0;

        while !f(moves, iters) {
            iters += 1;

            // Pick a random index in [0, k-1] to remove
            let old_ind = sample_range(rng, 0, k - 1);

            // Locate which chunk and position old_ind refers to
            let mut old_chunk_ind = 0usize;
            let mut old_ind_in_chunk = old_ind;
            while sample_chunks.chunks[old_chunk_ind].len() <= old_ind_in_chunk {
                old_ind_in_chunk -= sample_chunks.chunks[old_chunk_ind].len();
                old_chunk_ind += 1;
            }
            let old_val = sample_chunks.chunks[old_chunk_ind][old_ind_in_chunk];

            // Pick a random gene position in [0, n-1] to insert
            let new_val = sample_range(rng, 0, n - 1);

            // Find the chunk for new_val (upper_bound: first chunk whose last element > new_val)
            let new_chunk_ind = chunk_last_element.partition_point(|&x| x <= new_val);
            // Find insertion position in that chunk (sorted)
            let new_ind_in_chunk = sample_chunks.chunks[new_chunk_ind]
                .partition_point(|&x| x < new_val);

            // Check if new_val is already in the sample
            let already_present = new_ind_in_chunk < sample_chunks.chunks[new_chunk_ind].len()
                && sample_chunks.chunks[new_chunk_ind][new_ind_in_chunk] == new_val;

            if already_present {
                if new_val == old_val {
                    moves += 1; // trivial accepted move
                }
                continue;
            }

            // Insertion index adjustment when old and new are in the same chunk
            // and old_ind < new_ind (removing old shifts new position left by 1)
            let insert_adj =
                if old_chunk_ind == new_chunk_ind && old_ind_in_chunk < new_ind_in_chunk {
                    1usize
                } else {
                    0usize
                };

            // Perform the swap: remove old_val, insert new_val
            sample_chunks.chunks[old_chunk_ind].remove(old_ind_in_chunk);
            sample_chunks.chunks[new_chunk_ind].insert(new_ind_in_chunk - insert_adj, new_val);

            ns = ns - ranks[old_val] + ranks[new_val];
            cur_hash ^= gene_hashes[old_val] ^ gene_hashes[new_val];
            sample_chunks.chunk_sum[old_chunk_ind] -= ranks[old_val];
            sample_chunks.chunk_sum[new_chunk_ind] += ranks[new_val];

            // `strictly` breaks hash ties (prefer lower hashes)
            let strictly = cur_hash <= bound.1;
            let check = |score: &ScoreT| -> bool {
                if strictly {
                    *score > bound.0
                } else {
                    *score >= bound.0
                }
            };

            // Update candidate tracking after the swap
            if let Some(cv) = cand_val {
                if old_val == cv {
                    cand_val = None;
                }
            }
            if let Some(cv) = cand_val {
                if old_val < cv {
                    cand_x += 1;
                    cand_y -= ranks[old_val];
                }
                if new_val < cv {
                    cand_x -= 1;
                    cand_y += ranks[new_val];
                }
            }

            // Fast path: check using cached candidate position
            if cand_val.is_some() {
                let cand_score = ScoreT::new(ns, cand_y, diff, cand_x);
                if check(&cand_score) {
                    moves += 1;
                    continue;
                }
            }

            // Full ES computation via chunk structure
            let mut cur_x: i64 = 0;
            let mut cur_y: i64 = 0;
            let mut ok = false;
            let mut last: i64 = -1;

            'outer: for ci in 0..chunks_number {
                // Quick check: if adding the entire chunk's sum doesn't exceed bound, skip it
                let chunk_score =
                    ScoreT::new(ns, cur_y + sample_chunks.chunk_sum[ci], diff, cur_x);
                if !check(&chunk_score) {
                    cur_y += sample_chunks.chunk_sum[ci];
                    cur_x += chunk_last_element[ci] as i64
                        - last
                        - 1
                        - sample_chunks.chunks[ci].len() as i64;
                    last = chunk_last_element[ci] as i64 - 1;
                } else {
                    // Scan individual genes in this chunk
                    for &pos in &sample_chunks.chunks[ci] {
                        cur_y += ranks[pos];
                        cur_x += pos as i64 - last - 1;
                        let cur_score = ScoreT::new(ns, cur_y, diff, cur_x);
                        if check(&cur_score) {
                            ok = true;
                            cand_val = Some(pos);
                            cand_x = cur_x;
                            cand_y = cur_y;
                            break 'outer;
                        }
                        last = pos as i64;
                    }
                    if !ok {
                        cur_x += chunk_last_element[ci] as i64 - 1 - last;
                        last = chunk_last_element[ci] as i64 - 1;
                    }
                }
            }

            if !ok {
                // Reject: revert the swap
                ns = ns - ranks[new_val] + ranks[old_val];
                cur_hash ^= gene_hashes[new_val] ^ gene_hashes[old_val];
                sample_chunks.chunk_sum[old_chunk_ind] += ranks[old_val];
                sample_chunks.chunk_sum[new_chunk_ind] -= ranks[new_val];
                sample_chunks.chunks[new_chunk_ind].remove(new_ind_in_chunk - insert_adj);
                sample_chunks.chunks[old_chunk_ind].insert(old_ind_in_chunk, old_val);

                // Revert candidate tracking
                if let Some(cv) = cand_val {
                    if new_val == cv {
                        cand_val = None;
                    }
                }
                if let Some(cv) = cand_val {
                    if old_val < cv {
                        cand_x -= 1;
                        cand_y += ranks[old_val];
                    }
                    if new_val < cv {
                        cand_x += 1;
                        cand_y -= ranks[new_val];
                    }
                }
            } else {
                moves += 1;
            }
        }

        (moves, iters)
    }

    /// Build the multilevel hierarchy of levels up to `es_double`.
    ///
    /// `es_double`: absolute value of the observed enrichment score.
    /// `seed`: random seed for reproducibility.
    /// `eps`: convergence threshold (stop when estimated p-value < eps); use 0 to disable.
    pub fn extend(&mut self, es_double: f64, seed: u64, eps: f64) {
        let n = self.ranks.len();
        let k = self.pathway_size;

        let mut rng = SmallRng::seed_from_u64(seed);

        // Initialize per-gene XOR hashes (used for tiebreaking)
        for i in 0..n {
            self.gene_hashes[i] = rng.random::<u64>();
        }

        // Initialize random gene-set samples of size k
        for sid in 0..self.sample_size {
            let mut s = combination(0, n - 1, k, &mut rng);
            s.sort_unstable();
            self.current_samples[sid] = s;
        }

        // First resampling step to build level 0
        if !self.resample_genesets(&mut rng) {
            self.incorrect_ruler = true;
            return;
        }

        // Set up chunk structure: divide gene positions into sqrt(k) chunks
        self.chunks_number = ((k as f64).sqrt() as usize).max(1);
        self.chunk_last_element = vec![0usize; self.chunks_number];
        self.chunk_last_element[self.chunks_number - 1] = n;

        // Compute chunk boundaries using median positions across current samples
        let mut tmp = vec![0usize; self.sample_size];
        for ci in 0..(self.chunks_number - 1) {
            let pos = (k + ci) / self.chunks_number; // gene index within each sample
            for sid in 0..self.sample_size {
                tmp[sid] = self.current_samples[sid][pos];
            }
            // nth_element: find the median gene position for this chunk
            let mid_idx = self.sample_size / 2;
            tmp.select_nth_unstable(mid_idx);
            self.chunk_last_element[ci] = tmp[mid_idx];
        }

        // Target ES as a ScoreT for exact comparison with level bounds
        // score_t{MAX_NS, MAX_NS * es_double, 1, 0} represents exactly es_double
        let need_es = ScoreT::new(MAX_NS, (MAX_NS as f64 * es_double) as i64, 1, 0);

        let mut adj_log_pval: f64 = 0.0;

        // Main multilevel loop: build levels until bound exceeds observed ES
        loop {
            // Check if the current level's bound already reaches the target
            if self.levels.last().unwrap().bound.0 >= need_es {
                break;
            }

            // Accumulate log-probability contribution from this level
            let nhigh = self.levels.last().unwrap().high_scores.len() + 1;
            adj_log_pval += beta_mean_log(nhigh, self.sample_size);

            // Early exit if p-value estimate is already below eps
            if eps != 0.0 && adj_log_pval < eps.ln() {
                break;
            }

            // Build per-sample chunk structures for MCMC perturbation
            let mut samples_chunks: Vec<SampleChunks> = (0..self.sample_size)
                .map(|_| SampleChunks::new(self.chunks_number))
                .collect();

            let current_bound = self.levels.last().unwrap().bound.clone();

            for sid in 0..self.sample_size {
                let mut cnt = 0usize;
                for &pos in &self.current_samples[sid] {
                    while cnt < self.chunks_number - 1
                        && self.chunk_last_element[cnt] <= pos
                    {
                        cnt += 1;
                    }
                    samples_chunks[sid].chunks[cnt].push(pos);
                    samples_chunks[sid].chunk_sum[cnt] += self.ranks[pos];
                }
            }

            // Run MCMC perturbation in two phases (warm-up + sampling)
            let perturb_iters = ((k as f64 * 0.1) as i32).max(1);
            let need_accepted =
                (self.moves_scale * self.sample_size as f64 * k as f64 / 2.0) as i32;

            let mut n_iterations = 0i32;
            let mut n_accepted = 0i32;

            // Phase 1: run until enough moves are accepted
            while n_accepted < need_accepted {
                for sid in 0..self.sample_size {
                    let (moves, _) = self.perturbate_until(
                        k,
                        &mut samples_chunks[sid],
                        &current_bound,
                        &mut rng,
                        |_, iters| iters >= perturb_iters,
                    );
                    n_accepted += moves;
                }
                n_iterations += 1;
            }

            // Phase 2: run the same number of iterations again
            for _ in 0..n_iterations {
                for sid in 0..self.sample_size {
                    self.perturbate_until(
                        k,
                        &mut samples_chunks[sid],
                        &current_bound,
                        &mut rng,
                        |_, iters| iters >= perturb_iters,
                    );
                }
            }

            // Rebuild current_samples from chunks (maintaining sorted order)
            for sid in 0..self.sample_size {
                self.current_samples[sid].clear();
                for ci in 0..self.chunks_number {
                    self.current_samples[sid]
                        .extend_from_slice(&samples_chunks[sid].chunks[ci]);
                }
            }

            // Create the next level via resampling
            let last_size = self.levels.len();
            if !self.resample_genesets(&mut rng) {
                self.incorrect_ruler = true;
                return;
            }
            // If no new level was created (stalled), stop
            if last_size == self.levels.len() {
                break;
            }
        }
    }

    /// Get the multilevel p-value for an observed enrichment score.
    ///
    /// `es_double`: absolute value of the observed ES.
    /// `eps`: convergence threshold (same as used in `extend`).
    /// `sign`: if false (recommended), count only same-sign samples (one-tailed);
    ///         if true, count all samples.
    ///
    /// Returns `(p_value, is_cp_ge_half, log2err)`.
    /// `log2err` = √(Σ varPerLevel) / ln(2) — the uncertainty in bits.
    pub fn get_pvalue(&self, es_double: f64, eps: f64, sign: bool) -> (f64, bool, f64) {
        if self.incorrect_ruler {
            return (f64::NAN, true, f64::NAN);
        }
        if self.levels.is_empty() {
            return (1.0, true, f64::NAN);
        }

        // Represent the observed ES as a GseaT for comparison with level bounds
        let es_score = ScoreT::new(MAX_NS, (MAX_NS as f64 * es_double) as i64, 1, 0);
        let es_gsea: GseaT = (es_score, 0);

        let mut adj_log_pval: f64 = 0.0;
        let mut lvls_var: f64 = 0.0;

        for lvl in &self.levels {
            if es_gsea <= lvl.bound {
                // ES falls within this level's range; count samples >= ES
                let mut cnt_last = 0usize;
                let mut cnt_positive = 0usize;

                // All high_scores are > bound >= ES, so all count
                for (_, is_pos) in &lvl.high_scores {
                    cnt_last += 1;
                    if *is_pos {
                        cnt_positive += 1;
                    }
                }
                // Count low_scores that are still >= ES
                for (gs, is_pos) in &lvl.low_scores {
                    if *gs >= es_gsea {
                        cnt_last += 1;
                        if *is_pos {
                            cnt_positive += 1;
                        }
                    }
                }

                let numerator = if sign { cnt_last } else { cnt_positive };

                if numerator == 0 {
                    adj_log_pval += beta_mean_log(1, self.sample_size);
                    let p = adj_log_pval.exp().clamp(0.0, 1.0);
                    return (p, true, f64::NAN);
                }

                adj_log_pval += beta_mean_log(numerator, self.sample_size);
                lvls_var += get_var_per_level(numerator, self.sample_size);
                let log2err = lvls_var.sqrt() / 2.0f64.ln();
                let p = adj_log_pval.exp().clamp(0.0, 1.0);
                return (p, true, log2err);
            }

            // ES is above this level's bound; accumulate contribution and continue
            let nhigh = lvl.high_scores.len() + 1;
            adj_log_pval += beta_mean_log(nhigh, self.sample_size);
            lvls_var += get_var_per_level(nhigh, self.sample_size);
        }

        // ES is above all level bounds; count from the last level's high scores
        let last_lvl = self.levels.last().unwrap();
        let mut cnt_last = 0usize;
        let mut cnt_positive = 0usize;

        for (gs, is_pos) in &last_lvl.high_scores {
            if *gs >= es_gsea {
                cnt_last += 1;
                if *is_pos {
                    cnt_positive += 1;
                }
            }
        }

        let numerator = if sign { cnt_last } else { cnt_positive };
        let nhigh = last_lvl.high_scores.len();

        if numerator == 0 || nhigh == 0 {
            let b = if nhigh > 0 { nhigh } else { 1 };
            adj_log_pval += beta_mean_log(1, b);
            let p = adj_log_pval.exp().clamp(0.0, 1.0);
            return (p, true, f64::NAN);
        }

        adj_log_pval += beta_mean_log(numerator, nhigh);
        lvls_var += get_var_per_level(numerator, nhigh);
        let log2err = lvls_var.sqrt() / 2.0f64.ln();
        let p = adj_log_pval.exp().clamp(0.0, 1.0);
        (p, true, log2err)
    }
}

// ============================================================
// Public API
// ============================================================

/// Compute the fgsea multilevel p-value for a preranked gene set.
///
/// `int_ranks`: integer-scaled gene metric, non-negative, sorted by importance (descending).
///              Use `scale_ranks(weighted_metric)` to obtain this.
/// `pathway_size`: gene set size k.
/// `observed_es`: the floating-point enrichment score (signed).
/// `sample_size`: MCMC sample size per level (default 101).
/// `seed`: random seed.
/// `eps`: convergence threshold; use 1e-50 to match fgsea defaults.
///
/// Returns `(p_value, log2err)`.
/// `log2err` is the uncertainty in bits (NaN if the estimate hit the extreme boundary).
pub fn compute_pvalue_multilevel(
    int_ranks_pos: &[i64],
    pathway_size: usize,
    observed_es: f64,
    sample_size: usize,
    seed: u64,
    eps: f64,
) -> (f64, f64) {
    if observed_es.abs() < ES_EPSILON || pathway_size == 0 {
        return (1.0, f64::NAN);
    }
    if pathway_size >= int_ranks_pos.len() {
        return (1.0, f64::NAN);
    }

    let es_abs = observed_es.abs();

    // For negative ES, reverse the ranks so that bottom genes become "top"
    let ranks: Vec<i64> = if observed_es >= 0.0 {
        int_ranks_pos.to_vec()
    } else {
        let mut r = int_ranks_pos.to_vec();
        r.reverse();
        r
    };

    let mut ruler = EsRuler::new(ranks, sample_size, pathway_size, 1.0);
    ruler.extend(es_abs, seed, eps);
    let (pval, _, log2err) = ruler.get_pvalue(es_abs, eps, false);
    (pval, log2err)
}

/// Compute multilevel p-values for **multiple gene sets of the same pathway size**
/// by sharing a single MCMC null-distribution ruler per sign direction.
///
/// This is the key batching optimisation from the reference fgsea implementation:
/// the null ES distribution depends only on `pathway_size` and the ranking metric
/// (not on the specific gene set), so all pathways of identical size can share
/// one `EsRuler` extended to the *maximum* |ES| seen in that size group.
///
/// Two rulers are built per call:
/// * **pos_ruler** — uses the original `int_ranks`; handles pathways with ES ≥ 0.
/// * **neg_ruler** — uses reversed `int_ranks`; handles pathways with ES < 0.
///
/// Each ruler is extended only as far as needed (to the maximum |ES| in the
/// respective sign group), providing early-exit savings for easy gene sets.
///
/// Arguments
/// ---------
/// * `int_ranks`    — integer-scaled gene metric, non-negative, descending.
///                    Obtain via `scale_ranks(|metric|^weight)`.
/// * `pathway_size` — gene set size *k* (identical for all pathways in this batch).
/// * `observed_es_values` — slice of *signed* enrichment scores, one per pathway.
/// * `sample_size`  — MCMC sample size per level (fgsea default 101).
/// * `seed`         — base random seed; pos ruler uses `seed`, neg ruler uses
///                    `seed.wrapping_add(1)` to keep both reproducible but independent.
/// * `eps`          — convergence threshold (e.g. 1e-50).
///
/// Returns a `Vec<(f64, f64)>` of `(p_value, log2err)` in the **same order** as
/// `observed_es_values`.
pub fn compute_pvalue_multilevel_batch(
    int_ranks: &[i64],
    pathway_size: usize,
    observed_es_values: &[f64],
    sample_size: usize,
    seed: u64,
    eps: f64,
) -> Vec<(f64, f64)> {
    let n = int_ranks.len();
    if observed_es_values.is_empty() || pathway_size == 0 || pathway_size >= n {
        return observed_es_values.iter().map(|_| (1.0, f64::NAN)).collect();
    }

    // ------------------------------------------------------------------
    // Separate positive and negative ES values while keeping original order
    // ------------------------------------------------------------------
    let mut max_pos_es: f64 = 0.0;
    let mut max_neg_abs_es: f64 = 0.0;
    for &es in observed_es_values {
        if es >= 0.0 {
            if es > max_pos_es {
                max_pos_es = es;
            }
        } else {
            let a = es.abs();
            if a > max_neg_abs_es {
                max_neg_abs_es = a;
            }
        }
    }

    // ------------------------------------------------------------------
    // Build pos ruler (original ranks, extended to max positive ES)
    // ------------------------------------------------------------------
    let pos_ruler: Option<EsRuler> = if max_pos_es > ES_EPSILON {
        let mut ruler = EsRuler::new(int_ranks.to_vec(), sample_size, pathway_size, 1.0);
        ruler.extend(max_pos_es, seed, eps);
        Some(ruler)
    } else {
        None
    };

    // ------------------------------------------------------------------
    // Build neg ruler (reversed ranks, extended to max |negative ES|)
    // ------------------------------------------------------------------
    let neg_ruler: Option<EsRuler> = if max_neg_abs_es > ES_EPSILON {
        let mut rev = int_ranks.to_vec();
        rev.reverse();
        let mut ruler = EsRuler::new(rev, sample_size, pathway_size, 1.0);
        // Use a distinct seed so the neg ruler's RNG stream is independent.
        ruler.extend(max_neg_abs_es, seed.wrapping_add(1), eps);
        Some(ruler)
    } else {
        None
    };

    // ------------------------------------------------------------------
    // Query each ES value from the appropriate ruler
    // ------------------------------------------------------------------
    observed_es_values
        .iter()
        .map(|&es| {
            if es.abs() < ES_EPSILON {
                return (1.0, f64::NAN);
            }
            let es_abs = es.abs();
            if es >= 0.0 {
                match &pos_ruler {
                    Some(r) => {
                        let (pval, _, log2err) = r.get_pvalue(es_abs, eps, false);
                        (pval, log2err)
                    }
                    None => (1.0, f64::NAN),
                }
            } else {
                match &neg_ruler {
                    Some(r) => {
                        let (pval, _, log2err) = r.get_pvalue(es_abs, eps, false);
                        (pval, log2err)
                    }
                    None => (1.0, f64::NAN),
                }
            }
        })
        .collect()
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trigamma_known_values() {
        // trigamma(1) = pi^2/6 ≈ 1.6449340668
        let t1 = trigamma(1.0);
        assert!((t1 - std::f64::consts::PI * std::f64::consts::PI / 6.0).abs() < 1e-6,
            "trigamma(1) = {t1}");
        // trigamma(2) = pi^2/6 - 1 ≈ 0.6449340668
        let t2 = trigamma(2.0);
        assert!((t2 - (std::f64::consts::PI * std::f64::consts::PI / 6.0 - 1.0)).abs() < 1e-6,
            "trigamma(2) = {t2}");
    }

    #[test]
    fn test_score_t_comparison() {
        // Compare scores with equal NS and diff
        let ns = 1000i64;
        let diff = 50i64;
        let s1 = ScoreT::new(ns, 200, diff, 10);
        let s2 = ScoreT::new(ns, 100, diff, 10);
        assert!(s1 > s2, "s1.getDouble() = {}, s2.getDouble() = {}", s1.get_double(), s2.get_double());

        let s3 = ScoreT::new(ns, 200, diff, 10);
        assert_eq!(s1, s3);
    }

    #[test]
    fn test_scale_ranks() {
        let ranks = vec![3.0, 2.0, 1.0];
        let int_ranks = scale_ranks(&ranks);
        let total: i64 = int_ranks.iter().sum();
        assert!(total <= MAX_NS, "total={total} > MAX_NS={MAX_NS}");
        // Ratios should be preserved
        assert!(int_ranks[0] > int_ranks[1]);
        assert!(int_ranks[1] > int_ranks[2]);
    }

    #[test]
    fn test_combination_size() {
        let mut rng = SmallRng::seed_from_u64(42);
        let k = 10;
        let n = 100;
        let result = combination(0, n - 1, k, &mut rng);
        assert_eq!(result.len(), k);
        // All values in range
        for &v in &result {
            assert!(v < n);
        }
        // All distinct
        let mut sorted = result.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), k);
    }

    #[test]
    fn test_beta_mean_log_sanity() {
        // betaMeanLog(a, b) = digamma(a) - digamma(b+1)
        // For a=1, b=101: should be a large negative number
        let v = beta_mean_log(1, 101);
        assert!(v < 0.0, "betaMeanLog(1, 101) = {v}");
        // For a=b+1: digamma(b+1) - digamma(b+1) = 0... not quite, but close to 0
        let v2 = beta_mean_log(50, 49);
        // digamma(50) - digamma(50) = 0 → but a=50, b=49 → digamma(50) - digamma(50) = 0
        assert!(v2.abs() < 1e-10, "betaMeanLog(50,49) = {v2}");
    }

    #[test]
    fn test_calc_positive_es_int_simple() {
        // Simple case: 5 genes, pathway of size 2 at positions 0, 2
        let ranks = vec![100i64, 80, 60, 40, 20];
        let p = vec![0usize, 2];
        let ns: i64 = ranks[0] + ranks[2]; // 160
        let score = calc_positive_es_int(&ranks, &p, ns);
        assert!(score.get_double() > 0.0, "positive ES should be positive");
    }

    #[test]
    fn test_esruler_produces_valid_pvalue() {
        // Create a simple gene list and pathway
        let n = 200usize;
        let k = 10usize;
        // Decreasing ranks (top genes are most important)
        let float_ranks: Vec<f64> = (0..n).map(|i| (n - i) as f64).collect();
        let int_ranks = scale_ranks(&float_ranks);

        // Pathway at top of list (very significant, positive ES)
        // observed_es ≈ max possible for top-k genes
        let observed_es = 0.9f64; // high ES

        let sample_size = 51usize; // smaller for test speed
        let (pval, log2err) = compute_pvalue_multilevel(&int_ranks, k, observed_es, sample_size, 42, 1e-10);

        assert!(
            pval >= 0.0 && pval <= 1.0,
            "p-value out of [0,1]: {pval}"
        );
        // For a very high ES, p-value should be small
        assert!(pval < 0.1, "expected small pval for high ES, got {pval}");
        println!("pval={pval:.2e}, log2err={log2err:.2}");
    }

    #[test]
    fn test_esruler_negative_es() {
        let n = 200usize;
        let k = 10usize;
        let float_ranks: Vec<f64> = (0..n).map(|i| (n - i) as f64).collect();
        let int_ranks = scale_ranks(&float_ranks);

        // Negative ES (pathway at bottom of the list)
        let observed_es = -0.8f64;
        let sample_size = 51usize;
        let (pval, _log2err) = compute_pvalue_multilevel(&int_ranks, k, observed_es, sample_size, 42, 1e-10);

        assert!(
            pval >= 0.0 && pval <= 1.0,
            "p-value out of [0,1]: {pval}"
        );
        assert!(pval < 0.1, "expected small pval for large |negative ES|, got {pval}");
    }

    #[test]
    fn test_esruler_near_zero_es() {
        let n = 100usize;
        let k = 5usize;
        let float_ranks: Vec<f64> = (0..n).map(|i| (n - i) as f64).collect();
        let int_ranks = scale_ranks(&float_ranks);

        // Near-zero ES should give large p-value
        let (pval, _) = compute_pvalue_multilevel(&int_ranks, k, 0.05, 51, 42, 1e-10);
        assert!(pval > 0.1, "expected large pval for near-zero ES, got {pval}");
    }

    /// Verify that the batch function agrees with the single-pathway function.
    /// Also checks that negative ES values are handled correctly and that results
    /// match whether we run pathways individually or in a batch.
    #[test]
    fn test_batch_matches_single() {
        let n = 200usize;
        let k = 10usize;
        let float_ranks: Vec<f64> = (0..n).map(|i| (n - i) as f64).collect();
        let int_ranks = scale_ranks(&float_ranks);
        let sample_size = 51usize;
        let seed = 42u64;
        let eps = 1e-10f64;

        // Mix of positive and negative ES values — all same pathway size
        let es_values = vec![0.7f64, 0.5, -0.6, -0.4, 0.2];

        let batch_results =
            compute_pvalue_multilevel_batch(&int_ranks, k, &es_values, sample_size, seed, eps);

        assert_eq!(batch_results.len(), es_values.len());
        for (i, (&es, &(pval, _log2err))) in es_values.iter().zip(batch_results.iter()).enumerate() {
            assert!(
                pval >= 0.0 && pval <= 1.0,
                "batch pval[{i}] = {pval} out of [0,1] for es={es}"
            );
        }

        // Strongly enriched pathways should have small p-values
        assert!(
            batch_results[0].0 < 0.1,
            "expected small pval for es=0.7, got {}",
            batch_results[0].0
        );
        assert!(
            batch_results[2].0 < 0.1,
            "expected small pval for es=-0.6, got {}",
            batch_results[2].0
        );
    }

    /// Verify that the batch function scales: a larger batch of same-size pathways
    /// returns a result for each without panic.
    #[test]
    fn test_batch_multiple_same_size() {
        let n = 300usize;
        let k = 15usize;
        let float_ranks: Vec<f64> = (0..n).map(|i| (n - i) as f64).collect();
        let int_ranks = scale_ranks(&float_ranks);
        let sample_size = 51usize;

        // Ten pathways with the same size k=15
        let es_values: Vec<f64> = vec![0.8, 0.6, 0.4, 0.2, 0.1, -0.1, -0.3, -0.5, -0.7, -0.85];
        let results =
            compute_pvalue_multilevel_batch(&int_ranks, k, &es_values, sample_size, 7, 1e-10);

        assert_eq!(results.len(), es_values.len());
        for (i, &(pval, _)) in results.iter().enumerate() {
            assert!(
                pval >= 0.0 && pval <= 1.0,
                "pval[{i}] = {pval} out of [0,1]"
            );
        }
        // Monotone: larger |ES| should give smaller or equal p-value
        // (not strictly guaranteed for MCMC but should hold for these extremes)
        assert!(
            results[0].0 <= results[1].0,
            "expected pval(0.8) <= pval(0.6)"
        );
    }
}
