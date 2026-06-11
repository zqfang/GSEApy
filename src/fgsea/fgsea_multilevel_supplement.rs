//! Translated from `fgsea/src/fgseaMultilevelSupplement.cpp` + `fgsea/src/fgseaMultilevelSupplement.h`.
//!
//! `std::function<bool(int,int)>` in `perturbate_until` becomes a generic closure
//! parameter `F: Fn(i32, i32) -> bool`. The C++ `check` lambda inside
//! `perturbate_until` is kept as a local closure (a faithful 1:1 of the original
//! lambda, not an added helper).
//!
//! Vendored from the faithful fgsea-rs translation.

use special::Gamma;

use crate::fgsea::es_calculation::{calcES, calcES_auto, calcPositiveES, calcPositiveES_auto, score_t};
use crate::fgsea::util::{betaMeanLog, combination, random_engine_t, uid_wrapper};

pub type hash_t = u64;
pub type gsea_t = (score_t, hash_t);

/// One level of the multilevel split: the samples partitioned around the level's
/// bound into those at or below it (`lowScores`) and those strictly above it
/// (`highScores`), each tagged with whether the underlying ES is positive.
///
/// C++: `struct Level` (nested in `EsRuler`).
pub struct Level {
    pub lowScores: Vec<(gsea_t, bool)>,
    pub highScores: Vec<(gsea_t, bool)>,
    pub bound: gsea_t,
}

/// Outcome of a perturbation run: the number of accepted `moves` and the number of
/// attempted `iters`.
///
/// C++: `struct PerturbateResult` (nested in `EsRuler`).
pub struct PerturbateResult {
    pub moves: i32,
    pub iters: i32,
}

/// A single sample's gene indices split into contiguous chunks, with the per-chunk
/// rank sums cached in `chunkSum` to speed up score recomputation during perturbation.
///
/// C++: `struct SampleChunks` (nested in `EsRuler`).
pub struct SampleChunks {
    pub chunkSum: Vec<i64>,
    pub chunks: Vec<Vec<i32>>,
}

impl SampleChunks {
    /// Creates empty chunks and zeroed chunk sums for `chunksNumber` chunks.
    ///
    /// C++: `EsRuler::SampleChunks::SampleChunks(int chunksNumber)`
    pub fn new(chunksNumber: i32) -> SampleChunks {
        SampleChunks {
            chunkSum: vec![0i64; chunksNumber as usize],
            chunks: vec![Vec::new(); chunksNumber as usize],
        }
    }
}

/// The multilevel split Monte-Carlo engine.
///
/// `EsRuler` builds an empirical "ruler" of enrichment-score levels by repeatedly
/// sampling random gene sets of a fixed size, splitting them at their median score,
/// and perturbating the surviving high-scoring samples upward. The resulting ladder
/// of levels lets [`EsRuler::getPvalue`] estimate very small GSEA p-values together
/// with their log2 error.
///
/// C++: `class EsRuler`
pub struct EsRuler {
    logStatus: bool,
    ranks: Vec<i64>,
    geneHashes: Vec<hash_t>,
    sampleSize: u32,
    pathwaySize: u32,
    movesScale: f64,
    incorrectRuler: bool,
    currentSamples: Vec<Vec<i32>>,
    oldSamplesStart: i32,
    levels: Vec<Level>,
    chunkLastElement: Vec<i32>,
    chunksNumber: i32,
}

impl EsRuler {
    /// Constructs an empty ruler over the given `inpRanks`, sized for `inpSampleSize`
    /// samples per level and a gene set of `inpPathwaySize`. No sampling happens here;
    /// call [`EsRuler::extend`] to build the levels.
    ///
    /// C++: `EsRuler::EsRuler(const vector<int64_t>&, unsigned, unsigned, double, bool)`
    pub fn new(
        inpRanks: &[i64],
        inpSampleSize: u32,
        inpPathwaySize: u32,
        inpMovesScale: f64,
        inpLog: bool,
    ) -> EsRuler {
        EsRuler {
            logStatus: inpLog,
            ranks: inpRanks.to_vec(),
            geneHashes: vec![0u64; inpRanks.len()],
            sampleSize: inpSampleSize,
            pathwaySize: inpPathwaySize,
            movesScale: inpMovesScale,
            incorrectRuler: false,
            currentSamples: vec![Vec::new(); inpSampleSize as usize],
            oldSamplesStart: 0,
            levels: Vec::new(),
            chunkLastElement: Vec::new(),
            chunksNumber: 0,
        }
    }

    /// Rescales floating-point `ranks` to integers, multiplying by the largest factor
    /// that keeps their sum within `score_t::getMaxNS()` (rounded down so input
    /// integers stay integer).
    ///
    /// C++: `vector<int64_t> EsRuler::scaleRanks(vector<double> const& ranks)`
    fn scaleRanks(&self, ranks: &[f64]) -> Vec<i64> {
        let MAX_NS: i64 = score_t::getMaxNS();
        let curNS: f64 = ranks.iter().sum();
        let mut result = vec![0i64; ranks.len()];
        // keep input integers integer
        let scale = (MAX_NS as f64 / curNS).floor();
        if self.logStatus {
            eprintln!("Scaling ranks by {}", scale);
        }
        for i in 0..ranks.len() {
            result[i] = (ranks[i] * scale) as i64;
        }
        result
    }

    /// Sorts the current samples by score, records a new [`Level`] split at the median,
    /// and resamples the lower half by drawing replacements uniformly from the upper
    /// half. Returns `true` on success; the special "all equal values" case ends the
    /// multilevel process without adding a level (also returning `true`).
    ///
    /// C++: `bool EsRuler::resampleGenesets(random_engine_t& rng)`
    fn resampleGenesets(&mut self, rng: &mut random_engine_t) -> bool {
        let sampleSize = self.sampleSize as usize;
        //  <score, is_positive, ind>
        let mut stats: Vec<(gsea_t, i32, i32)> = Vec::with_capacity(sampleSize);

        for sampleId in 0..sampleSize {
            let sampleEsPos = calcPositiveES_auto(&self.ranks, &self.currentSamples[sampleId]);
            let sampleEs = calcES_auto(&self.ranks, &self.currentSamples[sampleId]);
            let sampleHash = self.calcHash(&self.currentSamples[sampleId]);
            stats.push((
                (sampleEsPos, sampleHash),
                (sampleEs.getNumerator() >= 0) as i32,
                sampleId as i32,
            ));
        }
        stats.sort();

        let mut startFrom: usize = 0;
        let centralValue = stats[sampleSize / 2].0;
        for sampleId in 0..sampleSize {
            if stats[sampleId].0 >= centralValue {
                startFrom = sampleId;
                break;
            }
        }

        if startFrom == 0 {
            while startFrom < sampleSize && stats[startFrom].0 == stats[0].0 {
                startFrom += 1;
            }
        }

        if startFrom == sampleSize {
            if self.logStatus {
                eprintln!("Got all equal values. Ending multilevel process");
            }
            return true;
        }

        let mut level = Level {
            lowScores: Vec::new(),
            highScores: Vec::new(),
            bound: stats[startFrom - 1].0, //  greater
        };
        for i in 0..startFrom {
            level.lowScores.push((stats[i].0, stats[i].1 != 0));
        }
        for i in startFrom..sampleSize {
            level.highScores.push((stats[i].0, stats[i].1 != 0));
        }
        self.levels.push(level);

        let uid = uid_wrapper::new(0, (sampleSize - startFrom - 1) as i32, rng);

        // gen_new_sample lambda (inlined): int ind = uid() + startFrom;
        //                                  return currentSamples[get<2>(stats[ind])];
        let mut new_sets: Vec<Vec<i32>> = Vec::new();
        for _ in 0..startFrom {
            let ind = uid.call(rng) as usize + startFrom;
            new_sets.push(self.currentSamples[stats[ind].2 as usize].clone());
        }
        for i in startFrom..sampleSize {
            new_sets.push(self.currentSamples[stats[i].2 as usize].clone());
        }

        self.oldSamplesStart = startFrom as i32;
        std::mem::swap(&mut self.currentSamples, &mut new_sets);
        true
    }

    /// Computes the XOR of the per-gene hashes for the genes in `curSample`, giving an
    /// order-independent fingerprint used to break ties between equal scores.
    ///
    /// C++: `EsRuler::hash_t EsRuler::calcHash(const vector<int>& curSample)`
    fn calcHash(&self, curSample: &[i32]) -> hash_t {
        let mut res: hash_t = 0;
        for &i in curSample {
            res ^= self.geneHashes[i as usize];
        }
        res
    }

    /// Builds the ladder of levels up to the target enrichment score `ES_double`.
    ///
    /// Starting from random gene sets, it repeatedly resamples and perturbates the
    /// high-scoring samples (running enough MCMC iterations, scaled by `movesScale`)
    /// until the latest level's bound reaches `ES_double` or progress stalls. When
    /// `eps` is non-zero the build exits early once the accumulated log p-value drops
    /// below `ln(eps)`.
    ///
    /// The target score is floored to an exact `score_t` ("flooring ES an exact score,
    /// but 1.0 should stay 1.0"; works best when `getMaxNS()` is a power of two).
    ///
    /// C++: `void EsRuler::extend(double ES_double, int seed, double eps)`
    pub fn extend(&mut self, ES_double: f64, seed: i32, eps: f64) {
        let mut gen = random_engine_t::new(seed as u32);
        let n = self.ranks.len() as i32;
        let k = self.pathwaySize as i32;
        let sampleSize = self.sampleSize as usize;

        for i in 0..n {
            self.geneHashes[i as usize] = gen.next_u32() as hash_t;
        }

        for sampleId in 0..sampleSize {
            self.currentSamples[sampleId] =
                combination(0, self.ranks.len() as i32 - 1, self.pathwaySize as i32, &mut gen);
            self.currentSamples[sampleId].sort();
        }

        if !self.resampleGenesets(&mut gen) {
            if self.logStatus {
                eprintln!("Could not advance in the start");
            }
            self.incorrectRuler = true;
            return;
        }

        self.chunksNumber = std::cmp::max(1, (self.pathwaySize as f64).sqrt() as i32);
        self.chunkLastElement = vec![0i32; self.chunksNumber as usize];
        self.chunkLastElement[(self.chunksNumber - 1) as usize] = self.ranks.len() as i32;
        let mut tmp = vec![0i32; sampleSize];
        let mut samplesChunks: Vec<SampleChunks> = (0..sampleSize)
            .map(|_| SampleChunks::new(self.chunksNumber))
            .collect();

        // flooring ES an exact score, but 1.0 should stay 1.0
        // works best if getMaxNS is a power of two
        let NEED_ES = score_t {
            NS: score_t::getMaxNS(),
            coef_NS: (score_t::getMaxNS() as f64 * ES_double) as i64,
            diff: 1,
            coef_const: 0,
        };

        let mut adjLogPval = 0.0;
        let mut levelNum = 1;
        while self.levels.last().unwrap().bound.0 < NEED_ES {
            adjLogPval += betaMeanLog(
                (self.levels.last().unwrap().highScores.len() + 1) as u64,
                self.sampleSize as u64,
            );
            if eps != 0.0 && adjLogPval < eps.ln() {
                break;
            }

            if self.logStatus {
                eprintln!(
                    "Iteration {}: score={:.15}, hash={}",
                    levelNum,
                    self.levels.last().unwrap().bound.0.getDouble(),
                    self.levels.last().unwrap().bound.1
                );
            }

            let mut pos = 0i32;
            for i in 0..(self.chunksNumber - 1) {
                pos += (self.pathwaySize as i32 + i) / self.chunksNumber;
                for j in 0..sampleSize {
                    tmp[j] = self.currentSamples[j][pos as usize];
                }
                let mid = sampleSize / 2;
                tmp.select_nth_unstable(mid);
                self.chunkLastElement[i as usize] = tmp[mid];
            }

            for i in 0..sampleSize {
                for s in samplesChunks[i].chunkSum.iter_mut() {
                    *s = 0;
                }
                for j in 0..self.chunksNumber as usize {
                    samplesChunks[i].chunks[j].clear();
                }
                let mut cnt = 0usize;
                for idx in 0..self.currentSamples[i].len() {
                    let pos = self.currentSamples[i][idx];
                    while self.chunkLastElement[cnt] <= pos {
                        cnt += 1;
                    }
                    samplesChunks[i].chunks[cnt].push(pos);
                    samplesChunks[i].chunkSum[cnt] += self.ranks[pos as usize];
                }
            }

            let mut nIterations = 0;
            let mut nAccepted = 0;
            let needAccepted =
                (self.movesScale * self.sampleSize as f64 * self.pathwaySize as f64 / 2.0) as i32;
            while nAccepted < needAccepted {
                for sampleId in 0..sampleSize {
                    let bound = self.levels.last().unwrap().bound;
                    let perturbResult =
                        self.perturbate(&self.ranks, k, &mut samplesChunks[sampleId], bound, &mut gen);
                    nAccepted += perturbResult.moves;
                }
                nIterations += 1;
            }
            for _ in 0..nIterations {
                for sampleId in 0..sampleSize {
                    let bound = self.levels.last().unwrap().bound;
                    self.perturbate(&self.ranks, k, &mut samplesChunks[sampleId], bound, &mut gen);
                }
            }

            for i in 0..sampleSize {
                self.currentSamples[i].clear();
                for j in 0..self.chunksNumber as usize {
                    for idx in 0..samplesChunks[i].chunks[j].len() {
                        let pos = samplesChunks[i].chunks[j][idx];
                        self.currentSamples[i].push(pos);
                    }
                }
            }

            let lastSize = self.levels.len();
            if !self.resampleGenesets(&mut gen) {
                self.incorrectRuler = true;
                if self.logStatus {
                    eprintln!("Could not advance after level {}", levelNum);
                }
            }
            if lastSize == self.levels.len() {
                break;
            }
            levelNum += 1;
        }
    }

    /// Estimates `Pr[score >= ES]` from the built ruler, walking the levels and
    /// accumulating the mean log conditional probability along with the per-level
    /// variance used for the error estimate.
    ///
    /// The target score is floored to an exact `score_t` ("flooring ES an exact score,
    /// but 1.0 should stay 1.0"; works best when `getMaxNS()` is a power of two).
    /// `sign` selects whether the counts use the signed ES (`cntLast`) or only the
    /// positive ES (`cntPositive`). If the ruler is incorrect, `(NaN, true, NaN)` is
    /// returned.
    ///
    /// # Returns
    /// A tuple `(p-value, true, error)`: the estimated p-value clamped to `[0, 1]`, a
    /// legacy boolean flag (always `true`), and the log2 estimation error (`NaN` when
    /// the numerator is zero).
    ///
    /// C++: `tuple<double, bool, double> EsRuler::getPvalue(double ES_double, double eps, bool sign)`
    pub fn getPvalue(&self, ES_double: f64, eps: f64, sign: bool) -> (f64, bool, f64) {
        if self.incorrectRuler {
            return (f64::NAN, true, f64::NAN);
        }

        // flooring ES an exact score, but 1.0 should stay 1.0
        // works best if getMaxNS is a power of two
        let ES_score = score_t {
            NS: score_t::getMaxNS(),
            coef_NS: (score_t::getMaxNS() as f64 * ES_double) as i64,
            diff: 1,
            coef_const: 0,
        };
        let ES: gsea_t = (ES_score, 0);
        //  Calculate Pr[<score, hash> >= ES]

        let mut adjLogPval = 0.0;
        let mut lvlsVar = 0.0;

        for lvl in &self.levels {
            if ES <= lvl.bound {
                let mut cntLast = 0;
                let mut cntPositive = 0;
                //  highScores > lvl.bound >= ES
                for &(_x, isPositive) in &lvl.highScores {
                    cntLast += 1;
                    cntPositive += isPositive as i32;
                }
                for &(x, isPositive) in &lvl.lowScores {
                    if x >= ES {
                        cntLast += 1;
                        cntPositive += isPositive as i32;
                    }
                }

                let numerator = if sign { cntLast } else { cntPositive };
                if numerator == 0 {
                    adjLogPval += betaMeanLog(1, self.sampleSize as u64);
                    return (
                        0.0_f64.max(1.0_f64.min(adjLogPval.exp())),
                        true,
                        f64::NAN,
                    );
                }

                adjLogPval += betaMeanLog(numerator as u64, self.sampleSize as u64);
                lvlsVar += getVarPerLevel(numerator as u64, self.sampleSize as u64);

                let log2err = lvlsVar.sqrt() / (2.0_f64).ln();
                return (0.0_f64.max(1.0_f64.min(adjLogPval.exp())), true, log2err);
            }

            let mut nhigh = lvl.highScores.len() as i32;
            nhigh += 1;
            adjLogPval += betaMeanLog(nhigh as u64, self.sampleSize as u64);
            lvlsVar += getVarPerLevel(nhigh as u64, self.sampleSize as u64);
        }

        let lastLevel = self.levels.last().unwrap();
        //  ES > lastLevel.bound
        //  highScores -> uniform [> lastLevel.bound]
        let mut cntLast = 0;
        let mut cntPositive = 0;
        for &(x, isPositive) in &lastLevel.highScores {
            if x >= ES {
                cntLast += 1;
                cntPositive += isPositive as i32;
            }
        }

        let numerator = if sign { cntLast } else { cntPositive };

        if numerator == 0 {
            adjLogPval += betaMeanLog(1, lastLevel.highScores.len() as u64);
            return (0.0_f64.max(1.0_f64.min(adjLogPval.exp())), true, f64::NAN);
        }

        adjLogPval += betaMeanLog(numerator as u64, lastLevel.highScores.len() as u64);
        lvlsVar += getVarPerLevel(numerator as u64, lastLevel.highScores.len() as u64);

        let log2err = lvlsVar.sqrt() / (2.0_f64).ln();
        (0.0_f64.max(1.0_f64.min(adjLogPval.exp())), true, log2err)
    }

    /// Returns the length of chunk `ind` (the span between consecutive
    /// `chunkLastElement` boundaries).
    ///
    /// C++: `int EsRuler::chunkLen(int ind)`
    fn chunkLen(&self, ind: i32) -> i32 {
        if ind == 0 {
            return self.chunkLastElement[0];
        }
        self.chunkLastElement[ind as usize] - self.chunkLastElement[(ind - 1) as usize]
    }

    /// Runs a default-sized perturbation pass on one sample: `max(1, k * 0.1)`
    /// iterations of swapping a member gene for a random non-member while keeping the
    /// sample's score at or above `bound`.
    ///
    /// C++: `EsRuler::PerturbateResult EsRuler::perturbate(...)`
    fn perturbate(
        &self,
        ranks: &[i64],
        k: i32,
        sampleChunks: &mut SampleChunks,
        bound: gsea_t,
        rng: &mut random_engine_t,
    ) -> PerturbateResult {
        let pertPrmtr = 0.1;
        let iters = std::cmp::max(1, (k as f64 * pertPrmtr) as i32);
        self.perturbate_iters(ranks, k, sampleChunks, bound, rng, iters)
    }

    /// Perturbates until at least `need_iters` iterations have been attempted.
    ///
    /// C++: `EsRuler::PerturbateResult EsRuler::perturbate_iters(...)`
    fn perturbate_iters(
        &self,
        ranks: &[i64],
        k: i32,
        sampleChunks: &mut SampleChunks,
        bound: gsea_t,
        rng: &mut random_engine_t,
        need_iters: i32,
    ) -> PerturbateResult {
        self.perturbate_until(ranks, k, sampleChunks, bound, rng, move |_moves, iters| {
            iters >= need_iters
        })
    }

    /// Perturbates until at least `need_successes` moves have been accepted.
    ///
    /// C++: `EsRuler::PerturbateResult EsRuler::perturbate_success(...)`
    fn perturbate_success(
        &self,
        ranks: &[i64],
        k: i32,
        sampleChunks: &mut SampleChunks,
        bound: gsea_t,
        rng: &mut random_engine_t,
        need_successes: i32,
    ) -> PerturbateResult {
        self.perturbate_until(ranks, k, sampleChunks, bound, rng, move |moves, _iters| {
            moves >= need_successes
        })
    }

    /// Core MCMC perturbation loop, shared by the iter- and success-bounded variants.
    ///
    /// On each iteration it picks a member gene to drop and a random gene to add; if
    /// the swap keeps the sample's score at or above `bound` (compared strictly or
    /// non-strictly depending on the hash tie-break), the move is accepted and the
    /// chunked representation, running rank sum `NS`, and hash are updated in place,
    /// otherwise the swap is rolled back. A cached candidate split point lets most
    /// acceptance checks short-circuit. The loop runs until the stopping predicate `f`,
    /// applied to `(moves, iters)`, returns `true`.
    ///
    /// C++: `EsRuler::PerturbateResult EsRuler::perturbate_until(..., std::function<bool(int,int)> const& f)`
    fn perturbate_until<F: Fn(i32, i32) -> bool>(
        &self,
        ranks: &[i64],
        k: i32,
        sampleChunks: &mut SampleChunks,
        bound: gsea_t,
        rng: &mut random_engine_t,
        f: F,
    ) -> PerturbateResult {
        let n = ranks.len() as i32;
        let uid_n = uid_wrapper::new(0, n - 1, rng);
        let uid_k = uid_wrapper::new(0, k - 1, rng);

        let mut NS: i64 = 0;
        let mut curHash: hash_t = 0;
        for i in 0..sampleChunks.chunks.len() {
            for &pos in &sampleChunks.chunks[i] {
                NS += ranks[pos as usize];
                curHash ^= self.geneHashes[pos as usize];
            }
        }
        let mut candVal: i32 = -1;
        let mut hasCand = false;
        let mut candX: i32 = 0;
        let mut candY: i64 = 0;

        let mut moves = 0;
        let mut iters = 0;
        while !f(moves, iters) {
            iters += 1;
            let oldInd = uid_k.call(rng);

            let mut oldChunkInd = 0i32;
            let oldIndInChunk: i32;
            let oldVal: i32;
            {
                let mut t = oldInd;
                while (sampleChunks.chunks[oldChunkInd as usize].len() as i32) <= t {
                    t -= sampleChunks.chunks[oldChunkInd as usize].len() as i32;
                    oldChunkInd += 1;
                }
                oldIndInChunk = t;
                oldVal = sampleChunks.chunks[oldChunkInd as usize][oldIndInChunk as usize];
            }

            let newVal = uid_n.call(rng);

            // upper_bound / lower_bound
            let newChunkInd =
                self.chunkLastElement.partition_point(|&x| x <= newVal) as i32;
            let newIndInChunk = sampleChunks.chunks[newChunkInd as usize]
                .partition_point(|&x| x < newVal) as i32;

            if newIndInChunk < sampleChunks.chunks[newChunkInd as usize].len() as i32
                && sampleChunks.chunks[newChunkInd as usize][newIndInChunk as usize] == newVal
            {
                if newVal == oldVal {
                    moves += 1;
                }
                continue;
            }

            sampleChunks.chunks[oldChunkInd as usize].remove(oldIndInChunk as usize);
            let ins = newIndInChunk
                - if oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk {
                    1
                } else {
                    0
                };
            sampleChunks.chunks[newChunkInd as usize].insert(ins as usize, newVal);

            NS = NS - ranks[oldVal as usize] + ranks[newVal as usize];
            curHash ^= self.geneHashes[oldVal as usize] ^ self.geneHashes[newVal as usize];
            sampleChunks.chunkSum[oldChunkInd as usize] -= ranks[oldVal as usize];
            sampleChunks.chunkSum[newChunkInd as usize] += ranks[newVal as usize];

            let strictly = curHash <= bound.1;

            let check = |score: score_t| -> bool {
                if strictly {
                    score > bound.0
                } else {
                    score >= bound.0
                }
            };

            if hasCand {
                if oldVal == candVal {
                    hasCand = false;
                }
            }

            if hasCand {
                if oldVal < candVal {
                    candX += 1;
                    candY -= ranks[oldVal as usize];
                }
                if newVal < candVal {
                    candX -= 1;
                    candY += ranks[newVal as usize];
                }
            }

            if hasCand
                && check(score_t {
                    NS,
                    coef_NS: candY,
                    diff: (n - k) as i64,
                    coef_const: candX as i64,
                })
            {
                moves += 1;
                continue;
            }

            let mut curX = 0i32;
            let mut curY: i64 = 0;
            let mut ok = false;
            let mut last = -1i32;

            for i in 0..sampleChunks.chunks.len() {
                if !check(score_t {
                    NS,
                    coef_NS: curY + sampleChunks.chunkSum[i],
                    diff: (n - k) as i64,
                    coef_const: curX as i64,
                }) {
                    curY += sampleChunks.chunkSum[i];
                    curX += self.chunkLastElement[i] - last - 1
                        - sampleChunks.chunks[i].len() as i32;
                    last = self.chunkLastElement[i] - 1;
                } else {
                    for idx in 0..sampleChunks.chunks[i].len() {
                        let pos = sampleChunks.chunks[i][idx];
                        curY += ranks[pos as usize];
                        curX += pos - last - 1;
                        if check(score_t {
                            NS,
                            coef_NS: curY,
                            diff: (n - k) as i64,
                            coef_const: curX as i64,
                        }) {
                            ok = true;
                            hasCand = true;
                            candX = curX;
                            candY = curY;
                            candVal = pos;
                            break;
                        }
                        last = pos;
                    }
                    if ok {
                        break;
                    }
                    curX += self.chunkLastElement[i] - 1 - last;
                    last = self.chunkLastElement[i] - 1;
                }
            }

            if !ok {
                NS = NS - ranks[newVal as usize] + ranks[oldVal as usize];
                curHash ^= self.geneHashes[newVal as usize] ^ self.geneHashes[oldVal as usize];

                sampleChunks.chunkSum[oldChunkInd as usize] += ranks[oldVal as usize];
                sampleChunks.chunkSum[newChunkInd as usize] -= ranks[newVal as usize];

                let rem = newIndInChunk
                    - if oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk {
                        1
                    } else {
                        0
                    };
                sampleChunks.chunks[newChunkInd as usize].remove(rem as usize);
                sampleChunks.chunks[oldChunkInd as usize].insert(oldIndInChunk as usize, oldVal);

                if hasCand {
                    if newVal == candVal {
                        hasCand = false;
                    }
                }
                if hasCand {
                    if oldVal < candVal {
                        candX -= 1;
                        candY += ranks[oldVal as usize];
                    }
                    if newVal < candVal {
                        candX += 1;
                        candY -= ranks[newVal as usize];
                    }
                }
            } else {
                moves += 1;
            }
        }
        PerturbateResult { moves, iters }
    }
}

/// Per-level variance contribution to the p-value estimate, computed as
/// `trigamma(k) - trigamma(n + 1)`.
///
/// C++: `double getVarPerLevel(unsigned long k, unsigned long n)`
pub fn getVarPerLevel(k: u64, n: u64) -> f64 {
    (k as f64).trigamma() - ((n + 1) as f64).trigamma()
}

/// Computes a single log conditional-probability correction term from the
/// `probCorrector` table, returning it alongside a flag indicating whether the
/// corresponding conditional probability is at least one half.
///
/// C++: `pair<double, bool> calcLogCorrection(const vector<unsigned int>& probCorrector, long probCorrIndx, unsigned int sampleSize)`
pub fn calcLogCorrection(probCorrector: &[u32], probCorrIndx: i64, sampleSize: u32) -> (f64, bool) {
    let mut result = 0.0f64;

    let halfSize: u64 = ((sampleSize + 1) / 2) as u64;
    let remainder: u64 = sampleSize as u64 - (probCorrIndx as u64) % halfSize;

    let condProb = betaMeanLog((probCorrector[probCorrIndx as usize] as u64) + 1, remainder);
    result += condProb;

    if condProb.exp() >= 0.5 {
        (result, true)
    } else {
        (result, false)
    }
}
