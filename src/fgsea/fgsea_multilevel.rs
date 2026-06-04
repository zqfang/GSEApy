//! Translated from `fgsea/src/fgseaMultilevel.cpp` + `fgsea/src/fgseaMultilevel.h`.
//!
//! `ranks` arrives as an R integer vector (`INTSXP`); modelled as `&[i32]`.
//! The Rcpp `DataFrame` return becomes the `FgseaMultilevelResult` struct.
//!
//! Vendored from the faithful fgsea-rs translation.

use crate::fgsea::fgsea_multilevel_supplement::EsRuler;

/// Result of [`fgseaMultilevelCpp`], standing in for the Rcpp `DataFrame` returned
/// by the original C++. Each field is a column of the table of p-values and
/// estimation errors:
///
/// * `cppMPval` - the estimated GSEA p-values.
/// * `cppIsCpGeHalf` - whether the conditional probability is greater than or equal
///   to one half (legacy flag carried over from the C++ tuple).
/// * `log2err` - the estimation error of each p-value, in log2 units.
pub struct FgseaMultilevelResult {
    pub cppMPval: Vec<f64>,
    pub cppIsCpGeHalf: Vec<bool>,
    pub log2err: Vec<f64>,
}

/// Calculates low GSEA p-values for a given gene set size using the multilevel split
/// Monte Carlo approach.
///
/// A positive and a negative [`EsRuler`] are built and extended up to the most
/// extreme enrichment scores; each score is then converted to a p-value via the
/// appropriate ruler depending on its sign.
///
/// # Arguments
/// * `enrichmentScores` - the enrichment scores for which p-values should be calculated.
/// * `ranks` - the gene-level statistics, a 1-based integer vector; absolute values
///   are taken internally. (In C++ this arrives as an `INTSXP`; the Rust type system
///   already enforces an integer vector, so the original runtime type check is dropped.)
/// * `pathwaySize` - the size of the gene set.
/// * `sampleSize` - number of samples per level (expected to be odd).
/// * `seed` - random seed.
/// * `eps` - p-values below `eps` are not calculated; pass `0` to disable the early exit.
/// * `sign` - controls whether the ES^+ or the ES score is used.
/// * `moveScale` - controls the number of MCMC iterations on each level.
/// * `logStatus` - controls whether debug output is shown.
///
/// # Returns
/// A [`FgseaMultilevelResult`] table with p-values (`cppMPval`), the conditional
/// probability flag (`cppIsCpGeHalf`), and estimation errors (`log2err`).
///
/// C++: `DataFrame fgseaMultilevelCpp(const NumericVector& enrichmentScores, const SEXP& ranks, int pathwaySize, int sampleSize, int seed, double eps, bool sign, double moveScale = 1.0, bool logStatus = false)`
pub fn fgseaMultilevelCpp(
    enrichmentScores: &[f64],
    ranks: &[i32],
    pathwaySize: i32,
    sampleSize: i32,
    seed: i32,
    eps: f64,
    sign: bool,
    moveScale: f64,
    logStatus: bool,
) -> FgseaMultilevelResult {
    // C++ checks TYPEOF(ranks) == INTSXP and stop()s otherwise; the Rust type
    // system enforces the integer ranks vector, so the check is unnecessary.
    let mut posRanks: Vec<i64> = ranks.iter().map(|&r| r as i64).collect();
    for i in 0..posRanks.len() {
        posRanks[i] = posRanks[i].abs();
    }
    let mut negRanks: Vec<i64> = posRanks.clone();
    negRanks.reverse();

    let esVector: Vec<f64> = enrichmentScores.to_vec();

    let mut esRulerPos = EsRuler::new(
        &posRanks,
        sampleSize as u32,
        pathwaySize as u32,
        moveScale,
        logStatus,
    );
    let mut esRulerNeg = EsRuler::new(
        &negRanks,
        sampleSize as u32,
        pathwaySize as u32,
        moveScale,
        logStatus,
    );

    let maxES = *esVector
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let minES = *esVector
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    if maxES >= 0.0 {
        esRulerPos.extend(maxES.abs(), seed, eps);
    }
    if minES < 0.0 {
        esRulerNeg.extend(minES.abs(), seed, eps);
    }

    let mut pvalRes: Vec<f64> = Vec::new();
    let mut isCpGeHalf: Vec<bool> = Vec::new();
    let mut log2err: Vec<f64> = Vec::new();

    let nrow = esVector.len();
    for i in 0..nrow {
        let currentES = esVector[i];
        let res = if currentES >= 0.0 {
            esRulerPos.getPvalue(currentES.abs(), eps, sign)
        } else {
            esRulerNeg.getPvalue(currentES.abs(), eps, sign)
        };

        pvalRes.push(res.0);
        isCpGeHalf.push(res.1);
        log2err.push(res.2);
    }

    // return vector with pvalues and vector with conditional probability result
    FgseaMultilevelResult {
        cppMPval: pvalRes,
        cppIsCpGeHalf: isCpGeHalf,
        log2err,
    }
}
