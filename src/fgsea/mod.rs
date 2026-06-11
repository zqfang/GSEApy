//! Faithful Rust translation of the `fgsea` C++ core (multilevel p-value +
//! calcGseaStat batch), vendored from the standalone fgsea-rs port.
//!
//! Source: `fgsea/src/*.cpp` + `*.h` (alserglab/fgsea). Each original C++ function
//! maps to exactly one Rust function; original camelCase names are preserved to keep
//! the audit 1:1, hence the module-wide allowances below.
//!
//! Module ↔ source mapping:
//!   util.rs                        <- util.cpp/.h
//!   es_calculation.rs              <- esCalculation.cpp/.h
//!   fast_gsea.rs                   <- fastGSEA.cpp/.h
//!   fgsea_multilevel_supplement.rs <- fgseaMultilevelSupplement.cpp/.h
//!   fgsea_multilevel.rs            <- fgseaMultilevel.cpp/.h
//!
//! Out of scope for GSEApy (geseca / ScoreRuler / ScoreCalculation are not ported).
//!
//! [`compute_pvalue_multilevel`] and [`scale_ranks`] are thin GSEApy-facing adapters
//! over the faithful core; `crate::stats` consumes exactly these two symbols, so the
//! prerank multilevel driver stays unchanged.

#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(unused_variables)]
#![allow(unused_imports)]

pub mod es_calculation;
pub mod fast_gsea;
pub mod fgsea_multilevel;
pub mod fgsea_multilevel_supplement;
pub mod util;

#[cfg(test)]
mod validation_tests;

// Faithful C++ entry points (the functions Rcpp exported to R).
pub use fast_gsea::{calcGseaStatBatchCpp, calcGseaStatCumulative, calcGseaStatCumulativeBatch};
pub use fgsea_multilevel::fgseaMultilevelCpp;

use es_calculation::score_t;

/// Scales a non-negative float metric to i64 integers whose sum stays within
/// `score_t::getMaxNS()` (`1 << 30`), matching fgsea's `EsRuler::scaleRanks`.
///
/// The scale factor is floored (so integer inputs stay integer). Two safety guards
/// beyond the C++ original handle degenerate inputs that fgsea never feeds in
/// practice: an all-zero metric returns zeros, and a tiny sum (scale < 1) is clamped
/// to 1 so individual elements never overflow.
///
/// Kept signature-compatible with `crate::stats`'s import.
pub fn scale_ranks(ranks: &[f64]) -> Vec<i64> {
    let max_ns: i64 = score_t::getMaxNS();
    let cur_ns: f64 = ranks.iter().sum();
    if cur_ns == 0.0 {
        return vec![0i64; ranks.len()];
    }
    let mut scale = (max_ns as f64 / cur_ns).floor();
    if scale < 1.0 {
        scale = 1.0;
    }
    ranks.iter().map(|&r| (r * scale) as i64).collect()
}

/// Computes the fgsea multilevel p-value and log2 error for one preranked gene set.
///
/// Adapter over the faithful [`fgseaMultilevelCpp`], preserving the exact signature
/// the prerank multilevel driver (`crate::stats::GSEAResult::prerank_multilevel`)
/// already calls.
///
/// * `int_ranks` — integer-scaled gene metric for ALL genes, in rank order
///   (descending), as produced by [`scale_ranks`].
/// * `pathway_size` — gene set size `k`.
/// * `observed_es` — the signed enrichment score of the gene set.
/// * `sample_size` — MCMC samples per level (fgsea default 101, should be odd).
/// * `seed` — random seed.
/// * `eps` — convergence threshold (fgsea default 1e-50; 0 disables early exit).
///
/// Returns `(p_value, log2err)`. Sign handling, abs of ranks, and positive/negative
/// ruler selection all happen inside `fgseaMultilevelCpp`.
pub fn compute_pvalue_multilevel(
    int_ranks: &[i64],
    pathway_size: usize,
    observed_es: f64,
    sample_size: usize,
    seed: u64,
    eps: f64,
) -> (f64, f64) {
    if observed_es.abs() < 1e-15 || pathway_size == 0 || pathway_size >= int_ranks.len() {
        return (1.0, f64::NAN);
    }

    // Scaled ranks sum to <= 1<<30, so every element fits in i32; fgseaMultilevelCpp
    // takes the R INTSXP rank vector (modelled as &[i32]) and casts back to i64.
    let ranks_i32: Vec<i32> = int_ranks.iter().map(|&r| r as i32).collect();

    let res = fgseaMultilevelCpp(
        &[observed_es],
        &ranks_i32,
        pathway_size as i32,
        sample_size as i32,
        seed as i32,
        eps,
        false, // sign: use the positive-ES conditional count (fgsea default)
        1.0,   // moveScale
        false, // logStatus
    );

    (res.cppMPval[0], res.log2err[0])
}
