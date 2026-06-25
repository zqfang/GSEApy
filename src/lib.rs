use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::HashMap;
// import own modules
mod algorithm;
mod fgsea;
mod gsva;
mod stats;
mod utils;
// export module fn, struct, trait ...
use algorithm::GseaStatResult;
use gsva::gsva;
use stats::{GSEAResult, GSEASummary};
use utils::{CorrelType, Metric, ScoreType};

/// Build a dedicated, per-call rayon thread pool sized to `threads`.
///
/// We intentionally do **not** touch the global pool or the `RAYON_NUM_THREADS`
/// environment variable. The env var is read only once, when rayon lazily
/// initializes its global pool, so mutating it from here would silently ignore
/// `threads` on every call after the first in a long-lived process (e.g. a
/// Python console). It is also process-global state and racy to mutate while
/// other threads read the environment. A fresh local pool, driven via
/// `pool.install(..)`, honors `threads` on every call and keeps the parallelism
/// fully scoped to this invocation. Nested `par_iter` calls inside the closure
/// automatically use the installed pool.
fn build_pool(threads: usize) -> PyResult<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyRuntimeError::new_err(format!("failed to build thread pool: {e}")))
}

/// Prerank RUST
/// Arguments:
/// - genes: vector of gene_names
/// - metrics: vector of ranking values
/// - gene_sets: a hashmap (dict) of GMT file
/// - weight
/// - min_size
/// - max_size
/// - nperm: number of permutation
/// - threads: number of threads
/// - seed: random seed
#[pyfunction]
fn prerank_rs(
    py: Python<'_>,
    genes: Vec<String>,
    metric: Vec<f64>,
    gene_sets: HashMap<String, Vec<String>>,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: usize,
    threads: usize,
    seed: u64,
) -> PyResult<GSEAResult> {
    let pool = build_pool(threads)?;
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    // Release the GIL: the compute is pure Rust and does not touch Python objects,
    // so other Python threads can run while we crunch numbers.
    py.detach(|| pool.install(|| gsea.prerank(&genes, &metric, &gmt)));
    Ok(gsea)
}

/// Prerank RUST
/// Arguments:
/// - genes: vector of gene_names
/// - metrics: 2d vector input with shape [N_gene, N_samples]. handy for multiple ranking metrics input
/// - gene_sets: a hashmap (dict) of GMT file
/// - weight
/// - min_size
/// - max_size
/// - nperm: number of permutation
/// - threads: number of threads
/// - seed: random seed
#[pyfunction]
fn prerank2d_rs(
    py: Python<'_>,
    genes: Vec<String>,
    metric: Vec<Vec<f64>>,
    gene_sets: HashMap<String, Vec<String>>,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: usize,
    threads: usize,
    seed: u64,
) -> PyResult<GSEAResult> {
    let pool = build_pool(threads)?;
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    py.detach(|| pool.install(|| gsea.prerank2(&genes, &metric, &gmt)));
    Ok(gsea)
}

/// GSEA RUST
/// Arguments:
/// - gene_name: vector of gene_names
/// - gene_exp: gene_expression table. each row is gene, each column is sample
/// - gene_sets: a hashmap (dict) of GMT file
/// - group: bool vector of the sample group
/// - weight
/// - min_size
/// - max_size
/// - nperm: number of permutation
/// - threads: number of threads
/// - seed: random seed
#[pyfunction]
fn gsea_rs(
    py: Python<'_>,
    gene_name: Vec<String>,
    gene_exp: Vec<Vec<f64>>,
    gene_sets: HashMap<String, Vec<String>>,
    group: Vec<bool>,
    method: Metric,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: usize,
    threads: usize,
    seed: u64,
) -> PyResult<GSEAResult> {
    let pool = build_pool(threads)?;

    // get gene sets dict
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }

    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    py.detach(|| pool.install(|| gsea.gsea(&gene_name, &group, &gene_exp, &gmt, method)));
    Ok(gsea)
}

/// ssGSEA RUST
/// Arguments:
/// - gene_name: vector of gene_names
/// - gene_exp: gene_expression table. each row is sample, each column is gene values
/// - gene_sets: a hashmap (dict) of GMT file
/// - sample_names: vector of sample names
/// - weight
/// - min_size
/// - max_size
/// - nperm: number of permutation
/// - threads: number of threads
/// - seed: random seed
/// #[pyo3(signature)]
#[pyfunction]
#[pyo3(signature = (gene_name, gene_exp, gene_sets, weight = 1.0, min_size = 5, max_size = 500, nperm = None, ctype = CorrelType::Rank, threads = 4,seed = 0))]
fn ssgsea_rs(
    py: Python<'_>,
    gene_name: Vec<String>,
    gene_exp: Vec<Vec<f64>>,
    gene_sets: HashMap<String, Vec<String>>,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: Option<usize>,
    ctype: CorrelType,
    threads: usize,
    seed: u64,
) -> PyResult<GSEAResult> {
    let pool = build_pool(threads)?;

    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let _nperm = nperm.unwrap_or(0);
    let mut gsea = GSEAResult::new(weight, max_size, min_size, _nperm, seed);
    py.detach(|| {
        pool.install(|| {
            if _nperm > 0 {
                gsea.ss_gsea_permuate(&gene_name, &gene_exp, &gmt, ctype);
            } else {
                gsea.ss_gsea(&gene_name, &gene_exp, &gmt, ctype);
            }
        })
    });
    Ok(gsea)
}
/// Run GSVA (Gene Set Variation Analysis)
///
/// Compute gene set variation scores from the input gene expression data
///
/// Parameters
/// ----------
/// gene_name: list-like of str
///     Gene names
/// gene_expr: array-like of float
///     Gene expression data
/// gene_sets: dict-like of str : list-like of str
///     Gene sets
/// kcdf: bool
///     Use kernel density estimation
/// rnaseq: bool
///     Use RNA-seq data
/// mx_diff: bool
///     Use maximum difference between genes
/// abs_rnk: bool
///     Use absolute rank
/// tau: float
///     Weight for pathway score
/// min_size: int
///     Minimum number of genes in a gene set
/// max_size: int
///     Maximum number of genes in a gene set
/// threads: int
///     Number of threads
///
/// Returns
/// -------
/// GSEAResult
///     A GSEAResult object containing the GSVA results
#[pyfunction]
fn gsva_rs(
    py: Python<'_>,
    gene_name: Vec<String>,
    gene_expr: Vec<Vec<f64>>,
    gene_sets: HashMap<String, Vec<String>>,
    kcdf: bool,
    rnaseq: bool,
    mx_diff: bool,
    abs_rnk: bool,
    tau: f64,
    min_size: usize,
    max_size: usize,
    threads: usize
) -> PyResult<GSEAResult> {
    let pool = build_pool(threads)?;
    let gs = py.detach(|| {
        pool.install(|| {
            gsva(
                gene_name, gene_expr, gene_sets, kcdf, rnaseq, mx_diff, abs_rnk, tau, min_size,
                max_size,
            )
        })
    });
    Ok(gs)
}

/// Preranked GSEA using the fgsea multilevel p-value algorithm.
///
/// Ports `fgseaMultilevel` from the fgsea R/C++ package. Two-phase computation:
///
/// 1. **NES**: `n_perm_simple` simple gene permutations build a per-gene-set null ES
///    distribution. `NES = ES / mean(positive null ESs)` for ES ≥ 0, or
///    `ES / |mean(negative null ESs)|` for ES < 0 — identical to fgsea.
/// 2. **p-value**: adaptive multilevel splitting + MCMC gives arbitrarily precise
///    values with quantified error bound `log2err`.
/// 3. **FDR**: Benjamini-Hochberg correction on the multilevel p-values.
///
/// Arguments:
/// - genes: gene names in rank-descending order
/// - metric: gene-level ranking metric (same order as genes)
/// - gene_sets: a hashmap (dict) of GMT gene sets
/// - weight: GSEA weight parameter (default 1.0)
/// - min_size: minimum gene set size to test
/// - max_size: maximum gene set size to test
/// - sample_size: MCMC sample size per level (fgsea default: 101)
/// - nperm: gene permutations for NES normalization (fgsea default: 1000;
///   set to 0 to skip NES normalization)
/// - eps: convergence threshold for multilevel algorithm (default 1e-50)
/// - threads: number of parallel threads
/// - seed: random seed
#[pyfunction]
#[pyo3(signature = (genes, metric, gene_sets, weight=1.0, min_size=15, max_size=500, sample_size=101, nperm = 1000, eps=1e-50, threads=4, seed=0))]
fn fgsea_rs(
    py: Python<'_>,
    genes: Vec<String>,
    metric: Vec<f64>,
    gene_sets: HashMap<String, Vec<String>>,
    weight: f64,
    min_size: usize,
    max_size: usize,
    sample_size: usize,
    nperm: usize,
    eps: f64,
    threads: usize,
    seed: u64,
) -> PyResult<GSEAResult> {
    let pool = build_pool(threads)?;
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    py.detach(|| {
        pool.install(|| gsea.prerank_multilevel(&genes, &metric, &gmt, sample_size, eps))
    });
    Ok(gsea)
}


/// Compute the GSEA enrichment score and optional curve data for a single gene set.
///
/// This is a direct port of fgsea's `calcGseaStat()` — an O(k) algorithm where k is
/// the gene set size.  It is useful both for computing the scalar ES and for generating
/// enrichment-plot curves (`tops` / `bottoms`) without iterating over all N genes.
///
/// Parameters
/// ----------
/// stats : list of float
///     Gene-level ranking metric, **sorted descending**, length N.
///     Values are raised to the power `gsea_param` internally.
/// selected_stats : list of int
///     **0-based** indices of gene-set members in `stats`, **sorted ascending**.
/// gsea_param : float, optional
///     Weight exponent p (default 1.0).  Use 0 for unweighted ES.
/// score_type : ScoreType, optional
///     Which extreme to return: ``Std`` (default, sign-aware max |ES|),
///     ``Pos`` (maximum), or ``Neg`` (minimum).
/// return_all_extremes : bool, optional
///     If True, populate ``tops`` and ``bottoms`` in the result (needed for
///     enrichment-plot curves).  Default False.
/// return_leading_edge : bool, optional
///     If True, populate ``leading_edge`` in the result.  Default False.
///
/// Returns
/// -------
/// GseaStatResult
///     Object with fields ``es``, ``tops``, ``bottoms``, ``leading_edge``.
// #[pyfunction]
// #[pyo3(signature = (stats, selected_stats, gsea_param=1.0, score_type=ScoreType::Std, return_all_extremes=false, return_leading_edge=false))]
// fn calc_gsea_stat_rs(
//     stats: Vec<f64>,
//     selected_stats: Vec<usize>,
//     gsea_param: f64,
//     score_type: ScoreType,
//     return_all_extremes: bool,
//     return_leading_edge: bool,
// ) -> PyResult<GseaStatResult> {
//     Ok(calc_gsea_stat(
//         &stats,
//         &selected_stats,
//         gsea_param,
//         score_type,
//         return_all_extremes,
//         return_leading_edge,
//     ))
// }

/// Python module for GSEA (Gene Set Enrichment Analysis) and ssGSEA
/// 
/// This module provides five functions:
/// 
/// - `gsea_rs`: performs GSEA
/// - `prerank_rs`: performs GSEA preranking
/// - `prerank2d_rs`: performs GSEA preranking with multiple input datasets
/// - `ssgsea_rs`: performs ssGSEA
/// - `gsva_rs`: performs GSVA
/// - `prerank_fgsea_rs`: performs preranked GSEA with the fgsea multilevel p-value algorithm
/// - `calc_gsea_stat_rs`: O(k) enrichment score + optional curve data for a single gene set
#[pymodule]
fn gse(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // m.add_class::<GSEAResult>()?;
    m.add_class::<GSEASummary>()?;
    m.add_class::<GSEAResult>()?;
    m.add_class::<Metric>()?;
    m.add_class::<CorrelType>()?;
    m.add_class::<ScoreType>()?;
    m.add_class::<GseaStatResult>()?;
    m.add_function(wrap_pyfunction!(gsea_rs, m)?)?;
    m.add_function(wrap_pyfunction!(prerank_rs, m)?)?;
    m.add_function(wrap_pyfunction!(prerank2d_rs, m)?)?;
    m.add_function(wrap_pyfunction!(ssgsea_rs, m)?)?;
    m.add_function(wrap_pyfunction!(gsva_rs, m)?)?;
    m.add_function(wrap_pyfunction!(fgsea_rs, m)?)?;
    // m.add_function(wrap_pyfunction!(calc_gsea_stat_rs, m)?)?;
    Ok(())
}

#[cfg(all(test, feature = "extension-module"))]
mod thread_pool_tests {
    //! Tests for the per-call thread-pool fix (replaces the old
    //! `env::set_var("RAYON_NUM_THREADS", ..)` global-state hack).
    use super::build_pool;
    use rayon::prelude::*;

    /// Every call gets a pool with *exactly* the requested number of threads.
    #[test]
    fn pool_honors_requested_thread_count() {
        for n in [1usize, 2, 3] {
            let pool = build_pool(n).expect("pool should build");
            let observed = pool.install(rayon::current_num_threads);
            assert_eq!(observed, n, "pool should expose exactly {n} threads");
        }
    }

    /// The core regression: two differently-sized pools coexist in one process.
    /// This is impossible with `build_global()` (callable once) or
    /// `RAYON_NUM_THREADS` (read once at lazy global-pool init) — the exact bug
    /// that made `threads=` a no-op after the first call in a Python session.
    #[test]
    fn pools_are_independent_per_call() {
        let small = build_pool(1).expect("pool should build");
        let big = build_pool(4).expect("pool should build");
        assert_eq!(small.install(rayon::current_num_threads), 1);
        assert_eq!(big.install(rayon::current_num_threads), 4);
        // ...and the first pool is unchanged by the second's creation.
        assert_eq!(small.install(rayon::current_num_threads), 1);
    }

    /// Nested `par_iter` inside `install` inherits the installed pool — this is
    /// why wrapping only the top-level compute call is sufficient, even though
    /// the `par_iter` calls live deep in stats.rs / algorithm.rs / gsva.rs.
    #[test]
    fn nested_parallel_iter_uses_installed_pool() {
        let pool = build_pool(2).expect("pool should build");
        let observed = pool.install(|| {
            (0..1_000)
                .into_par_iter()
                .map(|_| rayon::current_num_threads())
                .max()
                .unwrap()
        });
        assert_eq!(observed, 2, "nested par_iter should run on the 2-thread pool");
    }
}
