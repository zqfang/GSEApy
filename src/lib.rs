use pyo3::prelude::*;
use std::collections::HashMap;
use std::env;
// import own modules
mod algorithm;
mod gsva;
mod stats;
mod utils;
// export module fn, struct, trait ...
use gsva::gsva;
use stats::{GSEAResult, GSEASummary};
use utils::{CorrelType, Metric};

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
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    env::set_var("RAYON_NUM_THREADS", threads.to_string());
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    gsea.prerank(&genes, &metric, &gmt);
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
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    env::set_var("RAYON_NUM_THREADS", threads.to_string());
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    gsea.prerank2(&genes, &metric, &gmt);
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
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    // set number of threads
    env::set_var("RAYON_NUM_THREADS", threads.to_string());

    // get gene sets dict
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }

    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    gsea.gsea(&gene_name, &group, &gene_exp, &gmt, method);
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
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    env::set_var("RAYON_NUM_THREADS", threads.to_string());

    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter() {
        gmt.insert(k.as_str(), v.as_slice());
    }
    let _nperm = nperm.unwrap_or(0);
    let mut gsea = GSEAResult::new(weight, max_size, min_size, _nperm, seed);
    if _nperm > 0 {
        gsea.ss_gsea_permuate(&gene_name, &gene_exp, &gmt, ctype);
    } else {
        gsea.ss_gsea(&gene_name, &gene_exp, &gmt, ctype);
    }
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
    env::set_var("RAYON_NUM_THREADS", threads.to_string());
    let gs = gsva(
        gene_name, gene_expr, gene_sets, kcdf, rnaseq, mx_diff, abs_rnk, tau, min_size, max_size
    );
    Ok(gs)
}

/// Python module for GSEA (Gene Set Enrichment Analysis) and ssGSEA
/// 
/// This module provides four functions:
/// 
/// - `gsea_rs`: performs GSEA
/// - `prerank_rs`: performs GSEA preranking
/// - `prerank2d_rs`: performs GSEA preranking with multiple input datasets
/// - `ssgsea_rs`: performs ssGSEA
/// - `gsva_rs`: performs GSVA
#[pymodule]
fn gse(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // m.add_class::<GSEAResult>()?;
    m.add_class::<GSEASummary>()?;
    m.add_class::<GSEAResult>()?;
    m.add_class::<Metric>()?;
    m.add_class::<CorrelType>()?;
    m.add_function(wrap_pyfunction!(gsea_rs, m)?)?;
    m.add_function(wrap_pyfunction!(prerank_rs, m)?)?;
    m.add_function(wrap_pyfunction!(prerank2d_rs, m)?)?;
    m.add_function(wrap_pyfunction!(ssgsea_rs, m)?)?;
    m.add_function(wrap_pyfunction!(gsva_rs, m)?)?;
    Ok(())
}
