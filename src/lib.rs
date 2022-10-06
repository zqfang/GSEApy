
use std::env;
use std::collections::HashMap;
use pyo3::prelude::*;
// import own modules
mod algorithm;
mod stats;
mod utils;
// export module fn, struct, trait ...
use stats::{GSEAResult, GSEASummary};
use utils::{Metric};

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
    seed: u64
) -> PyResult<Vec<GSEASummary>> {
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    env::set_var("RAYON_NUM_THREADS", threads.to_string());
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter(){
        gmt.insert(k.as_str(), v.as_slice());
    }
    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    gsea.prerank(&genes, &metric, &gmt);
    Ok(gsea.summaries)
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
) -> PyResult<Vec<GSEASummary>> {
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    // set number of threads
    env::set_var("RAYON_NUM_THREADS", threads.to_string());

    // get gene sets dict
    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter(){
        gmt.insert(k.as_str(), v.as_slice());
    }

    let mut gsea = GSEAResult::new(weight, max_size, min_size, nperm, seed);
    gsea.gsea(&gene_name, &group, &gene_exp, &gmt, method);
    Ok(gsea.summaries)
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
#[pyfunction]
fn ssgsea_rs(
    gene_name: Vec<String>,
    gene_exp: Vec<Vec<f64>>,
    gene_sets: HashMap<String, Vec<String>>,
    sample_names: Vec<String>,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: Option<usize>,
    threads: usize,
    seed: u64,
) -> PyResult<Vec<GSEASummary>> {
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .unwrap_or_default();
    env::set_var("RAYON_NUM_THREADS", threads.to_string());

    let mut gmt = HashMap::<&str, &[String]>::new();
    for (k, v) in gene_sets.iter(){
        gmt.insert(k.as_str(), v.as_slice());
    }
    let _nperm = nperm.unwrap_or(0);
    let mut gsea = GSEAResult::new(weight, max_size, min_size, _nperm, seed);
    if _nperm > 0 {
        gsea.ss_gsea_permuate(&gene_name, &sample_names, &gene_exp, &gmt);
    } 
    else 
    {
        gsea.ss_gsea(&gene_name, &sample_names, &gene_exp, &gmt);
    }
    Ok(gsea.summaries)
}



#[pymodule]
fn gse(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // m.add_class::<GSEAResult>()?;
    m.add_class::<GSEASummary>()?;
    m.add_class::<Metric>()?;
    m.add_function(wrap_pyfunction!(gsea_rs, m)?)?;
    m.add_function(wrap_pyfunction!(prerank_rs, m)?)?;
    m.add_function(wrap_pyfunction!(ssgsea_rs, m)?)?;
    Ok(())
}
