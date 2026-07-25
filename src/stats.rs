#![allow(dead_code, unused)]

use crate::algorithm::{EnrichmentScore, EnrichmentScoreTrait};
use crate::fgsea::util::multilevelError;
use crate::fgsea::{calcGseaStatCumulativeBatch, compute_pvalue_multilevel, scale_ranks};
use crate::utils::{CorrelType, DynamicEnum, Metric, Statistic};
use itertools::{izip, Itertools};
use pyo3::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use special::Gamma;
use statrs::function::beta::beta_reg;
use std::collections::{BTreeMap, HashMap};

/// The `p`-quantile of a Beta(`a`, `b`) distribution -- R's `qbeta(p, a, b)`.
///
/// Found by bisection on the regularised incomplete beta function (Beta's CDF),
/// which is monotone on `[0, 1]`. 200 iterations take the bracket well below f64
/// resolution. Degenerate shapes match R's boundary conventions: `qbeta` is 0 when
/// `a <= 0` and 1 when `b <= 0`.
fn qbeta(p: f64, a: f64, b: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    if b <= 0.0 {
        return 1.0;
    }
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if beta_reg(a, b, mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Half-width (in log2 units) of the 95% interval around the simple permutation
/// p-value estimate `(n_more_extreme + 1) / (n_perm + 1)`.
///
/// Mirrors the `leftBorder`/`rightBorder`/`crudeEstimator` block of `fgseaMultilevel`.
/// With `n_more_extreme == 0` the left border is `log2(0) = -inf`, so the error is
/// infinite and the multilevel estimator always wins -- which is the intent: zero
/// more-extreme permutations tells you nothing about how small the p-value is.
fn simple_error_log2(n_more_extreme: f64, n_perm: f64) -> f64 {
    let left_border = qbeta(0.025, n_more_extreme, n_perm - n_more_extreme + 1.0).log2();
    let right_border = qbeta(0.975, n_more_extreme + 1.0, n_perm - n_more_extreme).log2();
    let crude = ((n_more_extreme + 1.0) / (n_perm + 1.0)).log2();
    0.5 * (crude - left_border).max(right_border - crude)
}

/// Expected log2 error of the multilevel estimate of `pval` -- fgsea's
/// `multilevelError(pval, sampleSize)`.
///
/// The multilevel algorithm needs `floor(-log2(pval) + 1)` splitting levels to reach
/// `pval`, and each level contributes independent sampling error.
fn multilevel_error_pval(pval: f64, sample_size: usize) -> f64 {
    let level = (-pval.log2() + 1.0).floor();
    multilevelError(level as i32, sample_size as i32)
}

#[pyclass(from_py_object)]
#[allow(dead_code)]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GSEASummary {
    #[pyo3(get, set)]
    pub term: String,
    #[pyo3(get, set)]
    pub es: f64,
    #[pyo3(get, set)]
    pub nes: f64,
    #[pyo3(get, set)]
    pub pval: f64, // Nominal Pvalue
    #[pyo3(get, set)]
    pub fwerp: f64, // FWER Pvalue
    #[pyo3(get, set)]
    pub fdr: f64, // FDR q value. adjusted FDR
    #[pyo3(get, set)]
    pub log2err: f64, // log2 error bound for multilevel p-value (NaN when not applicable)
    #[pyo3(get, set)]
    pub run_es: Vec<f64>,
    #[pyo3(get, set)]
    pub hits: Vec<usize>, // indices of genes that matches
    #[pyo3(get, set)]
    pub esnull: Vec<f64>,
    #[pyo3(get, set)]
    pub index: Option<usize>, // sample index
}

impl GSEASummary {
    pub fn new(
        &mut self,
        term: &str,
        es: f64,
        nes: f64,
        pval: f64,
        fwerpval: f64,
        fdr: f64,
        run_es: &[f64],
        hits: &[usize],
        esnull: &[f64],
        index: usize,
    ) -> Self {
        GSEASummary {
            term: term.to_string(),
            es: es,
            nes: nes,
            pval: pval,
            fwerp: fwerpval,
            fdr: fdr,
            log2err: f64::NAN,
            run_es: run_es.to_vec(),
            hits: hits.to_vec(),
            esnull: esnull.to_vec(),
            index: Some(index),
        }
    }

    /// for default values, you can then init the struct with
    /// let g = GSEASummary { es: 0.5, ..Default::default() };
    /// need trait bound #[derive(Default)]
    #[allow(dead_code)]
    fn default() -> GSEASummary {
        GSEASummary {
            term: "".to_string(),
            es: 0.0,
            nes: 0.0,
            pval: 1.0,
            fwerp: 1.0,
            fdr: 1.0,
            log2err: f64::NAN,
            run_es: Vec::<f64>::new(),
            hits: Vec::<usize>::new(),
            esnull: Vec::<f64>::new(),
            index: None,
        }
    }
    /// see the normalizatin code from
    /// https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    fn normalize(&mut self) -> Vec<f64> {
        let e: f64 = self.es;
        // n_mean = esnull[esnull>= 0].mean()
        let pos_phi: Vec<f64> = self
            .esnull
            .iter()
            .filter_map(|&x| if x >= 0.0 { Some(x) } else { None })
            .collect();

        // n_mean = esnull[esnull< 0].mean()
        let neg_phi: Vec<f64> = self
            .esnull
            .iter()
            .filter_map(|&x| if x < 0.0 { Some(x) } else { None })
            .collect();

        // FIXME: Potential NaN number here
        // When input a rare causes of an extreamly screwed null distribution. e.g.
        // es = - 27, esnull = [13, 24, 57, 88]
        // nes will be NaN. You have to increased the permutation number for safe
        // a tricky fixed here: set n_mean as itself
        // so esnull = [-27, 13, 24, 57, 88]
        let pos_mean = if pos_phi.len() > 0 {
            pos_phi.as_slice().mean()
        } else {
            e
        };

        let neg_mean = if neg_phi.len() > 0 {
            neg_phi.as_slice().mean()
        } else {
            e
        };

        self.nes = if e >= 0.0 {
            e / pos_mean
        } else {
            e / neg_mean.abs()
        };

        let nesnull: Vec<f64> = self
            .esnull
            .iter()
            .map(|&e| {
                if e >= 0.0 {
                    e / pos_mean
                } else {
                    e / neg_mean.abs()
                }
            })
            .collect();
        // store normalized esnull temporatory.
        nesnull
    }

    fn pval(&mut self) {
        let deno: usize;
        let nomi: usize;
        // When input a rare causes of an extreamly screwed null distribution. e.g.
        // es = - 27, esnull = [13, 24, 57, 88]
        // pval will be NaN.
        if self.es < 0.0 {
            deno = self.esnull.iter().filter(|&x| *x < 0.0).count();
            nomi = self.esnull.iter().filter(|&x| x <= &self.es).count();
        } else {
            deno = self.esnull.iter().filter(|&x| *x >= 0.0).count();
            nomi = self.esnull.iter().filter(|&x| x >= &self.es).count();
        }

        if deno == 0 {
            self.pval = 1.0;
            return;
        }
        self.pval = (nomi as f64) / (deno as f64);
    }
}

#[pyclass(skip_from_py_object)]
#[allow(dead_code)]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GSEAResult {
    #[pyo3(get, set)]
    pub summaries: Vec<GSEASummary>,
    #[pyo3(get, set)]
    pub weight: f64,
    #[pyo3(get, set)]
    pub min_size: usize,
    #[pyo3(get, set)]
    pub max_size: usize,
    #[pyo3(get, set)]
    pub nperm: usize,
    nes_concat: Vec<f64>,
    nesnull_concat: Vec<f64>,
    #[pyo3(get, set)]
    pub seed: u64,
    #[pyo3(get, set)]
    pub rankings: Vec<Vec<f64>>,
    #[pyo3(get, set)]
    pub indices: Vec<Vec<usize>>, // indices after ranking
}

impl GSEAResult {
    pub fn new(weight: f64, max_size: usize, min_size: usize, nperm: usize, seed: u64) -> Self {
        GSEAResult {
            summaries: Vec::<GSEASummary>::new(),
            weight: weight,
            max_size: max_size,
            min_size: min_size,
            nperm: nperm,
            nes_concat: Vec::<f64>::new(),
            nesnull_concat: Vec::<f64>::new(),
            seed: seed,
            rankings: Vec::<Vec<f64>>::new(),
            indices: Vec::<Vec<usize>>::new(),
        }
    }
    pub fn default() -> GSEAResult {
        GSEAResult {
            summaries: Vec::<GSEASummary>::new(),
            weight: 1.0,
            max_size: 1000,
            min_size: 3,
            nperm: 1000,
            nes_concat: Vec::<f64>::new(),
            nesnull_concat: Vec::<f64>::new(),
            seed: 0,
            rankings: Vec::<Vec<f64>>::new(),
            indices: Vec::<Vec<usize>>::new(),
        }
    }
    pub fn stat(&mut self, summary: &mut [GSEASummary]) {
        // Group summaries by library prefix to compute FDR/FWER per library.
        // When multiple gene set libraries are combined (terms prefixed as "library__term"),
        // FDR should be computed separately for each library so that results for a given
        // term are not affected by how many other libraries are included in the run.
        let mut lib_groups: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, g) in summary.iter().enumerate() {
            let prefix = match g.term.find("__") {
                Some(idx) => g.term[..idx].to_string(),
                None => String::new(),
            };
            lib_groups.entry(prefix).or_default().push(i);
        }

        for (_lib, indices) in &lib_groups {
            // clear vector incase you re-run this command
            self.nes_concat.clear();
            self.nesnull_concat.clear();

            for &i in indices {
                let g = &mut summary[i];
                // calculate stats here
                g.pval();
                let mut nesnull = g.normalize(); // update esnull to normalized nesnull
                self.nes_concat.push(g.nes);
                self.nesnull_concat.append(&mut nesnull);
                // g.esnull.clear();
            }
            // FWER p
            let fwerps: Vec<f64> = self.fwer_pval();
            // FDR q
            let fdrs = self.fdr();

            for (j, &i) in indices.iter().enumerate() {
                summary[i].fdr = fdrs[j];
                summary[i].fwerp = fwerps[j];
            }
        }
        // clear vector to save some space
        self.nes_concat.clear();
        self.nesnull_concat.clear();
    }

    /// see line 844 - 876: https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    pub fn _fdr(&mut self) -> Vec<f64> {
        let nes_idx = self.nes_concat.iter().filter(|&x| *x < 0.0).count();
        // let mut nesnull_concat: Vec<&f64> = nesnull.iter().flatten().collect(); // nesnull.concat(); // concat items
        let fdrs: Vec<f64> = self
            .nes_concat
            .iter()
            .enumerate()
            .map(|(i, &e)| {
                let mut phi_norm: f64;
                let mut phi_obs: f64;
                let mut nes_higher: usize;
                let mut all_higher: usize;
                let mut all_pos: usize;
                let mut nes_pos: usize;
                let mut fdrs_all: Vec<f64> = Vec::new();
                for j in i..self.nperm {
                    let indexes = (j..self.nesnull_concat.len())
                        .step_by(self.nperm)
                        .into_iter();
                    let nesnull: Vec<f64> = indexes.map(|m| self.nesnull_concat[m]).collect();

                    if e < 0.0 {
                        nes_higher = self.nes_concat.iter().filter(|&x| *x <= e).count();
                        all_higher = nesnull.iter().filter(|&x| *x <= e).count();
                        all_pos = nesnull.iter().filter(|&x| *x < 0.0).count();
                        nes_pos = nes_idx;
                    } else {
                        nes_higher = self.nes_concat.iter().filter(|&x| *x >= e).count();
                        all_higher = nesnull.iter().filter(|&x| *x >= e).count();
                        all_pos = nesnull.iter().filter(|&x| *x >= 0.0).count();
                        nes_pos = self.nes_concat.len() - nes_idx;
                    }
                    // println!("neg_higher {}, all_higher {}, all_pos {}, nes_pos {}", nes_higher, all_higher, all_pos, all_higher);
                    phi_norm = if all_pos > 0 {
                        (all_higher as f64) / (all_pos as f64)
                    } else {
                        0.0
                    }; // count.col
                    phi_obs = if nes_pos > 0 {
                        (nes_higher as f64) / (nes_pos as f64)
                    } else {
                        0.0
                    }; // obs.count.col
                       // FDR is undefined when phi_obs == 0 (0/0 -> NaN); by
                       // convention set it to the maximum (1.0). Otherwise clamp
                       // to the valid [0, 1] range.
                    let fdr = if phi_obs == 0.0 {
                        1.0
                    } else {
                        (phi_norm / phi_obs).clamp(0.0, 1.0)
                    };
                    fdrs_all.push(fdr);
                }
                fdrs_all.as_slice().mean()
            })
            .collect();
        return fdrs;
    }

    /// see line 844 - 876: https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    /// To speed up the FDR computation, I used an expectation approximate to estimate the FDR.mean in the R code.
    pub fn fdr(&mut self) -> Vec<f64> {
        // let mut nesnull_concat: Vec<&f64> = nesnull.iter().flatten().collect(); // nesnull.concat(); // concat items

        // To speedup, sort f64 in acending order in place, then do a binary search
        self.nesnull_concat
            .sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)); // NaN-safe; if descending -> b.partial_cmp(a)
        let (indices, nes_sorted) = self.nes_concat.as_slice().argsort(true); // ascending order

        // binary_search assumes that the elements are sorted in less-to-greater order.
        // partition_point return the index of the first element of the second partition)
        // since partition_point is just a wrapper of self.binary_search_by(|x| if pred(x) { Less } else { Greater }).unwrap_or_else(|i| i)
        let all_idx = self.nesnull_concat.partition_point(|x| *x < 0.0);
        let nes_idx = nes_sorted.partition_point(|x| *x < 0.0);

        // fdr
        let fdrs: Vec<f64> = nes_sorted
            .iter()
            .map(|&e| {
                let phi_norm: f64;
                let phi_obs: f64;
                let nes_higher: usize;
                let all_higher: usize;
                let all_pos: usize;
                let nes_pos: usize;
                if e < 0.0 {
                    // let nes_higher = nes_concat.iter().filter(|&x| *x < e).count();
                    // let all_higher = nesnull_concat.iter().filter(|&x| *x < e).count();
                    nes_higher = nes_sorted.partition_point(|x| *x <= e); // left side
                    all_higher = self.nesnull_concat.partition_point(|x| *x <= e); // left side
                    all_pos = all_idx;
                    nes_pos = nes_idx;
                } else {
                    // let nes_higher = self.nes_concat.iter().filter(|&x| *x >= e).count();
                    // let all_higher = self.nesnull_concat.iter().filter(|&x| *x >= e).count();
                    nes_higher = nes_sorted.len() - nes_sorted.partition_point(|x| *x < e); // right side
                    all_higher =
                        self.nesnull_concat.len() - self.nesnull_concat.partition_point(|x| *x < e); // right side; count.col ( /count.col.norm)
                    all_pos = self.nesnull_concat.len() - all_idx; // right side; count.col.norm
                    nes_pos = nes_sorted.len() - nes_idx; // right side; obs.count.col.norm
                }
                // println!("neg_higher {}, all_higher {}, all_pos {}, nes_pos {}", nes_higher, all_higher, all_pos, all_higher);
                phi_norm = if all_pos > 0 {
                    (all_higher as f64) / (all_pos as f64)
                } else {
                    0.0
                }; // count.col
                phi_obs = if nes_pos > 0 {
                    (nes_higher as f64) / (nes_pos as f64)
                } else {
                    0.0
                }; // obs.count.col
                   // FDR is undefined when phi_obs == 0 (0/0 -> NaN); by
                   // convention set it to the maximum (1.0), else clamp to [0, 1].
                if phi_obs == 0.0 {
                    1.0
                } else {
                    (phi_norm / phi_obs).clamp(0.0, 1.0)
                }
            })
            .collect();

        // by default, we'er no gnna adjusted fdr q value
        // self.adjust_fdr(&mut fdrs, nes_idx);
        let mut fdr_orig_order: Vec<f64> = vec![0.0; fdrs.len()];
        indices.iter().zip(fdrs.iter()).for_each(|(&i, &v)| {
            fdr_orig_order[i] = v;
        });
        return fdr_orig_order;
    }

    /// # adjust fdr q-values
    /// see line 880:  https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    /// - fdrs:  Corresponds to the ascending order of NES.
    /// - partition_point_idx: the index of the first element of the second partition
    /// This function updates fdr value inplace.
    #[allow(dead_code)]
    fn adjust_fdr(&self, fdrs: &mut [f64], partition_point_idx: usize) {
        // If NES is a so screwd distribution, e.g. all positive or negative numbers.
        // partition_point_idx will be either of 0 or fdrs.len(). Need to skip. example here:
        // let s1 = [1,3,4,5,6,9];
        // let s2 = [-10, -8, -7,-4,-1];
        // let s3 = [-9,-8,-2,-1,1,2,3];

        // let b1 = s1.partition_point(|x| *x < 0);
        // let b2 = s2.partition_point(|x| *x < 0); neg_nes on the left
        // let b3 = s3.partition_point(|x| *x < 0);
        // the partition_point_idx will be: b1 = 0, b2 = 5, b3 = 4

        // thus, the transver order is opposit to the R code since we'er using acsending order of nes
        let mut min_fdr: f64;
        if partition_point_idx < fdrs.len() {
            // pos_nes on the right side, if only have postive numbers, idx must be < .len()
            let nes_pos_idx = partition_point_idx + 1;
            min_fdr = fdrs[partition_point_idx];
            for k in nes_pos_idx..fdrs.len() {
                // if fdrs[k] < min_fdr {
                //     min_fdr = fdrs[k]
                // }
                // if min_fdr < fdrs[k] {
                //     fdrs[k] = min_fdr
                // }
                min_fdr = min_fdr.min(fdrs[k]);
                fdrs[k] = min_fdr.min(fdrs[k]);
            }
        }

        if partition_point_idx > 0 {
            // neg_nes on the left side, if only have negative numbers, idx must be > 0
            let nes_neg_idx = partition_point_idx - 1;
            min_fdr = fdrs[nes_neg_idx];
            for k in (0..partition_point_idx).rev() {
                min_fdr = min_fdr.min(fdrs[k]);
                fdrs[k] = min_fdr.min(fdrs[k]);
            }
        }
    }
    /// Compute FWER p-vals
    /// line 788: https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    fn fwer_pval(&self) -> Vec<f64> {
        // suppose a matrix of nesnull with shape [ n_genesets, n_perm ]
        // max_nes_pos = colMax(nesull) for nes >= 0;
        // min_nes_neg = colMin(nesnull) for nes < 0;
        let mut max_nes_pos = vec![0.0; self.nperm];
        let mut min_nes_neg = vec![0.0; self.nperm];
        self.nesnull_concat.iter().enumerate().for_each(|(i, &e)| {
            let idx = i % self.nperm;
            if e >= 0.0 {
                max_nes_pos[idx] = e.max(max_nes_pos[idx]);
            } else {
                min_nes_neg[idx] = e.min(min_nes_neg[idx]);
            }
        });

        let fwerp: Vec<f64> = self
            .nes_concat
            .par_iter()
            .map(|e| {
                if e < &0.0 {
                    (min_nes_neg.iter().filter(|&x| x < e).count() as f64)
                        / (min_nes_neg.iter().filter(|&x| x < &0.0).count() as f64)
                } else {
                    (max_nes_pos.iter().filter(|&x| x >= e).count() as f64)
                        / (max_nes_pos.len() as f64)
                }
            })
            .collect();
        fwerp
    }
}

/// impl pipelines
impl GSEAResult {
    pub fn gsea(
        &mut self,
        genes: &[String],
        group: &[bool],
        gene_exp: &[Vec<f64>],
        gmt: &BTreeMap<&str, &[String]>,
        method: Metric,
    ) {
        let mut es = EnrichmentScore::new(genes, self.nperm, self.seed, false, false);
        // let end = Instant::now();
        let mut sorted_metric: Vec<(Vec<usize>, Vec<f64>)> =
            es.phenotype_permutation(gene_exp, group, method, false);
        // save indices and rankings
        self.rankings.push(sorted_metric[0].1.to_owned());
        self.indices.push(sorted_metric[0].0.to_owned());
        // weight 
        sorted_metric.iter_mut().for_each(|x| {
            x.1.iter_mut().for_each(|x| {
                *x = x.abs().powf(self.weight);
            });
        });

        // let end1 = Instant::now();
        // println!("Permutation time: {:.2?}", end1.duration_since(end));

        let mut summ = Vec::<GSEASummary>::new();
        for (&term, &gset) in gmt.iter() {
            let tag = es.gene.isin(gset);
            // get es hit index of sorted array
            let tag_new: Vec<f64> = sorted_metric[0].0.iter().map(|&i| tag[i]).collect();
            let gidx = es.hit_index(&tag_new); // need update the sorted indices
            if gidx.len() > self.max_size || gidx.len() < self.min_size {
                continue;
            }
            // get es
            let run_es = es.running_enrichment_score(sorted_metric[0].1.as_slice(), &tag_new);
            //let es0: f64 = es.select_es(&run_es);
            let ess: Vec<f64> = sorted_metric
                .par_iter()
                .map(|(indices, gm)| {
                    // update tag_indicator since you've update metric
                    let tag_new: Vec<f64> = indices.iter().map(|&i| tag[i]).collect();
                    // calculate ES
                    let r = es.fast_random_walk(gm, &tag_new);
                    r
                })
                .collect();

            // let (ess, run_es) = es.enrichment_score_pheno(&weighted_metric, &tag);
            let esnull: Vec<f64> = if ess.len() > 1 { ess[1..].to_vec() } else { Vec::new() };
            let gss = GSEASummary {
                term: term.to_string(),
                es: ess[0],
                run_es: run_es, // run_es[0].to_owned(),
                hits: gidx,
                esnull: esnull,
                ..Default::default()
            };
            summ.push(gss);
        }

        // let end2 = Instant::now();
        if self.nperm > 0 {
            self.stat(&mut summ);
        }
        self.summaries = summ;
    }

    pub fn prerank(&mut self, genes: &[String], metric: &[f64], gmt: &BTreeMap<&str, &[String]>) {
        // NOTE: input must not contain duplcated genes

        let weighted_metric: Vec<f64> = metric.iter().map(|x| x.abs().powf(self.weight)).collect();
        // start to calculate
        let mut es = EnrichmentScore::new(genes, self.nperm, self.seed, false, false);
        // let end1 = Instant::now();
        let gperm = es.gene_permutation(); // gene permutation, only record gene idx here
        let mut summ = Vec::<GSEASummary>::new();

        for (&term, &gset) in gmt.iter() {
            // convert gene String --> Int
            let gtag = es.gene.isin(gset);
            let gidx = es.hit_index(&gtag);
            if gidx.len() > self.max_size || gidx.len() < self.min_size {
                continue;
            }
            let tag_indicators: Vec<Vec<f64>> = gperm.par_iter().map(|de| de.isin(&gidx)).collect();
            let (ess, run_es) = es.enrichment_score_gene(&weighted_metric, &tag_indicators);
            let esnull: Vec<f64> = if ess.len() > 1 {
                ess[1..].to_vec()
            } else {
                Vec::new()
            };
            let gss = GSEASummary {
                term: term.to_string(),
                es: ess[0],
                run_es: run_es,
                hits: gidx,
                esnull: esnull,
                ..Default::default()
            };
            summ.push(gss);
        }
        // let end3 = Instant::now();
        // println!("Calculation time: {:.2?}", end3.duration_since(end2));
        if self.nperm > 0 {
            self.stat(&mut summ);
        }
        self.summaries = summ;
        self.indices.push((0..genes.len()).collect_vec());
        self.rankings.push(metric.to_owned());
        // let end4 = Instant::now();
        // println!("Statistical time: {:.2?}", end4.duration_since(end3));
    }
    /// multi-preranking datasets input
    /// metric: 2d vector with shape: [N_genes, N_samples]
    pub fn prerank2(
        &mut self,
        genes: &[String],
        metric: &[Vec<f64>], // 2d vector [m_gene, n_sample];
        gmt: &BTreeMap<&str, &[String]>,
    ) {
        // Calculate number of samples
        let n_samples = metric[0].len();
        
        // Pre-allocate vectors to store results
        self.summaries.clear();
        self.indices.clear();
        self.indices.reserve(n_samples);
        
        // Process chunks of samples in parallel to control memory usage
        let chunk_size = (n_samples / rayon::current_num_threads()).max(1);
        
        // Create iterator over sample indices in chunks
        (0..n_samples).chunks(chunk_size).into_iter().for_each(|chunk| {
            let chunk_results: Vec<_> = chunk.collect::<Vec<_>>().into_par_iter().map(|sample_idx| {
                // Create sample vector without transposing entire matrix
                let mut sample_metric: Vec<f64> = Vec::with_capacity(metric.len());
                for row in metric.iter() {
                    sample_metric.push(row[sample_idx]);
                }
                
                // Sort and weight the metrics
                let (indices, sorted_metric) = sample_metric.as_slice().argsort(false);
                
                // Create ordered gene list (reuse allocation)
                let sample_genes: Vec<String> = indices.iter()
                    .map(|&j| genes[j].to_string())
                    .collect();
                
                // Run prerank for this sample
                let mut gsea_tmp = GSEAResult::new(
                    self.weight,
                    self.max_size,
                    self.min_size,
                    self.nperm,
                    self.seed,
                );
                gsea_tmp.prerank(&sample_genes, &sorted_metric, gmt);
                
                // Return only what we need
                (indices, gsea_tmp.summaries.into_iter().map(|mut s| {
                    s.index = Some(sample_idx);
                    s
                }).collect::<Vec<_>>())
            }).collect();
            
            // Process results for this chunk
            for (indices, summaries) in chunk_results {
                self.indices.push(indices);
                self.summaries.extend(summaries);
            }
        });
    }

    pub fn ss_gsea(
        &mut self,
        genes: &[String],
        gene_exp: &[Vec<f64>], // 2d vector [m_gene, n_sample];
        gmt: &BTreeMap<&str, &[String]>,
        correl_type: CorrelType,
    ) {
        // transpose [m_gene, n_sample] --> [n_sample, m_gene]
        let mut gene_metric: Vec<Vec<f64>> = vec![vec![]; gene_exp[0].len()];
        gene_exp.iter().for_each(|row| {
            row.iter().enumerate().for_each(|(j, e)| {
                gene_metric[j].push(*e);
            });
        });

        // sort first and then set weight,
        let weighted_sorted_metric: Vec<(Vec<usize>, Vec<f64>)> = gene_metric
            .into_par_iter()
            .map(|rank| {
                // https://github.com/GSEA-MSigDB/ssGSEA-gpmodule/blob/master/src/ssGSEA.Library.R: line 323
                // in ssGSEA we rank genes by their expression level rather than by a measure of correlation
                // between expression profile and phenotype.
                // transform the normalized expression data for a single sample into ranked (in decreasing order) expression values
                let mut tmp = rank.as_slice().argsort(false);
                // line 338
                if self.weight > 0.0 {
                    // see code
                    // https://github.com/broadinstitute/ssGSEA2.0/blob/f682082f62ae34185421545f25041bae2c78c89b/src/ssGSEA2.0.R#L396
                    match correl_type {
                        CorrelType::SymRank => {
                            let idx = (tmp.1.len() + 2 - 1) / 2; // ceiling div by 2
                            let mid = tmp.1.get(idx).unwrap().to_owned();
                            tmp.1.iter_mut().for_each(|x| {
                                if *x > mid {
                                    *x
                                } else {
                                    (*x) + (*x) - mid
                                };
                            })
                        }

                        CorrelType::ZScore => {
                            let (m, sd) = tmp.1.as_slice().stat(0);
                            tmp.1.iter_mut().for_each(|x| {
                                *x = (*x - m) / sd;
                            })
                        }
                        CorrelType::Rank => {} // do nothing
                        _ => panic!("unsupported method"),
                    }
                }
                // if weight == 0, ranked values turns to 1 automatically
                tmp.1.iter_mut().for_each(|x| {
                    *x = x.abs().powf(self.weight);
                });
                return tmp;
            })
            .collect();
        // save indices
        weighted_sorted_metric.iter().for_each(|(idx, m)| {
            self.indices.push(idx.to_owned());
        });
        let es = EnrichmentScore::new(genes, self.nperm, self.seed, true, false);
        // let end1 = Instant::now();
        for (&term, &gset) in gmt.iter() {
            let tag = es.gene.isin(gset);
            let hit = tag.iter().filter(|&x| x > &0.0).count();
            if hit > self.max_size || hit < self.min_size {
                continue;
            }
            let mut summ: Vec<GSEASummary> = weighted_sorted_metric
                .par_iter()
                .enumerate()
                .map(|(i, (indices, metric))| {
                    let tag_new: Vec<f64> = indices.iter().map(|&idx| tag[idx]).collect();
                    let gidx = es.hit_index(&tag_new);
                    // let run_es = es.running_enrichment_score(metric, &tag_new);
                    // let es1 = es.select_es(&run_es);
                    // faster version of ssGSEA
                    let es2 = es.fast_random_walk_ss(metric, &tag_new);
                    GSEASummary {
                        term: term.to_string(),
                        es: es2,
                        run_es: Vec::<f64>::new(), // run_es,
                        hits: gidx,                // gene hit idx of each sample after sorting
                        index: Some(i),
                        ..Default::default()
                    }
                })
                .collect();
            self.summaries.append(&mut summ);
        }
        // let end2 = Instant::now();
        // println!("Calculation time: {:.2?}", end2.duration_since(end1));
        // self.stat(); // NES
        let max = self
            .summaries
            .iter()
            .fold(std::f64::MIN, |a, b| a.max(b.es));
        let min = self
            .summaries
            .iter()
            .fold(std::f64::MAX, |a, b| a.min(b.es));
        let norm = max - min;
        self.summaries.iter_mut().for_each(|b| b.nes = b.es / norm);

        // let end3 = Instant::now();
        // println!("Statistical time: {:.2?}", end3.duration_since(end2));
    }

    /// single sample gsea
    pub fn ss_gsea_permuate(
        &mut self,
        genes: &[String],
        gene_exp: &[Vec<f64>], // 2d vector [m_gene, n_sample];
        gmt: &BTreeMap<&str, &[String]>,
        correl_type: CorrelType,
    ) {
        // transpose [m_gene, n_sample] --> [n_sample, m_gene]
        let mut gene_metric: Vec<Vec<f64>> = vec![vec![]; gene_exp[0].len()];
        gene_exp.iter().for_each(|row| {
            row.iter().enumerate().for_each(|(j, e)| {
                gene_metric[j].push(*e);
            });
        });

        // sort first and then set weight,
        let weighted_sorted_metric: Vec<(Vec<usize>, Vec<f64>)> = gene_metric
            .into_par_iter()
            .map(|rank| {
                let mut tmp = rank.as_slice().argsort(false);
                if self.weight > 0.0 {
                    // https://github.com/broadinstitute/ssGSEA2.0/blob/f682082f62ae34185421545f25041bae2c78c89b/src/ssGSEA2.0.R#L396
                    match correl_type {
                        CorrelType::SymRank => {
                            let idx = (tmp.1.len() + 2 - 1) / 2;
                            let mid = tmp.1.get(idx).unwrap().to_owned();
                            tmp.1.iter_mut().for_each(|x| {
                                if *x > mid {
                                    *x
                                } else {
                                    (*x) + (*x) - mid
                                };
                            })
                        }

                        CorrelType::ZScore => {
                            let (m, sd) = tmp.1.as_slice().stat(0);
                            tmp.1.iter_mut().for_each(|x| {
                                *x = (*x - m) / sd;
                            })
                        }
                        CorrelType::Rank => {}
                        _ => panic!("unsupported method"),
                    }
                }
                tmp.1.iter_mut().for_each(|x| {
                    *x = x.abs().powf(self.weight);
                });
                return tmp;
            })
            .collect();
        // save indices
        weighted_sorted_metric.iter().for_each(|(idx, m)| {
            self.indices.push(idx.to_owned());
        });
        // build genes permutations
        let mut es = EnrichmentScore::new(genes, self.nperm, self.seed, true, false);
        let gperm = es.gene_permutation(); // gene permutation
        let mut _all = Vec::<GSEASummary>::new();
        // let end1 = Instant::now();
        weighted_sorted_metric
            .into_iter()
            .enumerate()
            .for_each(|(i, (indices, metric))| {
                // update the order of genes
                let _genes: Vec<String> =
                    indices.into_iter().map(|j| genes[j].to_string()).collect();
                let od_genes = DynamicEnum::from(&_genes);
                // write summary
                let mut summ = Vec::<GSEASummary>::new();
                for (&term, &gset) in gmt.iter() {
                    // update tag indicator
                    let gtag = od_genes.isin(gset);
                    let gidx = es.hit_index(&gtag);
                    if gidx.len() > self.max_size || gidx.len() < self.min_size {
                        continue;
                    }
                    // note: update first element of gperm to get correct order of the gene ranking
                    let mut tag_indicators: Vec<Vec<f64>> =
                        gperm.par_iter().map(|de| de.isin(&gidx)).collect();
                    tag_indicators[0] = gtag; // update
                                              // get runing enrichment score
                    let run_es: Vec<f64> = es.running_enrichment_score(&metric, &tag_indicators[0]);
                    let ess: Vec<f64> = tag_indicators
                        .par_iter()
                        .map(|tag| {
                            // calculate ES
                            es.fast_random_walk_ss(&metric, tag)
                        })
                        .collect();

                    let esnull: Vec<f64> = if ess.len() > 1 {
                        ess[1..].to_vec()
                    } else {
                        Vec::new()
                    };
                    let gsu = GSEASummary {
                        term: term.to_string(),
                        es: ess[0],
                        run_es: run_es,
                        hits: gidx, // hit index of each sample after sorting
                        esnull: esnull,
                        index: Some(i),
                        ..Default::default()
                    };
                    summ.push(gsu);
                }
                // calculate nes, pval, fdr
                self.stat(&mut summ);
                _all.append(&mut summ);
            });

        self.summaries = _all;
    }
}
/// Extension of `GSEAResult` providing the fgsea multilevel p-value variant.
///
/// `prerank_multilevel()` ports `fgseaMultilevel` from the fgsea R/C++ package:
///
/// - **p-value**: adaptive multilevel splitting + MCMC algorithm (arbitrarily precise,
///   with quantified error bound `log2err`).
/// - **NES**: `ES / mean(null ES with same sign)`, derived from `n_perm_simple` simple
///   gene permutations — identical to the fgsea normalization.
/// - **FDR**: Benjamini-Hochberg correction on the multilevel p-values.
/// - **FWER**: not applicable (set to 1.0).
impl GSEAResult {
    /// Benjamini-Hochberg FDR correction.
    ///
    /// `pvals`: p-values in arbitrary order (each must be in `[0, 1]`).
    /// Returns a `Vec<f64>` of BH-adjusted q-values in the **same order** as the
    /// input, with each q-value guaranteed to be in `[0, 1]` and monotone with
    /// respect to the input p-value ordering (smaller p ⇒ smaller or equal q).
    fn bh_correction(pvals: &[f64]) -> Vec<f64> {
        let m = pvals.len();
        if m == 0 {
            return Vec::new();
        }
        // Build (pval, original_index) pairs sorted ascending by p-value
        let mut order: Vec<usize> = (0..m).collect();
        order.sort_unstable_by(|&a, &b| {
            pvals[a].partial_cmp(&pvals[b]).unwrap_or(std::cmp::Ordering::Equal)
        });

        // BH adjustment: q[rank] = p[rank] * m / rank
        // Enforce monotonicity by scanning from the largest rank downward and
        // keeping track of the running minimum.  All input p-values are in
        // [0, 1] so `p * m / rank` is also in [0, m] — the min() pass keeps
        // q within [0, 1].
        let mut adjusted = vec![1.0f64; m];
        let mut min_so_far = 1.0f64;
        for (rank_minus1, &orig_idx) in order.iter().enumerate().rev() {
            let rank = (rank_minus1 + 1) as f64;
            // p * m / rank is non-negative; cap at 1.0 to stay in [0,1]
            let q = (pvals[orig_idx] * m as f64 / rank).min(1.0);
            min_so_far = min_so_far.min(q);
            adjusted[orig_idx] = min_so_far;
        }
        adjusted
    }

    /// Preranked GSEA using the fgsea multilevel algorithm for p-value computation.
    ///
    /// This method ports `fgseaMultilevel` from the fgsea R/C++ package.
    ///
    /// **Two-phase computation (matching fgsea):**
    /// 1. *NES*: `n_perm_simple` simple gene permutations build a per-gene-set null ES
    ///    distribution. `NES = ES / mean(positive null ESs)` when ES ≥ 0, or
    ///    `ES / |mean(negative null ESs)|` when ES < 0.
    /// 2. *p-value*: adaptive multilevel splitting + MCMC yields arbitrarily precise
    ///    values with quantified error bound `log2err`.
    /// 3. *FDR*: Benjamini-Hochberg correction across all tested gene sets.
    ///
    /// Parameters
    /// ----------
    /// - `genes`: gene names in rank-descending order.
    /// - `metric`: gene-level ranking metric (descending, all values real).
    /// - `gmt`: gene set dictionary.
    /// - `sample_size`: MCMC sample size per level (fgsea default: 101).
    /// - `n_perm_simple`: number of simple gene permutations for NES normalization
    ///   (fgsea default: 1000; set to 0 to skip, in which case `nes = es`).
    /// - `eps`: convergence threshold for the multilevel algorithm (e.g. 1e-50).
    pub fn prerank_multilevel(
        &mut self,
        genes: &[String],
        metric: &[f64],
        gmt: &BTreeMap<&str, &[String]>,
        sample_size: usize,
        eps: f64,
    ) {
        // Build integer-scaled ranks for the multilevel algorithm.
        let weighted_metric: Vec<f64> =
            metric.iter().map(|x| x.abs().powf(self.weight)).collect();
        let int_ranks = scale_ranks(&weighted_metric);

        // Observed ES + running-ES curve per gene set. No gene-label permutation
        // null is built here: NES normalization uses fgsea's random-gene-set null
        // (calcGseaStatCumulativeBatch), computed in one batch after this loop.
        let es = EnrichmentScore::new(genes, 0, self.seed, false, false);

        // Collect per-gene-set results; BH correction requires all p-values.
        let mut summ = Vec::<GSEASummary>::new();
        // Parallel arrays (aligned 1:1 with `summ`) feeding the fgsea NES null.
        let mut pathway_sizes = Vec::<i32>::new();
        let mut pathway_scores = Vec::<f64>::new();

        for (&term, &gset) in gmt.iter() {
            let gtag = es.gene.isin(gset);
            let gidx = es.hit_index(&gtag);
            if gidx.len() > self.max_size || gidx.len() < self.min_size {
                continue;
            }

            // Observed enrichment score (signed) and running-ES curve.
            let observed_es = es.fast_random_walk(&weighted_metric, &gtag);
            let run_es = es.running_enrichment_score(&weighted_metric, &gtag);

            // p-value is deferred: fgsea needs the simple-permutation null (computed in
            // one batch below) both to normalise the multilevel estimate and to decide
            // which estimator to use for each gene set.
            pathway_sizes.push(gidx.len() as i32);
            pathway_scores.push(observed_es);
            summ.push(GSEASummary {
                term: term.to_string(),
                es: observed_es,
                pval: f64::NAN,
                log2err: f64::NAN,
                fwerp: 1.0,
                run_es,
                hits: gidx,
                ..Default::default()
            });
        }

        let mut n_unbalanced = 0usize;
        let mut n_overestimated = 0usize;

        // NES via fgsea's random-gene-set null (calcGseaStatCumulativeBatch):
        // for each of `nperm` random size-k samples, the per-size prefix ES is
        // accumulated, giving geZeroSum/geZero = mean positive null ES and
        // leZeroSum/leZero = mean negative null ES. NES = ES / mean(same-sign null),
        // exactly as fgsea's `fgseaMultilevel`. The RNG is boost-MT19937 (rand_mt),
        // so NES is bit-exact with fgsea for matching stats, sizes, nPermSimple
        // (= nperm), seed, weight (= gseaParam), and scoreType ("std").
        if self.nperm > 0 && !summ.is_empty() {
            let batch = calcGseaStatCumulativeBatch(
                metric, // signed gene-level stats, sorted descending
                self.weight, // gseaParam
                &pathway_scores,
                &pathway_sizes,
                self.nperm as i32,
                self.seed as i32,
                "std",
            );
            for (j, s) in summ.iter_mut().enumerate() {
                let ge_mean = if batch.geZero[j] > 0.0 {
                    batch.geZeroSum[j] / batch.geZero[j]
                } else {
                    0.0
                };
                let le_mean = if batch.leZero[j] > 0.0 {
                    batch.leZeroSum[j] / batch.leZero[j]
                } else {
                    0.0
                };
                // NES is normalised by the mean null ES of the SAME sign. `ES > 0`
                // (not `>= 0`) matches fgseaSimpleImpl, which sends a zero ES to the
                // negative branch.
                let nes_mean = if s.es > 0.0 { ge_mean } else { le_mean.abs() };
                s.nes = if nes_mean != 0.0 { s.es / nes_mean } else { f64::NAN };

                // ---- p-value, following fgseaMultilevel's two-estimator scheme ----
                //
                // `modeFraction` counts the permutations whose null ES has the same
                // sign as the observed ES; `nMoreExtreme` counts those at least as
                // extreme as it. (fgsea uses `ES >= 0` for the former and `ES > 0` for
                // the latter; both are reproduced verbatim.)
                let n_perm = self.nperm as f64;
                let mode_fraction = if s.es >= 0.0 {
                    batch.geZero[j]
                } else {
                    batch.leZero[j]
                };
                let n_more_extreme = if s.es > 0.0 {
                    batch.geEs[j]
                } else {
                    batch.leEs[j]
                };

                // No same-sign null mass to normalise against, or too little of it:
                // fgsea reports pval/NES/log2err as NA rather than a number it cannot
                // stand behind.
                if s.nes.is_nan() || mode_fraction < 10.0 {
                    n_unbalanced += 1;
                    s.pval = f64::NAN;
                    s.nes = f64::NAN;
                    s.log2err = f64::NAN;
                    continue;
                }

                // fgsea compares the two estimators' accuracy on the crude estimate
                // (nMoreExtreme + 1) / (nPerm + 1), not on the reported p-value.
                let crude = (n_more_extreme + 1.0) / (n_perm + 1.0);
                let simple_err = simple_error_log2(n_more_extreme, n_perm);
                let mult_err = multilevel_error_pval(crude, sample_size);

                if mult_err >= simple_err {
                    // The simple permutation estimate is at least as accurate here as
                    // the multilevel one would be, so fgsea keeps it. This is the
                    // common case: for all but the smallest p-values, `nperm` draws
                    // resolve the tail better than multilevel sampling does.
                    //
                    // Each side is normalised by its OWN sign's null count, and the
                    // smaller of the two is taken -- normalising by `n_perm` instead
                    // would bias every p-value low by the same-sign fraction (~2x).
                    s.pval = ((1.0 + batch.leEs[j]) / (1.0 + batch.leZero[j]))
                        .min((1.0 + batch.geEs[j]) / (1.0 + batch.geZero[j]));
                    s.log2err = ((n_more_extreme + 1.0).trigamma() - (n_perm + 1.0).trigamma())
                        .sqrt()
                        / std::f64::consts::LN_2;
                } else {
                    // Deep tail: fall back to adaptive multilevel splitting, which can
                    // resolve p-values far below 1/nperm.
                    let (cond_pval, is_cp_ge_half, _) = compute_pvalue_multilevel(
                        &int_ranks,
                        s.hits.len(),
                        s.es,
                        sample_size,
                        self.seed,
                        eps,
                    );
                    // The multilevel core returns a probability CONDITIONAL on the null
                    // ES sharing the observed sign. Dividing by that sign's frequency
                    // recovers the unconditional p-value; without this the result is
                    // biased low by ~1/denom_prob (roughly 2x for a balanced metric).
                    let denom_prob = (mode_fraction + 1.0) / (n_perm + 1.0);
                    s.pval = (cond_pval / denom_prob).min(1.0);
                    s.log2err = if is_cp_ge_half {
                        multilevel_error_pval(s.pval, sample_size)
                    } else {
                        // p-value likely overestimated; fgsea declines to report an error bound.
                        n_overestimated += 1;
                        f64::NAN
                    };
                }

                if s.pval < eps {
                    s.pval = eps;
                    s.log2err = f64::NAN;
                }
            }
        } else {
            // No simple-permutation null requested. Without it there is no
            // `denomProb` to normalise the multilevel estimate against, so fall back
            // to the raw conditional p-value and an unnormalised NES.
            for s in summ.iter_mut() {
                s.nes = s.es;
                let (cond_pval, _, ml_log2err) = compute_pvalue_multilevel(
                    &int_ranks,
                    s.hits.len(),
                    s.es,
                    sample_size,
                    self.seed,
                    eps,
                );
                s.pval = cond_pval;
                s.log2err = ml_log2err;
            }
        }

        if n_unbalanced > 0 {
            eprintln!(
                "Warning: {} gene sets had fewer than 10 same-sign null enrichment scores \
                 (unbalanced gene-level statistics); their pval, FDR, NES and log2err are NaN. \
                 Increasing the number of permutations may help.",
                n_unbalanced
            );
        }
        if n_overestimated > 0 {
            eprintln!(
                "Warning: for {} gene sets the p-value was likely overestimated; \
                 log2err is set to NaN.",
                n_overestimated
            );
        }

        // Benjamini-Hochberg correction for FDR. Gene sets with a NaN p-value are
        // excluded from the correction (and keep a NaN FDR), matching p.adjust's
        // handling of NAs -- leaving them in would corrupt the rank ordering.
        let scored: Vec<usize> = (0..summ.len()).filter(|&i| !summ[i].pval.is_nan()).collect();
        let pvals: Vec<f64> = scored.iter().map(|&i| summ[i].pval).collect();
        let fdrs = Self::bh_correction(&pvals);
        for s in summ.iter_mut() {
            s.fdr = f64::NAN;
        }
        for (&i, &fdr) in scored.iter().zip(fdrs.iter()) {
            summ[i].fdr = fdr;
        }

        self.summaries = summ;
        self.indices.push((0..genes.len()).collect_vec());
        self.rankings.push(metric.to_owned());
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;
    // use fastrand;
    use crate::stats::GSEAResult;
    use crate::utils::FileReader;
    #[test]
    fn test_prerank() {
        let cwd = std::env::current_dir().unwrap(); // prjoject root, directory to Cargo.toml
        let rnk_path = cwd.join("tests/data/mds.2k.rnk");
        let gmt_path = cwd.join("tests/data/hallmark.gmt");

        let start = Instant::now();
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()
            .unwrap();
        let mut rnk = FileReader::new();
        let _ = rnk.read_csv(rnk_path.to_str().unwrap(), b'\t', false, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table(gmt_path.to_str().unwrap(), '\t', false);

        // let gene: Vec<String> = vec!["A","B","C","D","E","F","G","H","J","K"].into_iter().map(|s| s.to_string()).collect();
        // let gene_set: Vec<String> = vec!["B","A","D","G"].into_iter().map(|s| s.to_string()).collect();
        // let gene_metric = vec![9.0,4.0,3.0,2.0,1.0,0.5,0.1,-0.1,-0.2,-0.5];
        let weight = 1.0;
        let mut gene: Vec<String> = Vec::new();
        // let mut gene_set: Vec<String> = Vec::new();
        let mut gene_metric: Vec<f64> = Vec::new();
        for r in rnk.record.iter() {
            gene.push(r[0].clone());
            gene_metric.push(r[1].parse::<f64>().unwrap());
        }

        // hashmap
        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {
            gmt2.insert(r[0].as_str(), &r[2..]);
        });

        // weighted then sort
        gene_metric
            .iter_mut()
            .for_each(|x| *x = x.abs().powf(weight));
        let (gidx, metric) = gene_metric.as_slice().argsort(false);
        gene = gidx.iter().map(|&i| gene[i].clone()).collect();
        // start to calculate
        let mut gsea = GSEAResult::new(weight, 500, 5, 10, 123);
        gsea.prerank(&gene, &metric, &gmt2);
        let end = Instant::now();
        println!("Overall run time: {:.2?}", end.duration_since(start));
        println!("This version 1");
        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });
    }
    #[test]
    fn test_gsea() {
        let cwd = std::env::current_dir().unwrap(); // prjoject root, directory to Cargo.toml
        let gct_path = cwd.join("tests/extdata/Leukemia_hgu95av2.trim2.txt");
        let gmt_path = cwd.join("tests/extdata/h.all.v7.0.symbols.gmt");
        let cls_path = cwd.join("tests/extdata/Leukemia.cls");
        let start = Instant::now();
        // set number of threads of rayon at the main()
        // rayon::ThreadPoolBuilder::new()
        //     .num_threads(1)
        //     .build_global()ƒ
        //     .unwrap();

        let mut gct = FileReader::new();
        let _ = gct.read_csv(gct_path.to_str().unwrap(), b'\t', true, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table(gmt_path.to_str().unwrap(), '\t', false);
        let mut cls = FileReader::new();
        let _ = cls.read_table(cls_path.to_str().unwrap(), ' ', false);
        println!("{:?}", &cls.record[2]);
        let gboo: Vec<bool> = cls.record[2].iter().map(|x| x != "AML").collect();
        println!("{:?}", &gboo);
        let weight = 1.0;
        let mut gene: Vec<String> = Vec::new();
        // let mut gene_set: Vec<String> = Vec::new();
        let mut gene_exp: Vec<Vec<f64>> = Vec::new();
        for r in gct.record.iter() {
            gene.push(r[0].to_string());
            let mut vv: Vec<f64> = Vec::new();
            for v in &r[2..] {
                vv.push(v.parse::<f64>().unwrap());
            }
            gene_exp.push(vv);
        }

        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {
            gmt2.insert(r[0].as_str(), &r[2..]);
        });

        let mut gsea = GSEAResult::new(weight, 1000, 3, 10, 123);
        gsea.gsea(&gene, &gboo, &gene_exp, &gmt2, Metric::Signal2Noise);

        let end = Instant::now();
        println!("Overall run time: {:.2?}", end.duration_since(start));
        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.5?}, fdr: {:.5?}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });
    }

    #[test]
    fn test_ssgsea() {
        let cwd = std::env::current_dir().unwrap(); // prjoject root, directory to Cargo.toml
        let gct_path = cwd.join("tests/extdata/Leukemia_hgu95av2.trim2.txt");
        let gmt_path = cwd.join("tests/extdata/h.all.v7.0.symbols.gmt");
        let cls_path = cwd.join("tests/extdata/Leukemia.cls");

        let mut gct = FileReader::new();
        let _ = gct.read_csv(gct_path.to_str().unwrap(), b'\t', true, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table(gmt_path.to_str().unwrap(), '\t', false);
        let mut cls = FileReader::new();
        let _ = cls.read_table(cls_path.to_str().unwrap(), ' ', false);
        println!("{:?}", &cls.record[2]);
        let gboo: Vec<bool> = cls.record[2].iter().map(|x| x != "AML").collect();
        println!("{:?}", &gboo);
        let weight = 1.0;
        let mut gene: Vec<String> = Vec::new();
        // let mut gene_set: Vec<String> = Vec::new();
        let mut gene_exp: Vec<Vec<f64>> = Vec::new();
        for r in gct.record.iter() {
            gene.push(r[0].to_string());
            let mut vv: Vec<f64> = Vec::new();
            for v in &r[2..] {
                vv.push(v.parse::<f64>().unwrap());
            }
            gene_exp.push(vv);
        }

        let sample_names = &gct.header.get_vec()[2..];

        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {
            gmt2.insert(r[0].as_str(), &r[2..]);
        });

        let nperm = 10;
        let mut gsea = GSEAResult::new(weight, 500, 3, nperm, 123);
        if nperm > 0 {
            gsea.ss_gsea_permuate(&gene, &gene_exp, &gmt2, CorrelType::Rank);
        } else {
            gsea.ss_gsea(&gene, &gene_exp, &gmt2, CorrelType::Rank);
        }

        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });

        let _gx: Vec<f64> = gene_exp.iter().map(|g| g[0]).collect();
        let (_gx_idx, _gx2) = _gx.as_slice().argsort(false);
        let _ge: Vec<String> = _gx_idx.into_iter().map(|i| gene[i].to_owned()).collect();
        gsea.prerank(&_ge, &_gx2, &gmt2);
        println!("\n\n\nThis version prerank version 1");
        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });
        println!("\n\n\nThis version prerank version 2");
        gsea.prerank2(&gene, &gene_exp, &gmt2);
        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });
    }

    #[test]
    fn test_bh_correction() {
        // Floating-point tolerance used across the BH tests
        const EPSILON: f64 = 1e-12;
        // Known BH correction example
        // p-values from Benjamini & Hochberg 1995 example
        let pvals = vec![0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216];
        let adj = GSEAResult::bh_correction(&pvals);
        // All adjusted values should be ≤ 1.0 and ≥ original p-values
        for (p, q) in pvals.iter().zip(adj.iter()) {
            assert!(*q >= *p - EPSILON, "BH q {q} < p {p}");
            assert!(*q <= 1.0 + EPSILON, "BH q {q} > 1");
        }
        // Smallest p-value should produce the smallest q-value
        let min_q = adj.iter().cloned().fold(f64::MAX, f64::min);
        assert!(adj[0] == min_q, "smallest p should have smallest q");
        // Monotonicity: q-values sorted in the same rank order as p-values (ascending p → ascending q)
        let mut order: Vec<usize> = (0..pvals.len()).collect();
        order.sort_unstable_by(|&a, &b| pvals[a].partial_cmp(&pvals[b]).unwrap());
        for w in order.windows(2) {
            assert!(adj[w[0]] <= adj[w[1]] + EPSILON,
                "BH monotonicity violated: adj[{}]={} > adj[{}]={}", w[0], adj[w[0]], w[1], adj[w[1]]);
        }
    }

    #[test]
    fn test_prerank_multilevel() {
        const EPSILON: f64 = 1e-12;
        let cwd = std::env::current_dir().unwrap();
        let rnk_path = cwd.join("tests/data/mds.2k.rnk");
        let gmt_path = cwd.join("tests/data/c2.cp.kegg.v7.5.1.symbols.gmt");
        let mut rnk = FileReader::new();
        let _ = rnk.read_csv(rnk_path.to_str().unwrap(), b'\t', false, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table(gmt_path.to_str().unwrap(), '\t', false);
        let weight = 1.0;
        let mut gene: Vec<String> = Vec::new();
        let mut gene_metric: Vec<f64> = Vec::new();
        for r in rnk.record.iter() {
            gene.push(r[0].clone());
            gene_metric.push(r[1].parse::<f64>().unwrap());
        }
        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {
            gmt2.insert(r[0].as_str(), &r[2..]);
        });
        gene_metric.iter_mut().for_each(|x| *x = x.abs().powf(weight));
        let (gidx, metric) = gene_metric.as_slice().argsort(false);
        gene = gidx.iter().map(|&i| gene[i].clone()).collect();

        // n_perm_simple=100 for fast testing; use 1000 in practice
        let n_perm_simple = 100usize;
        let mut gsea = GSEAResult::new(weight, 500, 5, n_perm_simple, 123);
        gsea.prerank_multilevel(&gene, &metric, &gmt2, 51, 1e-10);

        // Basic sanity checks
        assert!(!gsea.summaries.is_empty(), "expected some gene sets");
        for s in &gsea.summaries {
            assert!(s.pval >= 0.0 && s.pval <= 1.0, "pval out of [0,1]: {}", s.pval);
            assert!(s.fdr >= 0.0 && s.fdr <= 1.0, "fdr out of [0,1]: {}", s.fdr);
            // esnull should be cleared from output after NES is computed
            assert!(s.esnull.is_empty(), "esnull should be empty in output");
            // fwerp should be 1.0 (not applicable without permutation null)
            assert!((s.fwerp - 1.0).abs() < EPSILON, "fwerp should be 1.0");
        }
        // BH property: FDR >= p-value for each gene set
        for s in &gsea.summaries {
            assert!(s.fdr >= s.pval - EPSILON, "fdr {:.2e} < pval {:.2e}", s.fdr, s.pval);
        }
        // NES should differ from raw ES when n_perm_simple > 0 (permutation normalisation)
        // (at least one gene set should have NES != ES due to normalisation)
        let any_nes_differs = gsea.summaries.iter().any(|s| (s.nes - s.es).abs() > 1e-6);
        assert!(any_nes_differs, "NES should differ from ES when n_perm_simple > 0");

        // NES sign should match ES sign for all gene sets
        for s in &gsea.summaries {
            assert!(s.nes * s.es >= 0.0, "NES and ES should have the same sign: nes={}, es={}", s.nes, s.es);
        }

        println!("prerank_multilevel: {} gene sets", gsea.summaries.len());
        gsea.summaries.iter().take(5).for_each(|g| {
            println!("  term={}, es={:.4}, nes={:.4}, pval={:.2e}, fdr={:.2e}, log2err={:.2}",
                g.term, g.es, g.nes, g.pval, g.fdr, g.log2err);
        });

        // Also test n_perm_simple=0 fallback: NES should equal ES
        let mut gsea0 = GSEAResult::new(weight, 500, 5, 0, 123);
        gsea0.prerank_multilevel(&gene, &metric, &gmt2, 51, 1e-10);
        for s in &gsea0.summaries {
            assert!((s.nes - s.es).abs() < EPSILON, "with n_perm_simple=0, nes should equal es");
        }
    }

    /// Load the bundled prerank fixtures, ready for `prerank_multilevel`.
    ///
    /// The metric keeps its SIGN and is sorted descending, which is what
    /// `prerank_multilevel` expects -- it applies `|x|^weight` internally for the
    /// enrichment score and needs the signed values for the permutation null. Passing
    /// a pre-`abs()`ed metric makes every null enrichment score one-signed, which
    /// drives `denomProb` to 1 and hides any error in the sign normalisation.
    fn load_prerank_fixture() -> (Vec<String>, Vec<f64>, FileReader, FileReader) {
        let cwd = std::env::current_dir().unwrap();
        let mut rnk = FileReader::new();
        let _ = rnk.read_csv(
            cwd.join("tests/data/mds.2k.rnk").to_str().unwrap(), b'\t', false, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table(
            cwd.join("tests/data/c2.cp.kegg.v7.5.1.symbols.gmt").to_str().unwrap(), '\t', false);
        let mut gene: Vec<String> = Vec::new();
        let mut gene_metric: Vec<f64> = Vec::new();
        for r in rnk.record.iter() {
            gene.push(r[0].clone());
            gene_metric.push(r[1].parse::<f64>().unwrap());
        }
        let (gidx, metric) = gene_metric.as_slice().argsort(false);
        let gene = gidx.iter().map(|&i| gene[i].clone()).collect();
        (gene, metric, rnk, gmt)
    }

    /// Under the null, p-values must be approximately Uniform(0, 1).
    ///
    /// This is the property that a missing `denomProb` normalisation destroys. The
    /// multilevel core returns a probability CONDITIONAL on the null ES sharing the
    /// observed ES's sign; reporting it directly makes every p-value roughly
    /// `1 / denomProb` (~2x) too small. Shuffling the gene labels removes any real
    /// enrichment, so a correctly normalised implementation gives mean(p) ~ 0.5 and
    /// P(p < 0.05) ~ 0.05, while the unnormalised one gives ~0.29 and ~0.12.
    #[test]
    fn test_prerank_multilevel_null_calibration() {
        let weight = 1.0;
        let (gene, metric, _rnk, gmt) = load_prerank_fixture();
        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| { gmt2.insert(r[0].as_str(), &r[2..]); });

        // Deterministic label shuffle (LCG): keeps the metric distribution intact but
        // destroys any association between gene sets and the ranking.
        let mut shuffled = gene.clone();
        let mut state: u64 = 0x2545F491_4F6CDD1D;
        for i in (1..shuffled.len()).rev() {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            shuffled.swap(i, (state >> 33) as usize % (i + 1));
        }

        let mut gsea = GSEAResult::new(weight, 500, 15, 1000, 42);
        gsea.prerank_multilevel(&shuffled, &metric, &gmt2, 101, 1e-50);
        let p: Vec<f64> = gsea.summaries.iter().map(|s| s.pval).filter(|v| !v.is_nan()).collect();
        assert!(p.len() > 50, "expected a decent number of gene sets, got {}", p.len());

        let mean = p.iter().sum::<f64>() / p.len() as f64;
        let frac_sig = p.iter().filter(|&&v| v < 0.05).count() as f64 / p.len() as f64;
        println!("null calibration: n={} mean(p)={:.4} P(p<0.05)={:.4}", p.len(), mean, frac_sig);

        // Generous bounds: these must catch a ~2x bias, not police Monte-Carlo noise.
        // The unnormalised implementation lands at mean ~0.29 / P ~0.12.
        assert!((0.42..=0.58).contains(&mean),
                "null p-values are not uniform: mean(p)={:.4}, expected ~0.5. \
                 A low mean indicates the conditional multilevel p-value is not being \
                 divided by denomProb.", mean);
        assert!(frac_sig < 0.10,
                "type-I error inflated: P(p<0.05)={:.4}, expected ~0.05", frac_sig);
    }

    /// `eps` floors the reported p-value, and fgsea reports no error bound once the
    /// floor is hit (`result[pval < eps, c("pval","log2err") := list(eps, NA)]`).
    #[test]
    fn test_prerank_multilevel_eps_floor() {
        let weight = 1.0;
        let (gene, metric, _rnk, gmt) = load_prerank_fixture();
        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| { gmt2.insert(r[0].as_str(), &r[2..]); });

        let eps = 1e-3;
        let mut gsea = GSEAResult::new(weight, 500, 15, 1000, 42);
        gsea.prerank_multilevel(&gene, &metric, &gmt2, 101, eps);

        let mut clamped = 0;
        for s in &gsea.summaries {
            if s.pval.is_nan() { continue; }
            assert!(s.pval >= eps, "p-value {:.3e} fell below eps={:.0e}", s.pval, eps);
            if s.pval == eps {
                clamped += 1;
                assert!(s.log2err.is_nan(),
                        "log2err must be NaN once the eps floor is hit, got {}", s.log2err);
            }
        }
        assert!(clamped > 0, "expected at least one gene set to hit the eps floor");
    }

    /// Classical `prerank()` and `prerank_multilevel()` must agree exactly on ES, and
    /// only approximately on NES.
    ///
    /// ES is a deterministic function of the sorted metric and gene-set membership, so
    /// it is bit-for-bit equal between the two paths.
    ///
    /// NES is **not**, and must not be asserted equal: the two paths normalise against
    /// deliberately different null distributions. `prerank` divides by the mean
    /// same-sign ES of a *gene-label permutation* null, while `prerank_multilevel`
    /// divides by the mean same-sign ES of fgsea's *random-gene-set* null
    /// (`calcGseaStatCumulativeBatch`), which is what makes it fgsea-faithful. They
    /// estimate the same quantity from different samples, so they agree in sign and
    /// track each other closely, but never to the bit — measured over these 177 gene
    /// sets, zero have equal NES while the correlation is 0.997.
    ///
    /// An earlier version of this test asserted exact NES equality. That was valid when
    /// `prerank_multilevel` reused the classical `gperm` matrix, but stopped being true
    /// once the fgsea random-gene-set null landed; the assertion was left behind and had
    /// been failing ever since, reporting an arbitrary gene set each run because the
    /// comparison iterated a `HashMap`.
    #[test]
    fn test_prerank_vs_fgsea() {
        let weight = 1.0f64;
        // The shared fixture keeps the metric SIGNED. An earlier version of this test
        // sorted by |metric| instead, which made every null enrichment score one-signed
        // and hid exactly the sign-handling this comparison is meant to exercise.
        let (gene, metric, _rnk, gmt_file) = load_prerank_fixture();
        let mut gmt2 = BTreeMap::<&str, &[String]>::new();
        gmt_file.record.iter().for_each(|r| {
            gmt2.insert(r[0].as_str(), &r[2..]);
        });

        // Use a small permutation count so the test runs quickly.
        let n_perm = 50usize;
        let seed = 42u64;

        // --- classical prerank -----------------------------------------------
        let mut gsea_classical = GSEAResult::new(weight, 500, 5, n_perm, seed);
        gsea_classical.prerank(&gene, &metric, &gmt2);

        // --- fgsea multilevel (same n_perm, same seed) -----------------------
        let mut gsea_fgsea = GSEAResult::new(weight, 500, 5, n_perm, seed);
        gsea_fgsea.prerank_multilevel(&gene, &metric, &gmt2, 51, 1e-10);

        // Build maps: term → summary
        let classical_map: HashMap<&str, &GSEASummary> = gsea_classical
            .summaries.iter().map(|s| (s.term.as_str(), s)).collect();
        let fgsea_map: HashMap<&str, &GSEASummary> = gsea_fgsea
            .summaries.iter().map(|s| (s.term.as_str(), s)).collect();

        // Only compare terms present in both results (size filters are identical).
        // Sorted so that any failure names the same gene set on every run.
        let mut common_terms: Vec<&str> = classical_map.keys()
            .filter(|&&t| fgsea_map.contains_key(t)).copied().collect();
        common_terms.sort_unstable();

        assert!(!common_terms.is_empty(), "no gene sets in common between prerank and fgsea");

        // NES pairs, collected for the distribution-level checks below. Gene sets whose
        // NES is NaN are skipped: prerank_multilevel reports NaN when fewer than 10
        // permutations share the observed ES's sign, which fgsea also treats as
        // uninformative rather than as a value to compare.
        let mut pairs: Vec<(f64, f64)> = Vec::new();

        for term in &common_terms {
            let c = classical_map[term];
            let f = fgsea_map[term];

            // ES is deterministic — must be bit-for-bit equal.
            assert_eq!(c.es, f.es,
                "ES mismatch for '{}': classical={:.8}, fgsea={:.8}", term, c.es, f.es);

            if c.nes.is_nan() || f.nes.is_nan() {
                continue;
            }

            // Both nulls are same-sign normalised, so NES cannot flip sign relative to
            // ES. A disagreement here means one path picked the wrong side of the null.
            assert_eq!(c.nes.is_sign_positive(), f.nes.is_sign_positive(),
                "NES sign mismatch for '{}': classical={:.8}, fgsea={:.8}", term, c.nes, f.nes);

            pairs.push((c.nes, f.nes));
        }

        assert!(pairs.len() > 100, "expected >100 comparable gene sets, got {}", pairs.len());

        // The two NES estimates must track each other tightly. Bounds are deliberately
        // loose -- they exist to catch one path regressing to a different statistic, not
        // to police Monte-Carlo noise at only 50 permutations. Measured here: median
        // relative difference 0.033, correlation 0.997.
        let mut rel: Vec<f64> = pairs.iter()
            .map(|(c, f)| (f - c).abs() / c.abs())
            .collect();
        rel.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median_rel = rel[rel.len() / 2];

        let n = pairs.len() as f64;
        let mean_c = pairs.iter().map(|p| p.0).sum::<f64>() / n;
        let mean_f = pairs.iter().map(|p| p.1).sum::<f64>() / n;
        let cov: f64 = pairs.iter().map(|(c, f)| (c - mean_c) * (f - mean_f)).sum::<f64>();
        let var_c: f64 = pairs.iter().map(|(c, _)| (c - mean_c).powi(2)).sum::<f64>();
        let var_f: f64 = pairs.iter().map(|(_, f)| (f - mean_f).powi(2)).sum::<f64>();
        let corr = cov / (var_c.sqrt() * var_f.sqrt());

        println!("test_prerank_vs_fgsea: {} terms compared, {} with comparable NES; \
                  median rel diff {:.4}, correlation {:.6}",
                 common_terms.len(), pairs.len(), median_rel, corr);

        assert!(median_rel < 0.10,
            "classical and fgsea NES have diverged: median relative difference {:.4}, \
             expected < 0.10", median_rel);
        assert!(corr > 0.99,
            "classical and fgsea NES are no longer tracking: correlation {:.6}, \
             expected > 0.99", corr);

        // Print a few for visual inspection
        for term in common_terms.iter().take(5) {
            let c = classical_map[term];
            let f = fgsea_map[term];
            println!("  {} | classical es={:.6} nes={:.6} | fgsea es={:.6} nes={:.6}",
                term, c.es, c.nes, f.es, f.nes);
        }
    }
}
