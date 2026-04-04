#![allow(dead_code, unused)]

use crate::utils::DynamicEnum;
use crate::utils::{Metric, ScoreType, Statistic};
use pyo3::prelude::*;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;

/// Result of `calc_gsea_stat()`, mirroring fgsea's `calcGseaStat()` return value.
///
/// `tops` and `bottoms` are populated only when `return_all_extremes = true`.
/// `leading_edge` is populated only when `return_leading_edge = true`.
#[pyclass]
#[derive(Debug, Clone)]
pub struct GseaStatResult {
    /// Enrichment score (scalar summary statistic).
    #[pyo3(get)]
    pub es: f64,
    /// Per-hit running ES *after* each hit step (length k). Used for enrichment plots.
    #[pyo3(get)]
    pub tops: Vec<f64>,
    /// Per-hit running ES *before* each hit step, i.e. after the preceding miss run
    /// (length k). Used for enrichment plots.
    #[pyo3(get)]
    pub bottoms: Vec<f64>,
    /// 0-based indices (into the original `stats` slice) of leading-edge genes.
    #[pyo3(get)]
    pub leading_edge: Vec<usize>,
}

/// O(k) enrichment score computation — a direct port of fgsea's `calcGseaStat()`.
///
/// # Arguments
/// * `stats`               — gene-level ranking metric, **sorted descending**, length N.
///                           Values are raised to the power `gsea_param` internally.
/// * `selected`            — **0-based** indices of gene-set members in `stats`,
///                           **sorted ascending**, length k.
/// * `gsea_param`          — weight exponent p (1.0 for standard weighted ES).
/// * `score_type`          — which extreme to return:
///                           `Std` = sign-aware max |ES|, `Pos` = max, `Neg` = min.
/// * `return_all_extremes` — populate `tops` / `bottoms` in the result (for plots).
/// * `return_leading_edge` — populate `leading_edge` in the result.
///
/// # Panics
/// Panics if `selected` is empty or contains an index ≥ `stats.len()`.
pub fn calc_gsea_stat(
    stats: &[f64],
    selected: &[usize],
    gsea_param: f64,
    score_type: ScoreType,
    return_all_extremes: bool,
    return_leading_edge: bool,
) -> GseaStatResult {
    let n = stats.len();
    let m = selected.len();
    assert!(m > 0, "calc_gsea_stat: gene set must not be empty");
    assert!(m < n, "calc_gsea_stat: gene set must be smaller than stats");

    // Adjusted weights at each hit position: |stats[S[i]]|^p
    let r_adj: Vec<f64> = selected
        .iter()
        .map(|&i| stats[i].abs().powf(gsea_param))
        .collect();

    let nr: f64 = r_adj.iter().sum();

    // Cumulative sum of adjusted weights, normalised by NR.
    // Degenerate case (all weights zero) → uniform weights (equivalent to r_adj = 1/m each).
    let (r_cumsum, r_adj_norm): (Vec<f64>, Vec<f64>) = if nr == 0.0 {
        let uniform = 1.0 / m as f64;
        (
            (1..=m).map(|i| i as f64 * uniform).collect(),
            vec![uniform; m],
        )
    } else {
        let mut cum = 0.0;
        let cumsum: Vec<f64> = r_adj.iter().map(|&v| { cum += v / nr; cum }).collect();
        let norm: Vec<f64> = r_adj.iter().map(|&v| v / nr).collect();
        (cumsum, norm)
    };

    // For each hit i at position selected[i]:
    //   number of misses before hit i  = selected[i] - i   (0-based arithmetic)
    //   tops[i]    = r_cumsum[i] - misses_before_i / (N - m)
    //   bottoms[i] = tops[i] - r_adj_norm[i]
    let inv_miss = 1.0 / (n - m) as f64;
    let tops: Vec<f64> = r_cumsum
        .iter()
        .enumerate()
        .map(|(i, &c)| c - (selected[i] as f64 - i as f64) * inv_miss)
        .collect();
    let bottoms: Vec<f64> = tops
        .iter()
        .zip(r_adj_norm.iter())
        .map(|(&t, &r)| t - r)
        .collect();

    let max_p = tops.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_p = bottoms.iter().cloned().fold(f64::INFINITY, f64::min);

    let es = match score_type {
        ScoreType::Pos => max_p,
        ScoreType::Neg => min_p,
        ScoreType::Std => {
            if (max_p + min_p).abs() < 1e-12 {
                0.0
            } else if max_p > -min_p {
                max_p
            } else {
                min_p
            }
        }
    };

    // Leading edge: genes from the original stats array that drive the enrichment.
    let leading_edge = if return_leading_edge {
        match score_type {
            ScoreType::Pos => {
                let peak = tops
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                selected[..=peak].to_vec()
            }
            ScoreType::Neg => {
                let peak = bottoms
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                let mut le: Vec<usize> = selected[peak..].to_vec();
                le.reverse();
                le
            }
            ScoreType::Std => {
                if max_p > -min_p + 1e-12 {
                    let peak = tops
                        .iter()
                        .enumerate()
                        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                        .map(|(i, _)| i)
                        .unwrap_or(0);
                    selected[..=peak].to_vec()
                } else if -min_p > max_p + 1e-12 {
                    let peak = bottoms
                        .iter()
                        .enumerate()
                        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                        .map(|(i, _)| i)
                        .unwrap_or(0);
                    let mut le: Vec<usize> = selected[peak..].to_vec();
                    le.reverse();
                    le
                } else {
                    vec![]
                }
            }
        }
    } else {
        vec![]
    };

    GseaStatResult {
        es,
        tops: if return_all_extremes { tops } else { vec![] },
        bottoms: if return_all_extremes { bottoms } else { vec![] },
        leading_edge,
    }
}


/// Ensures std is at least 0.2 * abs(mean); if still 0, use 0.2.
/// Matches the original R/GSEA implementation.
#[inline]
fn sigma_correction(std: f64, mean: f64) -> f64 {
    let s = std.max(0.2 * mean.abs());
    if s == 0.0 { 0.2 } else { s }
}

pub trait EnrichmentScoreTrait {
    /// get run es only
    fn running_enrichment_score(&self, metric: &[f64], tag_indicator: &[f64]) -> Vec<f64>;
    /// fast GSEA only ES value return
    fn fast_random_walk(&self, metric: &[f64], tag_indicator: &[f64]) -> f64;
    /// fast ssGSEA. only ES value return
    fn fast_random_walk_ss(&self, metric: &[f64], tag_indicator: &[f64]) -> f64;
    /// calucalte metric, not sorted
    fn calculate_metric(&self, data: &[Vec<f64>], group: &[bool], method: Metric) -> Vec<f64>;
}

#[allow(dead_code)]
#[derive(Debug)]
pub struct EnrichmentScore {
    // ranking metric
    // metric: Vec<f64>,
    //metric2d: Vec<Vec<f64>>,
    pub gene: DynamicEnum<String>, //  Vec<String>, // gene names
    nperm: usize,                  // number of permutations
    single: bool,                  // single sample GSEA
    scale: bool,                   // whether to scale ES value
    rng: SmallRng,
}
// trait methods for EnrichmentScore
impl EnrichmentScoreTrait for EnrichmentScore {
    fn running_enrichment_score(&self, metric: &[f64], tag_indicator: &[f64]) -> Vec<f64> {
        let n: f64 = tag_indicator.len() as f64;
        let n_hint: f64 = tag_indicator.iter().sum();
        let n_miss: f64 = n - n_hint;
        let norm_notag: f64 = 1.0 / n_miss;
        let no_tag_indicator: Vec<f64> = tag_indicator.iter().map(|&b| 1.0 - b).collect();
        let sum_correl_tag: Vec<f64> = tag_indicator
            .iter()
            .zip(metric.iter())
            .map(|(&b, &v)| b * v)
            .collect();
        let norm_tag: f64 = 1.0 / sum_correl_tag.iter().sum::<f64>();
        // cumsum()
        let run_es: Vec<f64> = sum_correl_tag
            .iter()
            .zip(no_tag_indicator.iter())
            .map(|(&b, &v)| b * norm_tag - v * norm_notag)
            .scan(0.0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect();
        return run_es;
    }

    /// see here: https://github.com/ctlab/fgsea/blob/master/src/esCalculation.cpp
    fn fast_random_walk(&self, metric: &[f64], tag_indicator: &[f64]) -> f64 {
        // tag_indicator and metric must be sorted
        let ns: f64 = tag_indicator
            .iter()
            .zip(metric.iter())
            .map(|(&b, &v)| b * v)
            .sum::<f64>();
        let n: f64 = metric.len() as f64;
        let k: f64 = tag_indicator.iter().sum::<f64>() as f64;
        let mut res: f64 = 0.0; // running_es
        let mut cur: f64 = 0.0;
        let q1: f64 = 1.0 / (n - k);
        let q2: f64 = 1.0 / ns;
        let mut last: f64 = -1.0;
        let p: Vec<f64> = tag_indicator
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| if t > 0.0 { Some(i as f64) } else { None })
            .collect();

        for pos in p {
            cur -= q1 * (pos - last - 1.0);
            if cur.abs() > res.abs() {
                res = cur;
            }
            cur += q2 * metric.get(pos as usize).unwrap();
            if cur.abs() > res.abs() {
                res = cur;
            }
            last = pos;
        }

        // for pos in p {
        //     cur += q2 * metric.get(pos as usize).unwrap() - q1 * (pos - last - 1);
        //     res = max(res, cur);
        //     last = pos;
        // }
        return res;
    }

    fn fast_random_walk_ss(&self, metric: &[f64], tag_indicator: &[f64]) -> f64 {
        let n: f64 = metric.len() as f64;
        let k: f64 = tag_indicator.iter().sum();
        // idxs = np.flatnonzero(np.in1d(gene_ranking_index, gset)) # indices low to high
        let idxs: Vec<usize> = tag_indicator
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| if t > 0.0 { Some(i) } else { None })
            .collect();
        let _rnk: Vec<f64> = idxs
            .iter()
            .map(|&i| metric.get(i).unwrap().to_owned())
            .collect();
        let n_idxs: Vec<f64> = idxs.iter().map(|&i| n - (i as f64)).collect();
        // np.sum(ranking[idxs] * (n - idxs)) / np.sum(ranking[idxs])
        let step_cdf_in: f64 = _rnk
            .iter()
            .zip(n_idxs.iter())
            .map(|(r, i)| r * i)
            .sum::<f64>()
            / _rnk.iter().sum::<f64>();
        //  (n * (n + 1) /2 - np.sum(n-idxs)) / (n -k)
        let step_cdf_out: f64 = (n * (n + 1.0) / 2.0 - n_idxs.iter().sum::<f64>()) / (n - k);
        step_cdf_in - step_cdf_out
    }

    fn calculate_metric(&self, data: &[Vec<f64>], group: &[bool], method: Metric) -> Vec<f64> {
        data.iter()
            .map(|vec| {
                //let end1 = Instant::now();
                let mut pos: Vec<f64> = Vec::new();
                let mut neg: Vec<f64> = Vec::new();
                vec.iter().zip(group.iter()).for_each(|(&x, &b)| {
                    if b {
                        pos.push(x);
                    } else {
                        neg.push(x);
                    }
                });
                //let end2 = Instant::now();
                let pos_len = pos.len() as f64;
                let neg_len = neg.len() as f64;
                let (pos_mean, pos_std) = pos.as_slice().stat(1);
                let (neg_mean, neg_std) = neg.as_slice().stat(1);
                // let end3 = Instant::now();
                match method {
                    Metric::Signal2Noise | Metric::AbsSignal2Noise => {
                        // Apply GeneCluster sigma correction (same as original R/GSEA implementation):
                        // std = max(std, 0.2 * abs(mean)); if std is still 0, set to 0.2
                        let pos_std_c = sigma_correction(pos_std, pos_mean);
                        let neg_std_c = sigma_correction(neg_std, neg_mean);
                        let s2n = (pos_mean - neg_mean) / (pos_std_c + neg_std_c);
                        if method == Metric::AbsSignal2Noise { s2n.abs() } else { s2n }
                    }
                    Metric::Ttest => {
                        (pos_mean - neg_mean)
                            / (pos_std * pos_std / pos_len + neg_std * neg_std / neg_len).sqrt()
                    }
                    Metric::RatioOfClasses => pos_mean / neg_mean,
                    Metric::Log2RatioOfClasses => (pos_mean / neg_mean).log2(),
                    Metric::DiffOfClasses => pos_mean - neg_mean,
                }
            })
            .collect()
    }
}

/// associate and memeber function
impl EnrichmentScore {
    pub fn new(gene: &[String], nperm: usize, seed: u64, single: bool, scale: bool) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        EnrichmentScore {
            // metric: gene_metric,
            gene: DynamicEnum::from(gene),
            nperm: nperm + 1, // add 1 to kept track of original record
            single: single,
            scale: scale,
            rng: rng,
        }
    }
    pub fn hit_index(&self, tag_indicator: &[f64]) -> Vec<usize> {
        tag_indicator
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| if t > 0.0 { Some(i) } else { None })
            .collect()
        // let mut hit_ind: Vec<usize> = Vec::new();
        // for (i, b) in tag_indicator.iter().enumerate() {
        //     if b > &0.0 {
        //         hit_ind.push(i);
        //     }
        // }
        // return hit_ind;
    }
    /// get ES
    pub fn select_es(&self, run_es: &[f64]) -> f64 {
        let n = self.gene.size() as f64;
        let mut es: f64;
        if self.single {
            es = run_es.iter().sum::<f64>(); // FIXME: Overflow ?
            if self.scale {
                es = es / n;
            }
            return es;
        }

        // let max = run_es.iter().max().unwrap(); // max is not working for &f64
        let max = run_es.iter().fold(std::f64::MIN, |a, b| a.max(*b));
        let min = run_es.iter().fold(std::f64::MAX, |a, b| a.min(*b));
        // let argmax = run_es.iter().position(|r|  r== &max).unwrap();
        es = if max.abs() > min.abs() { max } else { min };
        return es;
    }
    /// phenotype permutation procedure
    /// shuffling group labels and calculate the new ranking metric
    /// return a vector of argsorted tuples (indices, sorted_metric)
    /// data - 2d vector [m_genes, n_samples]
    pub fn phenotype_permutation(
        &mut self,
        data: &[Vec<f64>],
        group: &[bool],
        method: Metric,
        ascending: bool,
    ) -> Vec<(Vec<usize>, Vec<f64>)> {
        let mut group_nperm = vec![group.to_vec(); self.nperm]; // Note: self.nperm has been ++
                                                                // let mut g = group.to_vec();
        for i in 1..self.nperm {
            //fastrand::shuffle(&mut group_arc);
            group_nperm[i].shuffle(&mut self.rng);
        }
        let arr: Vec<(Vec<usize>, Vec<f64>)> = group_nperm
            .par_iter()
            .map(|group_rng| {
                // implement methods in trait so you can call self.member_function in self closure
                // parallel computation inside
                let m = self.calculate_metric(data, &group_rng, method);
                m.as_slice().argsort(ascending) // default is false
            })
            .collect();
        return arr;
    }

    // /// phenotype permutation procedurce
    // /// return (hit_index, es, runing_es) only. not null distribution return
    pub fn enrichment_score_pheno(
        &self,
        sorted_metric: &[Vec<f64>],
        tag_indicator: &[Vec<f64>],
    ) -> (Vec<f64>, Vec<f64>) {
        let run_es: Vec<f64> = self.running_enrichment_score(&sorted_metric[0], &tag_indicator[0]);
        let ess: Vec<f64> = sorted_metric
            .par_iter()
            .zip(tag_indicator.par_iter())
            .map(|(gm, tag)| {
                let r = self.fast_random_walk(gm, tag);
                r
            })
            .collect();
        return (ess, run_es);
    }

    /// gene permutation percedure.
    /// shufling the genes in the given ranking metric
    /// use for prerank, ssgsea
    /// only kept record of gene_indexes of original order to accerrlate the shuffling
    pub fn gene_permutation(&mut self) -> Vec<DynamicEnum<usize>> {
        let vec: Vec<usize> = (0..self.gene.size()).collect();
        let mut orig: DynamicEnum<usize> = DynamicEnum::from(&vec);
        let mut gperm: Vec<DynamicEnum<usize>> = Vec::new();
        // // now = Instant::now();
        gperm.push(orig.clone());
        for _ in 1..self.nperm {
            // inplace shuffle
            //fastrand::shuffle(&mut tags[i]); // fastrand is a little bit faster
            orig.shuffle(&mut self.rng);
            gperm.push(orig.clone());
        }
        return gperm;
    }
    /// gene_set permutation procedurce.
    /// use vec slice as gene_set argument. since it accepts &vec<T>
    /// return (hit_index, enrichment_score, running_es) with null distrition
    pub fn enrichment_score_gene(
        &mut self,
        metric: &[f64],
        tag_indicators: &[Vec<f64>],
    ) -> (Vec<f64>, Vec<f64>) {
        // pararell computing
        // let run_es: Vec<Vec<f64>> = tag_indicators
        //     .par_iter()
        //     .map(|tag| {
        //         // implement the function in trait enable you to call self.methods in struct member closure!!!
        //         self.running_enrichment_score(metric, tag)
        //     })
        //     .collect();
        // let es: Vec<f64> = run_es.par_iter().map(|r| self.select_es(r)).collect();

        let es: Vec<f64> = tag_indicators
            .par_iter()
            .map(|tag| self.fast_random_walk(metric, tag))
            .collect();
        let run_es = self.running_enrichment_score(metric, &tag_indicators[0]);

        return (es, run_es);
    }
}

// mod tests {

//     use super::*;
//     use std::time::Instant;
//     use crate::utils::FileReader;
//     #[test]
//     fn test_calcu_metric () {
//         let mut gct = FileReader::new();
//         let _ = gct.read_csv("data/P53.txt", b'\t', true, Some(b'#'));
//         let mut cls = FileReader::new();
//         let _ = cls.read_table("data/P53.cls", ' ', false);
//         // let mut group = cls.record.pop().unwrap();
//         println!("{:?}", &cls.record[2]);
//         let gboo: Vec<bool> = cls.record[2].iter().map(|x| x != "WT").collect();
//         println!("{:?}", &gboo);
//         let weight = 1.0;
//         let mut genes: Vec<String> = Vec::new();
//         // let mut gene_set: Vec<String> = Vec::new();
//         let mut gene_exp: Vec<Vec<f64>> = Vec::new();
//         for r in gct.record.iter() {
//             genes.push(r[0].to_string());
//             let mut vv: Vec<f64> = Vec::new();
//             for v in &r[2..] {
//                 vv.push(v.parse::<f64>().unwrap());
//             }
//             vv.iter_mut().for_each(|x| *x = x.abs().powf(weight)); // modified value in place
//             gene_exp.push(vv);
//         }
//         let end1 = Instant::now();
//         let es = EnrichmentScore::new(&genes, 0, 123, false, false);
//         //let gene_metric = es.phenotype_permutation(&gene_exp, &gboo, Metric::Signal2Noise);
//         let _metric = es.calculate_metric(&gene_exp, &gboo, Metric::Signal2Noise);
//         println!("presorted metric checked: {:.4?}", &_metric[..10]);
//         let weighted_sorted_metric = _metric.as_slice().argsort(false);
//         let end2 = Instant::now();
//         println!("metric length: {}", &weighted_sorted_metric.1.len());
//         let idx = weighted_sorted_metric.0;
//         println!("sorted metric checked: ");
//         for (i, v) in (&weighted_sorted_metric.1[..10]).iter().enumerate()
//         {
//             println!("{:?}, {:.4?}", es.gene.elt_of(idx[i]), v);
//         }
//         println!("Overall run time: {:.2?}", end2.duration_since(end1));
//     }
// }
//    #[test]
//    fn gene_permutation_fastrand()
//    {
//        let nperm = 1000;
//        let rang = (0..20000).collect();
//        let mut tags: Vec<Vec<i64>> = vec![rang; 1+nperm]; // 2d array
//        let mut notags: Vec<Vec<i64>> = Vec::new();
//        let mut rng = fastrand::Rng::with_seed(43);
//        let start = Instant::now();
//        for i in 1..=nperm {
//            // inplace shuffle
//            //fastrand::shuffle(&mut tags[i]);
//            rng.shuffle(&mut tags[i]);
//            let ntg: Vec<i64>  = tags[i].iter().map( |b| 1- (*b)).collect();
//            notags.push(ntg);
//        }

//        let end = Instant::now();
//        println!("permutation time: {:.2?}", end.duration_since(start));
//        //tags.iter().for_each(|x| println!("{:?}", x))
//    }

//    #[test]
//    fn gene_permutation_smallrng()
//    {
//        let nperm = 1000;
//        let rang = (0..20000).collect();
//        let mut tags: Vec<Vec<i64>> = vec![rang; 1+nperm]; // 2d array
//        let mut notags: Vec<Vec<i64>> = Vec::new();
//        let mut rng = SmallRng::seed_from_u64(43);
//        let start = Instant::now();
//        for i in 1..=nperm {
//            // inplace shuffle
//            //fastrand::shuffle(&mut tags[i]);
//            tags[i].shuffle(&mut rng);
//            let ntg: Vec<i64>  = tags[i].iter().map( |b| 1- (*b)).collect();
//            notags.push(ntg);
//        }

//        let end = Instant::now();
//        println!("permutation time: {:.2?}", end.duration_since(start));
//        //tags.iter().for_each(|x| println!("{:?}", x))
//    }

//    #[test]
//    fn gene_permutation_xoshiro256plus()
//    {
//        let nperm = 1000;
//        let rang = (0..20000).collect();
//        let mut tags: Vec<Vec<i64>> = vec![rang; 1+nperm]; // 2d array
//        let mut notags: Vec<Vec<i64>> = Vec::new();
//        let mut rng = Xoshiro256Plus::seed_from_u64(0);
//        let start = Instant::now();
//        for i in 1..=nperm {
//            // inplace shuffle
//            //fastrand::shuffle(&mut tags[i]);
//            tags[i].shuffle(&mut rng);
//            let ntg: Vec<i64>  = tags[i].iter().map( |b| 1- (*b)).collect();
//            notags.push(ntg);
//        }

//        let end = Instant::now();
//        println!("permutation time: {:.2?}", end.duration_since(start));
//        //tags.iter().for_each(|x| println!("{:?}", x))
//    }
//    #[test]
//    fn gene_permutation_mcg128xsl6()
//    {
//        let nperm = 1000;
//        let rang = (0..20000).collect();
//        let mut tags: Vec<Vec<i64>> = vec![rang; 1+nperm]; // 2d array
//        let mut notags: Vec<Vec<i64>> = Vec::new();
//        let mut rng = Mcg128Xsl64::seed_from_u64(0);
//        let start = Instant::now();
//        for i in 1..=nperm {
//            // inplace shuffle
//            //fastrand::shuffle(&mut tags[i]);
//            tags[i].shuffle(&mut rng);
//            let ntg: Vec<i64>  = tags[i].iter().map( |b| 1- (*b)).collect();
//            notags.push(ntg);
//        }

//        let end = Instant::now();
//        println!("permutation time: {:.2?}", end.duration_since(start));
//        //tags.iter().for_each(|x| println!("{:?}", x))
//   }

#[cfg(test)]
mod tests {
    use super::*;

    /// Reconstruct the full O(N) running-ES vector from the O(k) tops/bottoms arrays
    /// produced by `calc_gsea_stat`.
    ///
    /// The enrichment-score curve is a staircase:
    /// - At a hit position `selected[i]`:          run_es = tops[i]
    /// - At any miss position p after hit i−1:     run_es = tops[i-1] − (p − selected[i-1]) * inv_miss
    /// - Before the first hit (p < selected[0]):   run_es = −(p + 1) * inv_miss
    ///
    /// This is the inverse of the O(k) compression: `bottoms[i]` should equal the
    /// reconstructed value at `selected[i] − 1` (the last miss before hit i).
    fn reconstruct_run_es_from_tops(tops: &[f64], selected: &[usize], n: usize) -> Vec<f64> {
        let k = selected.len();
        let inv_miss = 1.0 / (n - k) as f64;
        let mut run_es = vec![0.0f64; n];
        let mut hit_i = 0usize;
        // last_top and last_hit_pos track the ES and position of the previous hit.
        // We use i64 for last_hit_pos so we can represent "before any hit" as -1.
        let mut last_top = 0.0f64;
        let mut last_hit_pos: i64 = -1;

        for p in 0..n {
            if hit_i < k && p == selected[hit_i] {
                run_es[p] = tops[hit_i];
                last_top = tops[hit_i];
                last_hit_pos = p as i64;
                hit_i += 1;
            } else {
                let misses = p as i64 - last_hit_pos;
                run_es[p] = last_top - misses as f64 * inv_miss;
            }
        }
        run_es
    }

    /// Verify that `calc_gsea_stat` (O(k)) produces the same ES as `fast_random_walk`
    /// for the standard score type with gsea_param = 1.0.
    ///
    /// `fast_random_walk` takes a pre-weighted tag_indicator array and computes max|ES|,
    /// which equals `calc_gsea_stat` with ScoreType::Std and gsea_param = 1.0.
    #[test]
    fn test_calc_gsea_stat_matches_fast_random_walk() {
        // stats sorted descending, N = 10, gene set at positions [0, 2, 5, 8]
        let stats: Vec<f64> = vec![10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.5, 1.0, 0.5];
        let selected: Vec<usize> = vec![0, 2, 5, 8];
        let gsea_param = 1.0;

        // Build tag_indicator as expected by fast_random_walk (already weighted by gsea_param=1)
        let genes: Vec<String> = (0..stats.len()).map(|i| format!("g{}", i)).collect();
        let es_obj = EnrichmentScore::new(&genes, 0, 0, false, false);
        let mut tag: Vec<f64> = vec![0.0; stats.len()];
        for &i in &selected {
            tag[i] = stats[i]; // fast_random_walk uses the metric * tag directly
        }
        let frw_es = es_obj.fast_random_walk(&stats, &{
            let mut t = vec![0.0f64; stats.len()];
            for &i in &selected {
                t[i] = 1.0;
            }
            t
        });

        let result = calc_gsea_stat(&stats, &selected, gsea_param, ScoreType::Std, true, true);

        assert!(
            (result.es - frw_es).abs() < 1e-12,
            "ES mismatch: calc_gsea_stat={}, fast_random_walk={}",
            result.es,
            frw_es
        );

        // tops and bottoms must be populated
        assert_eq!(result.tops.len(), selected.len());
        assert_eq!(result.bottoms.len(), selected.len());

        // tops[i] >= bottoms[i] always (hit step increases ES)
        for (t, b) in result.tops.iter().zip(result.bottoms.iter()) {
            assert!(t >= b, "tops[i] < bottoms[i]: {} < {}", t, b);
        }

        // ES must equal the maximum of all tops/bottoms extremes (std mode)
        let max_p = result.tops.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_p = result.bottoms.iter().cloned().fold(f64::INFINITY, f64::min);
        let expected_es = if max_p.abs() >= min_p.abs() { max_p } else { min_p };
        assert!(
            (result.es - expected_es).abs() < 1e-12,
            "ES does not match max|tops/bottoms|: es={}, expected={}",
            result.es,
            expected_es
        );
    }

    /// Verify tops/bottoms produce the correct running ES staircase shape.
    /// For a simple case with uniform weights, the curve should be symmetric.
    #[test]
    fn test_calc_gsea_stat_extremes_shape() {
        // N=6, gene set at positions [0, 3] (2 hits, 4 misses)
        // uniform weights (gsea_param=0 makes all |stats|^0 = 1)
        let stats: Vec<f64> = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.5];
        let selected: Vec<usize> = vec![0, 3];
        let result = calc_gsea_stat(&stats, &selected, 0.0, ScoreType::Std, true, false);

        // With uniform weights: NR = 2*1 = 2 (both weights = 1^0 = 1)
        // r_adj_norm = [0.5, 0.5]
        // tops[0]    = 0.5 - (0-0)/4 = 0.5
        // tops[1]    = 1.0 - (3-1)/4 = 1.0 - 0.5 = 0.5
        // bottoms[0] = 0.5 - 0.5 = 0.0
        // bottoms[1] = 0.5 - 0.5 = 0.0
        assert!((result.tops[0] - 0.5).abs() < 1e-12);
        assert!((result.tops[1] - 0.5).abs() < 1e-12);
        assert!((result.bottoms[0] - 0.0).abs() < 1e-12);
        assert!((result.bottoms[1] - 0.0).abs() < 1e-12);
        // ES = 0.5 (or 0.0, maxP == -minP → 0; but here min_p = 0.0 and max_p = 0.5 → 0.5)
        assert!((result.es - 0.5).abs() < 1e-12);
    }

    /// Compare the full O(N) running-ES vector from `running_enrichment_score` with the
    /// vector reconstructed from `calc_gsea_stat`'s tops/bottoms arrays.
    ///
    /// The two representations must agree at every position to within floating-point
    /// precision (ε = 1e-12).  We also explicitly verify that `bottoms[i]` equals the
    /// reconstructed running-ES at position `selected[i] - 1` (the last miss before
    /// each hit), ensuring the staircase knots are correctly placed.
    #[test]
    fn test_run_es_full_vector_matches_reconstruction() {
        // N = 12, k = 4 hits at positions [1, 4, 7, 10]; 8 misses.
        // stats are sorted descending and all positive (gsea_param = 1).
        let stats: Vec<f64> = vec![
            12.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.5, 1.5, 0.5,
        ];
        let selected: Vec<usize> = vec![1, 4, 7, 10];
        let n = stats.len();
        let gsea_param = 1.0;

        // --- O(N) reference: running_enrichment_score ---
        let genes: Vec<String> = (0..n).map(|i| format!("g{}", i)).collect();
        let es_obj = EnrichmentScore::new(&genes, 0, 0, false, false);
        let mut tag: Vec<f64> = vec![0.0f64; n];
        for &i in &selected {
            tag[i] = 1.0;
        }
        let ref_run_es = es_obj.running_enrichment_score(&stats, &tag);

        // --- O(k) compressed form ---
        let result = calc_gsea_stat(&stats, &selected, gsea_param, ScoreType::Std, true, false);

        // --- Reconstruct full O(N) vector from tops ---
        let recon_run_es = reconstruct_run_es_from_tops(&result.tops, &selected, n);

        // Element-wise comparison
        assert_eq!(
            ref_run_es.len(),
            recon_run_es.len(),
            "vector length mismatch"
        );
        for (p, (&r, &c)) in ref_run_es.iter().zip(recon_run_es.iter()).enumerate() {
            assert!(
                (r - c).abs() < 1e-12,
                "run_es mismatch at position {}: running_enrichment_score={:.15}, \
                 reconstructed={:.15}, diff={:.3e}",
                p,
                r,
                c,
                (r - c).abs()
            );
        }

        // Explicit spot checks for bottoms:
        // bottoms[i] = run_es just before hit i, i.e. at selected[i] - 1
        for (i, &hit_pos) in selected.iter().enumerate() {
            if hit_pos > 0 {
                let before_hit = ref_run_es[hit_pos - 1];
                assert!(
                    (result.bottoms[i] - before_hit).abs() < 1e-12,
                    "bottoms[{}] mismatch: bottoms={:.15}, run_es[{}]={:.15}",
                    i,
                    result.bottoms[i],
                    hit_pos - 1,
                    before_hit
                );
            }
        }

        // Tops must equal run_es at hit positions
        for (i, &hit_pos) in selected.iter().enumerate() {
            assert!(
                (result.tops[i] - ref_run_es[hit_pos]).abs() < 1e-12,
                "tops[{}] mismatch: tops={:.15}, run_es[{}]={:.15}",
                i,
                result.tops[i],
                hit_pos,
                ref_run_es[hit_pos]
            );
        }
    }

    /// Verify the full-vector reconstruction with a larger, randomised example
    /// to guard against edge cases (hits at first/last position, dense gene sets).
    #[test]
    fn test_run_es_full_vector_large() {
        // N = 20, k = 6 hits at positions spread across the ranking.
        // Include hits at position 0 and 19 to test boundary conditions.
        let stats: Vec<f64> = (0..20u32)
            .rev()
            .map(|i| i as f64 + 1.0)
            .collect(); // [20,19,...,1]
        let selected: Vec<usize> = vec![0, 3, 7, 11, 15, 19];
        let n = stats.len();

        let genes: Vec<String> = (0..n).map(|i| format!("g{}", i)).collect();
        let es_obj = EnrichmentScore::new(&genes, 0, 0, false, false);
        let mut tag = vec![0.0f64; n];
        for &i in &selected {
            tag[i] = 1.0;
        }
        let ref_run_es = es_obj.running_enrichment_score(&stats, &tag);

        let result = calc_gsea_stat(&stats, &selected, 1.0, ScoreType::Std, true, false);
        let recon = reconstruct_run_es_from_tops(&result.tops, &selected, n);

        for (p, (&r, &c)) in ref_run_es.iter().zip(recon.iter()).enumerate() {
            assert!(
                (r - c).abs() < 1e-12,
                "large test: mismatch at p={}: ref={:.15} recon={:.15}",
                p,
                r,
                c
            );
        }
    }
}
