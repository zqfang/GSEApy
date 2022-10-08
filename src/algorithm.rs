#![allow(dead_code, unused)]

use crate::utils::DynamicEnum;
use crate::utils::{Metric, Statistic};
use rand::rngs::SmallRng; // use SmallRng intestad of StdRng to speedup shuffling
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;

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
                    Metric::Signal2Noise => (pos_mean - neg_mean) / (pos_std + neg_std),
                    Metric::AbsSignal2Noise => ((pos_mean - neg_mean) / (pos_std + neg_std)).abs(),
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
        // let rng = ThreadRng::default();
        let rng = SmallRng::seed_from_u64(seed);
        //let rng = thread_rng();
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
