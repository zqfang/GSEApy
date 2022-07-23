use serde::{Deserialize, Serialize};
use crate::algorithm::{EnrichmentScore, EnrichmentScoreTrait};
use crate::utils::{Metric, Statistic};
use pyo3::prelude::*;
use rayon::prelude::*;
use std::time::Instant;
use std::collections::HashMap;
use itertools::izip;
// serialize struct example
// let point = Point { x: 1, y: 2 };

// let serialized = serde_json::to_string(&point).unwrap();
// println!("serialized = {}", serialized);

// let deserialized: Point = serde_json::from_str(&serialized).unwrap();
// println!("deserialized = {:?}", deserialized);

#[pyclass]
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
    pub run_es: Vec<f64>,
    #[pyo3(get, set)]
    pub hits: Vec<usize>, // indices of genes that matches
    #[pyo3(get, set)]
    pub esnull: Vec<f64>,
    #[pyo3(get, set)]
    pub name: Option<String>,
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
        name: &str,
    ) -> Self {
        GSEASummary {
            term: term.to_string(),
            es: es,
            nes: nes,
            pval: pval,
            fwerp: fwerpval,
            fdr: fdr,
            run_es: run_es.to_vec(),
            hits: hits.to_vec(),
            esnull: esnull.to_vec(),
            name: Some(name.to_string())
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
            run_es: Vec::<f64>::new(),
            hits: Vec::<usize>::new(),
            esnull: Vec::<f64>::new(),
            name: None,
        }
    }

    fn normalize(&mut self) -> Vec<f64> {
        let n_mean: f64;
        let e: f64 = self.es;
        let tmp: Vec<f64>; 
        if e >= 0.0 {
            // n_mean = esnull[esnull>= 0].mean()
            tmp = self
                .esnull
                .iter()
                .filter_map(|&x| {if x >= 0.0 {Some(x)} else {None}})
                .collect();
        } else {
            // n_mean = esnull[esnull< 0].mean()
            tmp = self
                .esnull
                .iter()
                .filter_map(|&x| {if x < 0.0 {Some(x)} else {None}})
                .collect();    
        }

        // FIXME: Potential NaN number here
        // When input a rare causes of an extreamly screwed null distribution. e.g.
        // es = - 27, esnull = [13, 24, 57, 88]
        // nes will be NaN. You have to increased the permutation number for safe
        // a tricky fixed here: set n_mean as itself
        // so esnull = [-27, 13, 24, 57, 88]
        n_mean = if tmp.len() > 0 { tmp.iter().sum::<f64>() / (tmp.len() as f64) } else { e };
        self.nes = if e >= 0.0 {e / n_mean} else {-1.0 * e / n_mean};

        let nesnull: Vec<f64> = self
            .esnull
            .iter()
            .map(|&e| {
                if e >= 0.0 {
                    e / n_mean
                } else {
                    -1.0 * e / n_mean
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
            deno =  self.esnull.iter().filter(|&x| *x < 0.0).count();
            nomi = self.esnull.iter().filter(|&x| x < &self.es).count();
        } else {
            deno = self.esnull.iter().filter(|&x| *x >= 0.0).count();
            nomi = self.esnull.iter().filter(|&x| x >= &self.es).count();
        }

        if deno == 0 { 
            self.pval = 1.0; 
            return 
        }
        self.pval = (nomi as f64) / (deno as f64);
    }
}

#[derive(Debug, Clone, Default)]
pub struct GSEAResult {
    pub summaries: Vec<GSEASummary>,
    weight: f64,
    min_size: usize,
    max_size: usize,
    nperm: usize,
    nes_concat: Vec<f64>,
    nesnull_concat: Vec<f64>,
    seed: u64,
}

impl GSEAResult {
    pub fn new(
        weight: f64,
        max_size: usize,
        min_size: usize,
        nperm: usize,
        seed: u64,
    ) -> Self {
        GSEAResult {
            summaries: Vec::<GSEASummary>::new(),
            weight: weight,
            max_size: max_size,
            min_size: min_size,
            nperm: nperm,
            nes_concat: Vec::<f64>::new(),
            nesnull_concat: Vec::<f64>::new(),
            seed: seed,
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
        }
    }
    pub fn stat(&mut self, summary: &mut [GSEASummary] ) {
        // clear vector incase you re-run this command
        self.nes_concat.clear();
        self.nesnull_concat.clear();

        summary.iter_mut().for_each(|g| {
            // calculate stats here
            g.pval();
            let mut nesnull = g.normalize(); // update esnull to normalized nesnull
            self.nes_concat.push(g.nes);
            self.nesnull_concat.append(&mut nesnull);
            // g.esnull.clear();
        });
        // FWER p 
        let fwerps : Vec<f64> = self.fwer_pval();
        // FDR q
        let fdrs = self.fdr();


        for (p, q, g) in izip!(fwerps, fdrs, summary) {
            g.fdr = q;
            g.fwerp = p;
         }
        // clear vector to save some space
        self.nes_concat.clear();
        self.nesnull_concat.clear();
    }


    pub fn fdr(&mut self) -> Vec<f64> {
        // let mut nesnull_concat: Vec<&f64> = nesnull.iter().flatten().collect(); // nesnull.concat(); // concat items

        // To speedup, sort f64 in acending order in place, then do a binary search
        self.nesnull_concat
            .sort_unstable_by(|a, b| a.partial_cmp(b).unwrap()); // descending -> b.partial_cmp(a)        
        let (indices, nes_sorted) = self.nes_concat.as_slice().argsort(true); // ascending order

        // binary_search assumes that the elements are sorted in less-to-greater order.
        // partition_point return the index of the first element of the second partition)
        // since partition_point is just a wrapper of self.binary_search_by(|x| if pred(x) { Less } else { Greater }).unwrap_or_else(|i| i)
        let all_idx = self.nesnull_concat.partition_point(|x| *x < 0.0);
        let nes_idx = nes_sorted.partition_point(|x| *x < 0.0);

        // fdr
        let mut fdrs: Vec<f64> = nes_sorted
            .iter()
            .map(|&e| {
                let pi_norm: f64;
                let pi_obs: f64;
                let nes_higher: usize;
                let all_higher: usize;
                let all_pos: usize;
                let nes_pos: usize;
                if e < 0.0 {
                    // let nes_higher = nes_concat.iter().filter(|&x| *x < e).count();
                    // let all_higher = nesnull_concat.iter().filter(|&x| *x < e).count();
                    nes_higher = nes_sorted.partition_point(|x| *x < e); // left side
                    all_higher = self.nesnull_concat.partition_point(|x| *x < e); // left side
                    all_pos = all_idx;
                    nes_pos = nes_idx;
                } else {
                    // let nes_higher = self.nes_concat.iter().filter(|&x| *x >= e).count();
                    // let all_higher = self.nesnull_concat.iter().filter(|&x| *x >= e).count();
                    nes_higher = nes_sorted.len() - nes_sorted.partition_point(|x| *x < e); // right side
                    all_higher =
                        self.nesnull_concat.len() - self.nesnull_concat.partition_point(|x| *x < e); // right side
                    all_pos = self.nesnull_concat.len() - all_idx; // right side
                    nes_pos = nes_sorted.len() - nes_idx; // right side
                }
                // println!("neg_higher {}, all_higher {}, all_pos {}, nes_pos {}", nes_higher, all_higher, all_pos, all_higher);
                pi_norm = (all_higher as f64) / (all_pos as f64);
                pi_obs = (nes_higher as f64) / (nes_pos as f64);
                // FDR
                (pi_norm / pi_obs).clamp(f64::MIN, 1.0)
            })
            .collect();

        // adjusted fdr q value
        self.adjust_fdr(&mut fdrs, nes_idx);
        let mut fdr_orig_order: Vec<f64> = vec![0.0; fdrs.len()];
        indices.iter().zip(fdrs.iter()).for_each(|(&i, &v)| {fdr_orig_order[i] = v;} );
        return fdr_orig_order;
    }

    /// # adjust fdr q-values
    /// see line 880:  https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    /// - fdrs:  Corresponds to the ascending order of NES.
    /// - partition_point_idx: the index of the first element of the second partition
    /// This function updates fdr value inplace.
    fn adjust_fdr(&self, fdrs: &mut [f64], partition_point_idx: usize) 
    {
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
        if partition_point_idx < fdrs.len()  
        {
            // pos_nes on the right side, if only have postive numbers, idx must be < .len() 
            let nes_pos_idx =  partition_point_idx + 1;
            min_fdr = fdrs[partition_point_idx];
            for k in nes_pos_idx..fdrs.len()
            { 
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

        if partition_point_idx > 0
        {
            // neg_nes on the left side, if only have negative numbers, idx must be > 0
            let nes_neg_idx = partition_point_idx - 1;
            min_fdr = fdrs[nes_neg_idx];
            for k in (0..partition_point_idx).rev()
            { 
                min_fdr = min_fdr.min(fdrs[k]);
                fdrs[k] = min_fdr.min(fdrs[k]);
            }
        }
        
    }
    /// Compute FWER p-vals
    /// line 788: https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
    fn fwer_pval(&self) -> Vec<f64>
    {
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
        
        let fwerp: Vec<f64> = self.nes_concat.par_iter().map(|e| { 
            if e < &0.0 {
                (min_nes_neg.iter().filter(|&x| x < e ).count() as f64)
                    / (min_nes_neg.iter().filter(|&x| x < &0.0).count() as f64)
            } else {
                (max_nes_pos.iter().filter(|&x| x >= e).count() as f64)
                    / (max_nes_pos.len() as f64)
            }
        }).collect();
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
        gmt: &HashMap<&str, &[String]>,
        method: Metric,
    ) {
        let mut es = EnrichmentScore::new(genes, self.nperm, self.seed, false, false);
        // let end = Instant::now();
        let sorted_metric: Vec<(Vec<usize>, Vec<f64>)> = es.phenotype_permutation(gene_exp, group, method, false);
        // let end1 = Instant::now();
        // println!("Permutation time: {:.2?}", end1.duration_since(end));

        let mut summ = Vec::<GSEASummary>::new();
        for (&term, &gset) in gmt.iter() {
            let tag = es.gene.isin(gset);
            // get es hit index of sorted array
            let tag_new: Vec<f64> = sorted_metric[0].0.iter().map(|&i| tag[i]).collect();
            let gidx = es.hit_index(&tag_new); // need updat the sorted indices
            if gidx.len() > self.max_size || gidx.len() < self.min_size {
                continue;
            }
            // get running enrichment score
            let run_es: Vec<Vec<f64>> = sorted_metric
            .par_iter()
            .map(|(indices, gm)| {
                 // weight the metrics
                let weighted_gm: Vec<f64> = gm.iter().map(| x| {x.abs().powf(self.weight)}).collect();
                // update tag_indicator since you've update metric
                let tag_new: Vec<f64> = indices.iter().map(|&i| tag[i]).collect();
                // calculate ES
                let r = es.running_enrichment_score(&weighted_gm, &tag_new);
                r
            })
            .collect();

            // get es 
            let ess: Vec<f64> = run_es.par_iter().map(|r| es.select_es(r)).collect();
            // let (ess, run_es) = es.enrichment_score_pheno(&weighted_metric, &tag);

            let gss = GSEASummary {
                term: term.to_string(), 
                es: ess[0],
                run_es: run_es[0].to_owned(),
                hits: gidx,
                esnull: ess[1..].to_owned(),
                ..Default::default()
            };
            summ.push(gss);
        }

        // let end2 = Instant::now();
        // println!("Calculation time: {:.2?}", end2.duration_since(end1));
        self.stat(&mut summ);
        self.summaries = summ;

        // let end3 = Instant::now();
        // println!("Statistical time: {:.2?}", end3.duration_since(end2));

        // update es.gene attribute since we've kept track of sorted metric indices
        // let mut sorted_genes =  DynamicEnum::<String>::new();
        //  sorted_metric[0].0.iter().for_each(|&i| {
        //     let g = es.gene.elt_of(i).unwrap().to_string();
        //     sorted_genes.add_if_new(g);
        // });
        // es.gene = sorted_genes;

    }

    pub fn prerank(&mut self, 
        genes: &[String], 
        metric: &[f64], 
        gmt: &HashMap<&str, &[String]>) {
        // NOTE: input must not contain duplcated genes

        let weighted_metric: Vec<f64> = metric.iter().map(|x|  x.abs().powf(self.weight)).collect();
        // start to calculate
        let mut es = EnrichmentScore::new(genes, self.nperm, self.seed, false, false);
        // let end1 = Instant::now();
        let gperm = es.gene_permutation(); // gene permutation, only record gene idx here
        // let end2 = Instant::now();
        // println!("Permutation time: {:.2?}", end2.duration_since(end1));
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
            let gss = GSEASummary {
                term: term.to_string(),
                es: ess[0],
                run_es: run_es,
                hits: gidx,
                esnull: ess[1..].to_owned(),
                ..Default::default()
            };
            summ.push(gss);
        }
        // let end3 = Instant::now();
        // println!("Calculation time: {:.2?}", end3.duration_since(end2));

        self.stat(&mut summ);
        self.summaries = summ;
        // let end4 = Instant::now();
        // println!("Statistical time: {:.2?}", end4.duration_since(end3));
    }
    /// single sample gsea
    pub fn ss_gsea(
        &mut self,
        genes: &[String],
        samples: &[String],
        gene_exp: &[Vec<f64>], // 2d vector [m_gene, n_sample];
        gmt: &HashMap<&str, &[String]>,

    ) {
        // transpose [m_gene, n_sample] --> [n_sample, m_gene]
        let mut gene_metric: Vec<Vec<f64>> = vec![vec![]; samples.len()];
        gene_exp.iter().for_each(|row| {
            row.iter().enumerate().for_each(| (j, e)| {gene_metric[j].push(*e);});
        });

        // sort first and then set weight, 
        let weighted_sorted_metric: Vec<(Vec<usize>, Vec<f64>)> = gene_metric.into_par_iter().map(|rank| {
            let mut tmp = rank.as_slice().argsort(false);
            tmp.1.iter_mut().for_each(|x|  {*x = x.abs().powf(self.weight);});
            return tmp;
        }).collect();
        
        let es = EnrichmentScore::new(genes, self.nperm, self.seed, true, false);
        // let end1 = Instant::now();
        for (&term, &gset) in gmt.iter() {
            let tag = es.gene.isin(gset);
            let hit = tag.iter().filter(|&x| x > &0.0).count();
            if hit > self.max_size || hit < self.min_size {
                continue;
            }
            let mut summ: Vec<GSEASummary> = weighted_sorted_metric.par_iter().enumerate().map(|(i, (indices, metric))| 
            {
                let tag_new: Vec<f64> = indices.iter().map(|&idx| tag[idx]).collect();
                let gidx = es.hit_index(&tag_new);
                let run_es = es.running_enrichment_score(metric, &tag_new);
                let es = es.select_es(&run_es);
                GSEASummary {
                    term: term.to_string(),
                    es: es,
                    run_es: run_es,
                    hits: gidx,  // gene hit idx of each sample after sorting
                    name: Some(samples[i].to_string()),
                    ..Default::default()
                }
            }).collect();
            self.summaries.append(&mut summ);
        }
        // let end2 = Instant::now();
        // println!("Calculation time: {:.2?}", end2.duration_since(end1));
        // self.stat(); // NES
        let max = self.summaries.iter().fold(std::f64::MIN, |a, b| a.max(b.es));
        let min = self.summaries.iter().fold(std::f64::MAX, |a, b| a.min(b.es));
        let norm = max - min;
        self.summaries.iter_mut().for_each(|b| b.nes = b.es / norm);

        // let end3 = Instant::now();
        // println!("Statistical time: {:.2?}", end3.duration_since(end2));
    }



    /// single sample gsea
    pub fn ss_gsea_permuate(
        &mut self,
        genes: &[String],
        samples: &[String],
        gene_exp: &[Vec<f64>], // 2d vector [m_gene, n_sample];
        gmt: &HashMap<&str, &[String]>,

    ) {
        // transpose [m_gene, n_sample] --> [n_sample, m_gene]
        let mut gene_metric: Vec<Vec<f64>> = vec![vec![]; samples.len()];
        gene_exp.iter().for_each(|row| {
            row.iter().enumerate().for_each(| (j, e)| {gene_metric[j].push(*e);});
        });

        // sort first and then set weight, 
        let weighted_sorted_metric: Vec<(Vec<usize>, Vec<f64>)> = gene_metric.into_par_iter().map(|rank| {
            let mut tmp = rank.as_slice().argsort(false);
            tmp.1.iter_mut().for_each(|x|  {*x = x.abs().powf(self.weight);});
            return tmp;
        }).collect();

        let mut _all = Vec::<GSEASummary>::new();
        // let end1 = Instant::now();
        weighted_sorted_metric.into_iter().enumerate().for_each(|(i, (indices, metric))| 
        {
            // just run prerank procedure here
            // update the order of genes
            let _genes: Vec<String> = indices.into_iter().map(|j| {genes[j].to_string() }).collect(); 
            // init new ES object
            let mut es = EnrichmentScore::new(&_genes, self.nperm, self.seed, true, false);
            let gperm = es.gene_permutation(); // gene permutation, only record gene idx here
            let mut summ =  Vec::<GSEASummary>::new(); 
            for (&term, &gset) in gmt.iter() {

                let gtag = es.gene.isin(gset);
                let gidx = es.hit_index(&gtag);
                let tag_indicators: Vec<Vec<f64>> = gperm.par_iter().map(|de| de.isin(&gidx)).collect();
                if gidx.len() > self.max_size || gidx.len() < self.min_size { continue; }
                let (ess, run_es) = es.enrichment_score_gene(&metric, &tag_indicators);
                let gsu = GSEASummary {
                    term: term.to_string(),
                    es: ess[0],
                    run_es: run_es,
                    hits: gidx,  // hit index of each sample after sorting
                    esnull: ess[1..].to_owned(),
                    name: Some(samples[i].to_string()),
                    ..Default::default()
                };
                summ.push(gsu);
            }
            // calculate nes, pval, fdr
            self.stat(&mut summ);
            _all.append(&mut summ);
        });

        self.summaries = _all;
        // let end2 = Instant::now();
        // println!("Calculation time: {:.2?}", end2.duration_since(end1));
        // self.stat(); // NES
        //let max = self.summaries.iter().fold(std::f64::MIN, |a, b| a.max(b.es));
        //let min = self.summaries.iter().fold(std::f64::MAX, |a, b| a.min(b.es));
        //let norm = max - min;
        //self.summaries.iter_mut().for_each(|b| b.nes = b.es / norm);

        // let end3 = Instant::now();
        // println!("Statistical time: {:.2?}", end3.duration_since(end2));
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
        let start = Instant::now();
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()
            .unwrap();
        let mut rnk = FileReader::new();
        let _ = rnk.read_csv("data/temp.rnk", b'\t', false, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table("data/genes.gmt", '\t', false);

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
        let mut gmt2 = HashMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {gmt2.insert(r[0].as_str(), &r[2..]);});
        

        // weighted then sort
        gene_metric
            .iter_mut()
            .for_each(|x| *x = x.abs().powf(weight));
        let (gidx, metric) = gene_metric.as_slice().argsort(false);
        gene = gidx.iter().map(|&i| gene[i].clone()).collect();
        // start to calculate
        let mut gsea = GSEAResult::new(weight, 500, 3, 10, 123);
        gsea.prerank(&gene, &metric, &gmt2);
        let end = Instant::now();
        println!("Overall run time: {:.2?}", end.duration_since(start));

        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });
    }
    #[test]
    fn test_gsea() {
        let start = Instant::now();
        // set number of threads of rayon at the main()
        // rayon::ThreadPoolBuilder::new()
        //     .num_threads(1)
        //     .build_global()
        //     .unwrap();

        let mut gct = FileReader::new();
        let _ = gct.read_csv("data/P53.txt", b'\t', true, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table("data/genes.gmt", '\t', false);
        let mut cls = FileReader::new();
        let _ = cls.read_table("data/P53.cls", ' ', false);
        println!("{:?}", &cls.record[2]);
        let gboo: Vec<bool> = cls.record[2].iter().map(|x| x != "WT").collect();
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

        let mut gmt2 = HashMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {gmt2.insert(r[0].as_str(), &r[2..]);});

        let mut gsea = GSEAResult::new(weight, 1000, 3, 10, 123);
        gsea.gsea(&gene, &gboo, &gene_exp, &gmt2, Metric::Signal2Noise);

        let end = Instant::now();
        println!("Overall run time: {:.2?}", end.duration_since(start));
        gsea.summaries.iter().for_each(|g| {
            println!(
                "term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.term, g.es, g.nes, g.pval, g.fdr
            );
        });

    }

    #[test]
    fn test_ssgsea()
    {
        let mut gct = FileReader::new();
        let _ = gct.read_csv("data/P53.txt", b'\t', true, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table("data/genes.gmt", '\t', false);
        let mut cls = FileReader::new();
        let _ = cls.read_table("data/P53.cls", ' ', false);
        println!("{:?}", &cls.record[2]);
        let gboo: Vec<bool> = cls.record[2].iter().map(|x| x != "WT").collect();
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

        let mut gmt2 = HashMap::<&str, &[String]>::new();
        gmt.record.iter().for_each(|r| {gmt2.insert(r[0].as_str(), &r[2..]);});


        let nperm = 10;
        let mut gsea = GSEAResult::new(weight, 500, 3, nperm, 123);
        if nperm > 0 {
            gsea.ss_gsea_permuate(&gene, &sample_names, &gene_exp, &gmt2);
        } 
        else 
        {
            gsea.ss_gsea(&gene, &sample_names, &gene_exp, &gmt2);
        }

        gsea.summaries.iter().for_each(|g| {
            println!(
                "sample: {:?}, term: {:?}, es: {:.7?}, nes: {:.7?}, pval: {:.2e}, fdr: {:.2e}",
                g.name, g.term, g.es, g.nes, g.pval, g.fdr
            );
        });

    }
}
