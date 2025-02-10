#![allow(dead_code, unused)]

use crate::stats::{GSEAResult, GSEASummary};
use crate::utils::{DynamicEnum, Statistic};
use rayon::prelude::*;
use statrs::distribution::{ContinuousCDF, DiscreteCDF, Normal, Poisson};
use std::collections::HashMap;

pub struct GSVA {
    genes: DynamicEnum<String>,
    kcdf: bool,
    tau: f64,
    mx_diff: bool,
    abs_rnk: bool,
    rnaseq: bool,
    pre_res: usize,
    max_pre: usize,
    sigma: f64,
}

impl GSVA {
    pub fn default() -> Self {
        GSVA {
            genes: DynamicEnum::new(),
            kcdf: true,
            tau: 1.0,
            mx_diff: true,
            abs_rnk: false,
            rnaseq: false,
            pre_res: 10000,
            max_pre: 10,
            sigma: 4.0,
        }
    }
    pub fn new(
        genes: &[String],
        kcdf: bool,
        tau: f64,
        mx_dff: bool,
        abs_rnk: bool,
        rnaseq: bool,
    ) -> Self {
        GSVA {
            genes: DynamicEnum::from(genes),
            kcdf: kcdf,
            tau: tau,
            mx_diff: mx_dff,
            abs_rnk: abs_rnk,
            rnaseq: rnaseq,
            pre_res: 10000,
            max_pre: 10,
            sigma: 4.0,
        }
    }
    /// default: precomputed_resolution 10000, max_precompute: 10
    /// default: precomputed_resolution 10000, max_precompute: 10
    fn init_cdfs(&self, pre_res: usize, max_pre: usize) -> Vec<f64> {
        /// https://github.com/rcastelo/GSVA/blob/devel/src/kernel_estimation.c, line 123
        let norm = Normal::new(0.0, 1.0).unwrap();
        let divisor: f64 = pre_res as f64;
        /* standard normal distribution function, lower.tail=TRUE, log.p=FALSE */
        let res = (0..=pre_res)
            .into_iter()
            .map(|i| norm.cdf(((i * max_pre) as f64) / divisor))
            .collect();
        // println!("{:?}", &res);
        return res;
    }

    /// sigma factor: 4.0
    fn precomputed_cdfs(&self, x: f64, sigma: f64, pre_cdf: &[f64]) -> f64 {
        let v = x / sigma;
        let max_precompute = self.max_pre as f64;
        let precompute_resolution = self.pre_res as f64;
        if v < -1.0 * max_precompute {
            return 0.0;
        } else if v > max_precompute {
            return 1.0;
        } else {
            let idx: usize = (v.abs() / max_precompute * precompute_resolution) as usize;
            let cdf: f64 = pre_cdf[idx];
            if v < 0.0 {
                return 1.0 - cdf;
            } else {
                return cdf;
            }
        }
    }
    /// row: input are gene values of all samples.
    /// The empirical CDF is usually formally defined as
    /// CDF(x) = "number of samples <= x" / "number of samples"
    fn apply_ecdf(&self, row: &[f64]) -> Vec<f64> {
        let mut x0 = row.to_vec();
        let n = row.len() as f64;
        // To speedup, sort f64 in acending order in place, then do a binary search
        x0.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap()); // if descending -> b.partial_cmp(a)
                                                               // binary_search assumes that the elements are sorted in less-to-greater order.
                                                               // partition_point return the index of the first element of the second partition)
                                                               // since partition_point is just a wrapper of self.binary_search_by(|x| if pred(x) { Less } else { Greater }).unwrap_or_else(|i| i)
        row.iter()
            .map(|v| ((x0.partition_point(|x| x <= v)) as f64) / n)
            .map(|v| (v / (1.0 - v)).ln()) // https://github.com/rcastelo/GSVA/blob/devel/R/gsva.R, line 820
            .collect()
    }

    /// apply ecdf on th columns (genes)
    /// mat [n_samples, n_genes]
    fn mat_ecdf(&self, mat: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let mat3: Vec<Vec<f64>> = mat.into_iter().map(|v| self.apply_ecdf(v)).collect();
        return mat3;
    }

    fn row_d(&self, x: &[f64], y: &[f64], pre_cdf: &[f64]) -> Vec<f64> {
        let size = x.len();
        let bw = if self.rnaseq {
            0.5
        } else {
            x.stat(1).1 / self.sigma
        };
        let mut row = vec![0.0; size];
        for j in 0..y.len() {
            let mut left_tail = 0.0;
            for i in 0..x.len() {
                if self.rnaseq {
                    // ppois(y[j], x[i]+bw, TRUE, FALSE):
                    // this function returns the value of the Poisson cumulative density function
                    // NOTE: input has to be intergers
                    let pois = Poisson::new(x[i] + bw).unwrap();
                    left_tail += pois.cdf(y[j] as u64) as f64;
                } else {
                    left_tail += self.precomputed_cdfs(y[j] - x[i], bw, pre_cdf);
                }
            }
            left_tail = left_tail / (size as f64);
            row[j] = -1.0 * ((1.0 - left_tail) / left_tail).ln();
        }
        // println!("s1 {:?}, s2 {:?}, row: {:?}", x.len(), y.len(), &row);
        return row;
    }
    /// mat_density [n_genes, n_samples]
    /// return  gene_density matrix [n_genes, n_samples]
    fn mat_d(&self, mat_density: &[Vec<f64>], pre_cdf: &[f64]) -> Vec<Vec<f64>> {
        // let D = vec![vec![0;mat[0].lend()]; mat.len()];
        mat_density
            .par_iter()
            .enumerate()
            .map(|(i, vv)| {
                // print!("gene idx: {:?}", i);
                self.row_d(vv, vv, pre_cdf)
            })
            .collect()
    }

    /// mat: [n_genes, n_samples]
    /// return gene_density matrix [n_genes, n_samples]
    pub fn compute_density(&self, mat: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if self.kcdf {
            let pre_cdf = self.init_cdfs(self.pre_res, self.max_pre);
            let mat3 = self.mat_d(mat, &pre_cdf);
            return mat3;
        }
        return self.mat_ecdf(mat);
    }
    /// mat [n_samples, n_genes]
    pub fn compute_rank_score(&self, mat: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<usize>>) {
        let n = mat[0].len();
        let rev_idx: Vec<f64> = (1..=n)
            .into_iter()
            .rev()
            .map(|v| ((v as f64) - (n as f64) / 2.0).abs())
            .collect();
        // println!("{:?}", &rev_idx);
        let mut idxs: Vec<Vec<usize>> = vec![vec![0; n]; mat.len()];
        let mut mat2: Vec<Vec<f64>> = vec![vec![0.0; n]; mat.len()];
        // https://github.com/rcastelo/GSVA/blob/devel/R/gsva.R, line 865-874
        for i in 0..mat.len() {
            // ranking gene.density matrix in decending order
            idxs[i] = mat[i].as_slice().argsort(false).0;
            // get rank score, put it in original order
            let mut tmp: Vec<f64> = vec![0.0; n];
            // R code: tmp[sort_idx_vec] <- abs(seq(from=num_genes,to=1) - num_genes/2)
            idxs[i].iter().enumerate().for_each(|(ii, &j)| {
                tmp[j] = rev_idx[ii];
            });
            // println!("row index {:?}, {:?}", i, tmp);
            // get rank score
            mat2[i] = tmp;
        }
        (mat2, idxs)
    }
    /// compute rank score in parallel without matrix transposition
    fn compute_rank_score2(&self, mat: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<usize>>) {
        let n_samples = mat[0].len();
        let n_genes = mat.len();
    
        (0..n_samples).into_par_iter().map(|sample_idx| {
            let mut sample: Vec<f64> = Vec::with_capacity(n_genes);
            for gene in mat {
                sample.push(gene[sample_idx]);
            }
            
            let (sorted_idx, sorted_vals) = sample.as_slice().argsort(false);
            let rev_idx: Vec<f64> = (1..=n_genes)
                .rev()
                .map(|v| ((v as f64) - (n_genes as f64) / 2.0).abs())
                .collect();
    
            let mut tmp = vec![0.0; n_genes];
            sorted_idx.iter().enumerate().for_each(|(i, &j)| {
                tmp[j] = rev_idx[i];
            });
    
            (tmp, sorted_idx)
        }).unzip()
    }
    fn ks_sample(
        &self,
        gene_density: &[f64],
        sidx: &[usize],
        geneset_mask: &[usize],
        fset: &[usize],
    ) -> f64 {
        let tau = self.tau;
        let mut sum_gset = 0.0;
        let n_fset = fset.len();
        let n_genes = gene_density.len();
        let dec = 1.0 / ((n_genes - n_fset) as f64);

        for i in 0..n_fset {
            let tmp = gene_density[fset[i]];
            sum_gset += tmp.powf(tau);
        }

        let mut mx_value_sign = 0.0;
        let mut cum_sum = 0.0;
        let mut mx_pos = 0.0;
        let mut mx_neg = 0.0;
        for i in 0..n_genes {
            let idx = sidx[i]; // idx = sort_indxs[i];
            if geneset_mask[idx] == 1 {
                cum_sum += gene_density[idx].powf(tau) / sum_gset;
            } else {
                cum_sum -= dec;
            }
            if cum_sum > mx_pos {
                mx_pos = cum_sum;
            }
            if cum_sum < mx_neg {
                mx_neg = cum_sum;
            }
        }
        // mx value sign
        if self.mx_diff {
            mx_value_sign = mx_pos + mx_neg;
            if self.abs_rnk {
                mx_value_sign = mx_pos - mx_neg;
            }
        } else {
            mx_value_sign = if mx_pos > mx_neg.abs() {
                mx_pos
            } else {
                mx_neg
            }
        }
        return mx_value_sign;
    }
    /// fset : feature set
    /// d: [n_samples, n_genes], gene density scores
    /// idxs: sroted gene densities idxs
    pub fn ks_matrix(&self, d: &[Vec<f64>], sidxs: &[Vec<usize>], fset: &[usize]) -> Vec<f64> {
        // let n_sample = d.len();
        let n_genes = d[0].len();
        let mut geneset_mask = vec![0; n_genes]; // tag_indicator
        fset.iter().for_each(|i| {
            geneset_mask[*i] = 1;
        });
        // let dec = 1.0 / ((n_genes - n_geneset) as f64);
        d.par_iter()
            .zip(sidxs.par_iter())
            .map(|(vv, idx)| self.ks_sample(vv, idx, &geneset_mask, fset))
            .collect()
    }
}

fn transpose(mat: &[Vec<f64>]) -> Vec<Vec<f64>> {
    // transpose [m_gene, n_sample] --> [n_sample, m_gene]
    let mut mat2: Vec<Vec<f64>> = vec![vec![]; mat[0].len()];
    mat.iter().for_each(|row| {
        row.iter().enumerate().for_each(|(j, e)| {
            mat2[j].push(*e);
        });
    });
    return mat2;
}

pub fn gsva(
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
) -> GSEAResult {
    let es = GSVA::new(&gene_name, kcdf, tau, mx_diff, abs_rnk, rnaseq);
    let n_genes = gene_expr.len();
    let n_samples = gene_expr[0].len();

    // Precompute gene set hits in parallel
    let gene_set_hits: Vec<_> = gene_sets
        .into_par_iter()
        .filter_map(|(term, gset)| {
            let hits = es.genes.index_of_any(&gset);
            (!hits.is_empty() && hits.len() >= min_size && hits.len() <= max_size)
                .then(|| (term, hits))
        })
        .collect();

    // Process samples in parallel without transposing
    let (mat_score, sort_idxs) = es.compute_rank_score2(&es.compute_density(&gene_expr));

    // Process chunks of samples in parallel to control memory usage
    let chunk_size = (n_samples / rayon::current_num_threads()).max(1);

    // Parallel KS score calculation
    let summaries: Vec<_> = gene_set_hits
        .par_chunks(chunk_size)  // Process gene sets in chunks of 100
        .flat_map_iter(|chunk| {
            chunk.iter().flat_map(|(term, hits)| {
                let hits: Vec<usize> = hits.iter().copied().map(|&i| i).collect();
                es.ks_matrix(&mat_score, &sort_idxs, &hits)
                    .into_iter()
                    .enumerate()
                    .map(move |(i, es_val)| GSEASummary {
                        term: term.clone(),
                        es: es_val,
                        index: Some(i),
                        ..Default::default()
                    })
            })
        })
        .collect();

    let mut gs = GSEAResult::new(tau, max_size, min_size, 0, 0);
    gs.summaries = summaries;
    gs.indices = sort_idxs;
    gs.rankings = mat_score;
    return gs;
}

mod tests {
    use super::*;
    use crate::utils::FileReader;
    use csv::WriterBuilder;
    use std::fs::File;

    #[test]
    fn test_run_gsva() {
        let mut gene_name: Vec<String> = Vec::new();
        let mut gct = FileReader::new();
        let _ = gct.read_csv("tests/data/expr.gsva.csv", b',', true, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_table("tests/data/geneset.gsva.gmt", '\t', false);

        let mut gene_exp: Vec<Vec<f64>> = Vec::new();
        for r in gct.record.iter() {
            gene_name.push(r[0].to_string());
            let mut vv: Vec<f64> = Vec::new();
            for v in &r[1..] {
                vv.push(v.parse::<f64>().unwrap());
            }
            gene_exp.push(vv);
        }
        let mut gene_sets = HashMap::<String, Vec<String>>::new();
        gmt.record.iter().for_each(|r| {
            gene_sets.insert(r[0].to_string(), r[2..].to_vec());
        });
        let sample_names = &gct.header.get_vec()[1..];
        //let expr = transpose(&gene_exp);
        println!("{:?}", sample_names);
        let gsumm = gsva(
            gene_name, gene_exp, gene_sets, true, false, true, false, 1.0, 1, 1000,
        );
        let file = File::create("tests/data/gsva.rs.out").expect("Cannot write filename");
        let mut wrt = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);
        let _ = wrt.write_record(gct.header.get_vec().iter());
        let mut record = Vec::<String>::new();
        for (i, g) in gsumm.summaries.iter().enumerate() {
            if g.index.unwrap() == 0 {
                if !record.is_empty() {
                    let _ = wrt.write_record(&record);
                }
                record.clear();
                record.push(g.term.to_string());
            }
            record.push(g.es.to_string());
        }
        // don't forget write last record
        let _ = wrt.write_record(&record);
        wrt.flush();
    }
}
