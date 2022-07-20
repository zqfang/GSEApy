use rayon::prelude::*;
// use std::io;
use clap::Parser;
use gsea::*;
use std::time::Instant;
/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(version = "0.1.0", about = "Gene Enrichment Analysis in Rust", long_about = None)]
struct Args {
    /// Input file name
    #[clap(short, long)]
    input: String,
    #[clap(short, long)]
    gmt: String,
    #[clap(short, long)]
    cls: String,
    /// Output file name
    #[clap(short, long)]
    output: String,

    #[clap(short, long, default_value_t = 1.0)]
    weight_score: f64,
    #[clap(short, long, default_value_t = 1000)]
    permutation_num: usize,
    #[clap(long, default_value_t = 15)]
    min_size: u64,
    #[clap(long, default_value_t = 500)]
    max_size: u64,

    #[clap(short, long, default_value_t = 666)]
    seed: u64,
    #[clap(short, long, default_value_t = 1)]
    threads: usize,
}

fn main() {
    let args = Args::parse();
    // set thread number global
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    // read data
    let mut gct = FileReader::new();
    let _ = gct.read_csv(&args.input, b'\t', true, Some(b'#'));
    let mut gmt = FileReader::new();
    let _ = gmt.read_table(&args.gmt, '\t', false);
    let mut cls = FileReader::new();
    let _ = cls.read_table(&args.cls, ' ', false);
    // let mut group = cls.record.pop().unwrap();
    let group = &cls.record[1];
    let gboo: Vec<bool> = cls.record[2].iter().map(|x| *x == group[1]).collect();

    let weight = args.weight_score;
    let mut gene: Vec<String> = Vec::new();
    // let mut gene_set: Vec<String> = Vec::new();
    let mut gene_exp: Vec<Vec<f64>> = Vec::new();
    for r in gct.record.iter() {
        gene.push(r[0].clone());
        let mut vv: Vec<f64> = Vec::new();
        for v in &r[2..] {
            vv.push(v.parse::<f64>().unwrap());
        }
        vv.iter_mut().for_each(|x| *x = x.abs().powf(weight)); // modified value in place
        gene_exp.push(vv);
    }
    let mut es = EnrichmentScore::new(&gene, 100, 0, false, false);
    let start = Instant::now();
    let gene_metric = es.phenotype_permutation(&gene_exp, &gboo, Metric::Signal2Noise);
    let end1 = Instant::now();
    println!("permutation time: {:.2?}", end1.duration_since(start));
    let mut summ = Vec::<GSEASummary>::new();

}

mod tests {
    use super::*;
    use std::time::Instant;
    #[test]
    fn pipeline() {
        // set number of threads of rayon at the main()
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build_global()
            .unwrap();

    }
}
