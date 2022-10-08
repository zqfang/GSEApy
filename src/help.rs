/// phenotype permutation procedure
/// shuffling group labels and calculate the new ranking metric
/// return shuffled metric (not sorted)
pub fn phenotype_permutation(data: &[Vec<f64>], group: &[bool], method: Metric) -> Vec<Vec<f64>> {
    //let mut indices: Vec<Vec<usize>> = Vec::new();
    let mut arr: Vec<Vec<f64>> = Vec::new();
    let mut group_nperm = vec![group.to_vec(); self.nperm + 1];
    // let mut g = group.to_vec();
    for i in 1..=self.nperm {
        //fastrand::shuffle(&mut group_arc);
        group_nperm[i].shuffle(&mut self.rng);
    }
    // let metric = Arc::new(data.to_vec());
    let metric = data.to_vec();
    // thread pool
    let thread: usize = 8;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread)
        .build()
        .unwrap();
    let _method = Arc::new(method);
    // creat channel
    let (sender, receiver) = mpsc::channel();
    // The scope is guaranteed to exit only after all tasks launched inside it finished.
    // This essentially allows the tasks inside the scope to access variables that live at least as long as the scope:

    // TODO: need to keep order of the results
    pool.scope(move |s| {
        for i in 0..=self.nperm {
            let _m = Arc::clone(&_method);
            let _metric = Arc::clone(&metric);
            let group_rng = Arc::new(group_nperm[i].clone());
            let tx = sender.clone();

            s.spawn(move |s| {
                let arr: Vec<f64> = _metric
                    .iter()
                    .map(|vec| {
                        // select pos or neg data by a mask array
                        //let end1 = Instant::now();
                        // TODO: the interator is very slow here
                        let mut pos: Vec<f64> = Vec::new();
                        let mut neg: Vec<f64> = Vec::new();
                        vec.iter().zip(group_rng.iter()).for_each(|(&x, &b)| {
                            if b {
                                pos.push(x);
                            } else {
                                neg.push(x);
                            }
                        });
                        //let end2 = Instant::now();
                        let pos_len = pos.len() as f64;
                        let neg_len = neg.len() as f64;
                        let (pos_mean, pos_std) = pos.stat(1);
                        let (neg_mean, neg_std) = neg.stat(1);
                        // let end3 = Instant::now();

                        let m = *_m;
                        let r = match m {
                            Metric::Signal2Noise => (pos_mean - neg_mean) / (pos_std + neg_std),
                            Metric::AbsSignal2Noise => {
                                ((pos_mean - neg_mean) / (pos_std + neg_std)).abs()
                            }
                            Metric::Ttest => {
                                (pos_mean - neg_mean)
                                    / (pos_std * pos_std / pos_len + neg_std * neg_std / neg_len)
                                        .sqrt()
                            }
                            Metric::RatioOfClasses => pos_mean / neg_mean,
                            Metric::Log2RatioOfClasses => (pos_mean / neg_mean).log2(),
                            Metric::DiffOfClasses => pos_mean - neg_mean,
                        };

                        // println!("calulation time 1: {:.2?}", end2.duration_since(end1));
                        // println!("calulation time 2: {:.2?}", end3.duration_since(end2));
                        return r;
                    })
                    .collect();
                tx.send(arr).unwrap();
            });
        }
    });
    // drop the original sender, else the channel will remain open, causing the receiver to infinitely wait
    // mem::drop(sender);
    // The receiver will block when the iterator asks for the next value.
    // When the channel is closed, the iterator will return None and end.
    for received in receiver {
        //indices.push(received.0);
        arr.push(received);
    }
    return arr;
}
