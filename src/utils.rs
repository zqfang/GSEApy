#![allow(dead_code, unused)]

use pyo3::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
// use std::path::{Path};
use csv;
use itertools::Itertools;
use rand::seq::SliceRandom;
use rand::Rng;
use std::collections::HashMap;
use std::hash::Hash;

#[pyclass]
#[derive(Copy, Clone)]
pub enum Metric {
    Signal2Noise = 0,
    AbsSignal2Noise = 1,
    Ttest = 2,
    RatioOfClasses = 3,
    Log2RatioOfClasses = 4,
    DiffOfClasses = 5,
}

#[pyclass]
#[derive(Copy, Clone)]
pub enum CorrelType {
    Rank = 0,
    SymRank = 1,
    ZScore = 2,
}

pub trait Statistic {
    fn mean(&self) -> f64;
    fn stat(&self, ddof: usize) -> (f64, f64);
    fn argsort(&self, assending: bool) -> (Vec<usize>, Vec<f64>);
}

impl Statistic for &[f64] {
    /// caculate mean
    fn mean(&self) -> f64 {
        let sum = self.iter().sum::<f64>();
        let count = self.len() as f64;
        sum / count
    }
    /// return (mean, std), don't know why this is very slow
    fn stat(&self, ddof: usize) -> (f64, f64) {
        let sum: f64 = self.iter().sum();
        let count = self.len();
        let mean = sum / (count as f64);
        let variance = self
            .iter()
            .map(|&value| {
                let diff = mean - value;
                diff * diff
            })
            .sum::<f64>()
            / ((count - ddof) as f64);
        (mean, variance.sqrt())
    }
    fn argsort(&self, ascending: bool) -> (Vec<usize>, Vec<f64>) {
        let indices: Vec<usize> = (0..self.len()).collect();
        let sorted_col: Vec<(usize, &f64)> = indices
            .into_iter()
            .zip(self.iter())
            .sorted_by(|&a, &b| a.1.partial_cmp(b.1).unwrap())
            .collect();
        //.sorted_by(|&a, &b| a.1.partial_cmp(b.1).unwrap()).collect();

        let mut sidx: Vec<usize> = Vec::new();
        let mut sval: Vec<f64> = Vec::new();
        sorted_col.iter().for_each(|(i, &v)| {
            sidx.push(*i);
            sval.push(v);
        });
        if !ascending {
            sidx.reverse(); // inplace
            sval.reverse();
        }
        (sidx, sval)
    }
}

/// Dynamic Enum
#[derive(Debug, Clone)]
pub struct DynamicEnum<T> {
    _elt_to_idx: HashMap<T, usize>, // element to index
    _idx_to_elt: Vec<T>,            // index to element
    _num_indices: usize,            // size
}

impl<T> DynamicEnum<T>
where
    T: Eq + Hash + Clone,
{
    /// an empty object
    pub fn new() -> Self {
        DynamicEnum {
            _num_indices: 0,
            _idx_to_elt: Vec::<T>::new(),
            _elt_to_idx: HashMap::<T, usize>::new(),
        }
    }
    /// construct from vec
    pub fn from(vec: &[T]) -> Self {
        //let temp = vec.to_vec();
        let v2m: HashMap<T, usize> = vec
            .iter()
            .enumerate()
            .map(|(i, v)| (v.clone(), i))
            .collect();
        DynamicEnum {
            _num_indices: v2m.len(),
            _elt_to_idx: v2m,
            _idx_to_elt: vec.to_vec(),
        }
    }
    /// add element if new
    /// return indices whether new or not
    pub fn add_if_new(&mut self, element: T) -> usize {
        if self._elt_to_idx.contains_key(&element) {
            return *self._elt_to_idx.get(&element).unwrap();
        } 
        let key = element.clone();
        let idx = self._num_indices;
        self._idx_to_elt.push(element);
        self._elt_to_idx.insert(key, idx);
        self._num_indices += 1;
        return idx;
    }
    /// get index of element
    pub fn index_of(&self, element: &T) -> Option<&usize> {
        self._elt_to_idx.get(element)
    }
    pub fn index_of_any(&self, elements: &[T]) -> Vec<&usize> {
        elements.iter()
        .filter_map(|e| self.index_of(e))
        .collect()
    }
    #[allow(dead_code)]
    pub fn contain_elt(&self, element: &T) -> bool {
        self._elt_to_idx.contains_key(element)
    }
    /// get element at position of index
    pub fn elt_of(&self, idx: usize) -> Option<&T> {
        self._idx_to_elt.get(idx)
    }
    /// return indicator whether the self.elements in given elements (0: absent, 1: present)
    pub fn isin(&self, elements: &[T]) -> Vec<f64> {
        let mut _tag_indicator: Vec<f64> = vec![0.0; self._idx_to_elt.len()];
        elements.iter().for_each(|e| {
            if let Some(idx) = self.index_of(e) {
                _tag_indicator[*idx] = 1.0;
            }
        });
        return _tag_indicator;
    }
    pub fn size(&self) -> usize {
        return self._num_indices;
    }
    pub fn get_vec(&self) -> &Vec<T> {
        return &self._idx_to_elt;
    }

    /// inplace shuffle
    pub fn shuffle<R>(&mut self, rng: &mut R)
    where
        R: Rng + ?Sized,
    {
        self._idx_to_elt.shuffle(rng);
        self._idx_to_elt.iter().enumerate().for_each(|(i, e)| {
            self._elt_to_idx.insert(e.clone(), i);
        });
    }
}

pub struct FileReader {
    lineno: usize,
    pub header: DynamicEnum<String>,
    pub record: Vec<Vec<String>>,
}
impl FileReader {
    pub fn new() -> Self {
        FileReader {
            lineno: 0,
            header: DynamicEnum::<String>::new(),
            record: Vec::<Vec<String>>::new(), // 2d vec init
        }
    }
    /// delimiter, comment: byte char literal input, e.g.  b'\t', b'#'
    pub fn read_csv(
        &mut self,
        file_path: &str,
        delimiter: u8,
        header: bool,
        comment: Option<u8>,
    ) -> Result<(), Box<dyn Error>> {
        // u8: b'\t', b'#' ...
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .has_headers(header)
            .comment(comment)
            .from_path(file_path)
            .expect("Cannot read {file_path}");
        if header {
            let _h = rdr.headers()?; // Result<StringRecord, Error>, check
            for col in _h.iter() {
                self.header.add_if_new(col.to_string());
            }
        }
        for result in rdr.records() {
            // The iterator yields Result<StringRecord, Error>, so we check the error here.
            self.lineno += 1;
            let record = result?;
            //println!("{:?}", record);
            self.record
                .push(record.iter().map(|x| x.to_string()).collect());
        }
        Ok(())
    }
    pub fn to_csv(
        &mut self,
        file_path: &str,
        delimiter: u8,
        header: bool,
    ) -> Result<(), Box<dyn Error>> {
        let file = File::create(file_path).expect("Cannot write {file_path}");
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_writer(file);
        if header {
            wtr.write_record(self.header.get_vec())?;
        }
        for rec in self.record.iter() {
            wtr.write_record(rec)?;
        }
        wtr.flush()?;
        Ok(())
    }

    pub fn read_table(
        &mut self,
        file_path: &str,
        delimiter: char,
        header: bool,
        //comment: Option<u8>,
    ) -> Result<(), Box<dyn Error>> {
        let input = File::open(file_path)?;
        let mut buffered = BufReader::new(input);
        if header {
            let mut _header: String = String::new();
            let _ = buffered
                .read_line(&mut _header)
                .expect("Read header failed");
            for col in _header.split(delimiter) {
                self.header.add_if_new(col.to_string());
            }
        }
        for (_num, line) in buffered.lines().enumerate() {
            let record: Vec<String> = line?.split(delimiter).map(|x| x.to_string()).collect();
            self.record.push(record);
            self.lineno = _num;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_rdr() {
        let cwd = std::env::current_dir().unwrap(); // prjoject root, directory to Cargo.toml
        let rnk_path = cwd.join("tests/data/edb/gsea_data.gsea_data.rnk");
        let gmt_path = cwd.join("tests/data/edb/gene_sets.gmt");
        println!("{:?}", &rnk_path);
        let mut rnk = FileReader::new();
        let _ = rnk.read_csv(rnk_path.to_str().unwrap(), b'\t', false, Some(b'#'));
        let mut gmt = FileReader::new();
        let _ = gmt.read_csv(gmt_path.to_str().unwrap(), b'\t', false, None);
        let _wr = rnk.to_csv("example.output.txt", b'\t', true);
    }
    #[test]
    fn test_dynum() {
        //let dynum = DynamicEnum::new();
        let vec = vec!["A", "B", "C", "D"];
        let mut dynum = DynamicEnum::from(&vec);
        let x = dynum.add_if_new("E");
        let y = *dynum.index_of(&"C").unwrap();
        let z = *dynum.elt_of(1).unwrap();

        println!("{x}, {y}, {z}");
        assert_eq!(x, 4);
        assert_eq!(y, 2);
        assert_eq!(z, "B");
    }
}
