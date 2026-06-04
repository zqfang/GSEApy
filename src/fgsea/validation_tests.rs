//! End-to-end bit-exactness validation of the vendored fgsea core against
//! ground-truth outputs of the ORIGINAL fgsea C++ (generated via `Rcpp::sourceCpp`
//! with `boost::mt19937`). Inputs and reference outputs live under
//! `tests/data/fgsea/` (copied verbatim from the fgsea-rs `validation/` tree).
//!
//! Tolerance rationale: with `eps = 0`, `digamma`/`trigamma` (the only functions
//! that are not bit-exact vs Boost) affect ONLY the final p-value arithmetic, not
//! the RNG sampling trajectory. So integer-count outputs and the classic GSEA
//! statistics match EXACTLY; only `pval`/`log2err` may drift by a few ULPs.

use std::fs;
use std::path::PathBuf;

use crate::fgsea::{
    calcGseaStatBatchCpp, calcGseaStatCumulative, calcGseaStatCumulativeBatch, fgseaMultilevelCpp,
};

fn base(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/fgsea").join(rel)
}

fn read_to_string(rel: &str) -> String {
    fs::read_to_string(base(rel)).unwrap_or_else(|e| panic!("read {rel}: {e}"))
}

/// One value per line; optionally skip a one-line header.
fn load_f64s(rel: &str, has_header: bool) -> Vec<f64> {
    let s = read_to_string(rel);
    s.lines()
        .skip(if has_header { 1 } else { 0 })
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse::<f64>().unwrap())
        .collect()
}

fn load_i32s(rel: &str, has_header: bool) -> Vec<i32> {
    let s = read_to_string(rel);
    s.lines()
        .skip(if has_header { 1 } else { 0 })
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse::<i32>().unwrap())
        .collect()
}

/// Each line is a TAB-separated list of i32 (one pathway per line).
fn load_rows_i32(rel: &str) -> Vec<Vec<i32>> {
    let s = read_to_string(rel);
    s.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').map(|t| t.trim().parse::<i32>().unwrap()).collect())
        .collect()
}

/// Parse a ref TSV (with header) into (header, rows-of-fields).
fn load_tsv(rel: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let s = read_to_string(rel);
    let mut lines = s.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<String> = lines.next().unwrap().split('\t').map(|x| x.to_string()).collect();
    let rows: Vec<Vec<String>> = lines
        .map(|l| l.split('\t').map(|x| x.to_string()).collect())
        .collect();
    (header, rows)
}

/// |a-b| <= atol + rtol*|b|. Returns (passed, abs_diff, rel_diff).
fn close(a: f64, b: f64, atol: f64, rtol: f64) -> (bool, f64, f64) {
    let abs = (a - b).abs();
    let rel = abs / b.abs().max(f64::MIN_POSITIVE);
    (abs <= atol + rtol * b.abs(), abs, rel)
}

/// Assert two vectors agree within tolerance; report the worst element.
fn assert_vec_close(label: &str, got: &[f64], want: &[f64], atol: f64, rtol: f64) {
    assert_eq!(got.len(), want.len(), "{label}: length {} != {}", got.len(), want.len());
    let mut worst = (0usize, 0.0f64, 0.0f64);
    for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
        let (ok, abs, rel) = close(g, w, atol, rtol);
        if !ok {
            panic!("{label}[{i}]: got {g:.17e}, want {w:.17e} (abs={abs:.3e}, rel={rel:.3e})");
        }
        if abs > worst.1 {
            worst = (i, abs, rel);
        }
    }
    eprintln!(
        "{label}: OK ({} values), worst abs={:.3e} rel={:.3e} at [{}]",
        got.len(),
        worst.1,
        worst.2,
        worst.0
    );
}

fn assert_exact(label: &str, got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len(), "{label}: length mismatch");
    for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
        assert_eq!(g, w, "{label}[{i}]: got {g}, want {w} (expected EXACT equality)");
    }
    eprintln!("{label}: OK ({} values, exact)", got.len());
}

// ---------------------------------------------------------------------------
// 1. calcGseaStatCumulative — pure int/double, expect (near-)exact.
// ---------------------------------------------------------------------------
#[test]
fn validate_calcGseaStatCumulative() {
    let stats = load_f64s("inputs/stats.tsv", true);
    let selected = load_i32s("inputs/cumulative_selectedStats.tsv", false);

    let (header, rows) = load_tsv("refs/ref_calcGseaStatCumulative.tsv");
    assert_eq!(header, ["scoreType", "prefix_index0", "value"]);

    for score_type in ["std", "pos", "neg"] {
        let want: Vec<f64> = rows
            .iter()
            .filter(|r| r[0] == score_type)
            .map(|r| r[2].parse::<f64>().unwrap())
            .collect();
        let got = calcGseaStatCumulative(&stats, &selected, 1.0, score_type);
        assert_vec_close(&format!("cumulative[{score_type}]"), &got, &want, 1e-12, 1e-12);
    }
}

// ---------------------------------------------------------------------------
// 2. calcGseaStatBatchCpp — pure int/double, expect (near-)exact.
// ---------------------------------------------------------------------------
#[test]
fn validate_calcGseaStatBatchCpp() {
    let stats = load_f64s("inputs/batchcpp_stats.tsv", true);
    let gene_ranks = load_i32s("inputs/batchcpp_geneRanks.tsv", false);
    let selected_genes = load_rows_i32("inputs/pathways.tsv");

    let (header, rows) = load_tsv("refs/ref_calcGseaStatBatchCpp.tsv");
    assert_eq!(header[3], "ES");
    let want: Vec<f64> = rows.iter().map(|r| r[3].parse::<f64>().unwrap()).collect();

    let got = calcGseaStatBatchCpp(&stats, &selected_genes, &gene_ranks);
    assert_vec_close("batchCpp.ES", &got, &want, 1e-12, 1e-12);
}

// ---------------------------------------------------------------------------
// 3. calcGseaStatCumulativeBatch — RNG-driven; counts must be EXACT.
// ---------------------------------------------------------------------------
#[test]
fn validate_calcGseaStatCumulativeBatch() {
    let stats = load_f64s("inputs/stats.tsv", true);
    let pathway_scores = load_f64s("inputs/cumbatch_pathwayScores.tsv", true);
    let pathways_sizes = load_i32s("inputs/cumbatch_pathwaysSizes.tsv", false);

    // params per SPEC.md: iterations=1000, seed=42, gseaParam=1, scoreType="std"
    let got = calcGseaStatCumulativeBatch(&stats, 1.0, &pathway_scores, &pathways_sizes, 1000, 42, "std");

    let (header, rows) = load_tsv("refs/ref_calcGseaStatCumulativeBatch.tsv");
    assert_eq!(header, ["pathway_index0", "leEs", "geEs", "leZero", "geZero", "leZeroSum", "geZeroSum"]);
    let col = |j: usize| -> Vec<f64> { rows.iter().map(|r| r[j].parse::<f64>().unwrap()).collect() };

    // Integer counts: exact equality (any drift => the MT19937/combination stream diverged).
    assert_exact("cumBatch.leEs", &got.leEs, &col(1));
    assert_exact("cumBatch.geEs", &got.geEs, &col(2));
    assert_exact("cumBatch.leZero", &got.leZero, &col(3));
    assert_exact("cumBatch.geZero", &got.geZero, &col(4));
    // Sums of doubles over a matched RNG trajectory: expect ~exact.
    assert_vec_close("cumBatch.leZeroSum", &got.leZeroSum, &col(5), 1e-9, 1e-12);
    assert_vec_close("cumBatch.geZeroSum", &got.geZeroSum, &col(6), 1e-9, 1e-12);
}

// ---------------------------------------------------------------------------
// 4. fgseaMultilevelCpp — RNG trajectory exact (eps=0); pval/log2err ~ULP drift.
// ---------------------------------------------------------------------------
#[test]
fn validate_fgseaMultilevelCpp() {
    let ranks = load_i32s("inputs/multilevel_ranks.tsv", false);
    let es = load_f64s("inputs/multilevel_enrichmentScores.tsv", true);

    let (header, rows) = load_tsv("refs/ref_fgseaMultilevelCpp.tsv");
    assert_eq!(header, ["pathwaySize", "es_index0", "ES", "cppMPval", "cppIsCpGeHalf", "log2err"]);

    // params per SPEC.md: sampleSize=101, seed=42, eps=0, sign=false, moveScale=1
    for pathway_size in [25, 100] {
        let res = fgseaMultilevelCpp(&es, &ranks, pathway_size, 101, 42, 0.0, false, 1.0, false);

        let want_rows: Vec<&Vec<String>> =
            rows.iter().filter(|r| r[0].parse::<i32>().unwrap() == pathway_size).collect();
        assert_eq!(want_rows.len(), es.len());

        let want_pval: Vec<f64> = want_rows.iter().map(|r| r[3].parse::<f64>().unwrap()).collect();
        let want_iscp: Vec<f64> = want_rows.iter().map(|r| r[4].parse::<f64>().unwrap()).collect();
        let want_log2: Vec<f64> = want_rows.iter().map(|r| r[5].parse::<f64>().unwrap()).collect();

        let iscp: Vec<f64> = res.cppIsCpGeHalf.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect();
        assert_exact(&format!("multilevel[{pathway_size}].cppIsCpGeHalf"), &iscp, &want_iscp);
        assert_vec_close(&format!("multilevel[{pathway_size}].pval"), &res.cppMPval, &want_pval, 1e-12, 1e-6);
        assert_vec_close(&format!("multilevel[{pathway_size}].log2err"), &res.log2err, &want_log2, 1e-9, 1e-6);
    }
}
