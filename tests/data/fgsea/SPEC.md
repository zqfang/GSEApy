# fgsea C++ ground-truth validation spec

Ground-truth outputs of the ORIGINAL fgsea C++ (Rcpp + BH) for validating the
Rust port for bit-identical output.

- Build: `Rcpp::sourceCpp("validation/harness.cpp")` driven by
  `validation/gen_refs.R`. The harness unity-includes a verbatim copy of the
  package sources under `validation/srccopy/` (RcppExports.cpp omitted; the
  `// [[Rcpp::export]]` / `// [[Rcpp::plugins]]` attribute comment lines
  neutralized; files renamed `*.cppinc` so sourceCpp does not compile them as
  separate translation units — that prevented duplicate-symbol link errors).
  The originals under `fgsea/src/` were never modified.
- Thin wrappers `ref_<fn>` forward to the originals with identical signatures.
- All floating-point outputs are written with `sprintf("%.17g", x)` (17
  significant digits, round-trippable doubles). Booleans as 0/1 ints.
- Index convention in emitted files: `*_index0` columns are 0-based to match
  Rust. Inputs that feed the C++ (selectedStats, geneRanks, ranks, pathway
  index lists) remain 1-based, exactly as the C++ API expects.

All reference inputs are dumped under `validation/inputs/` so the Rust side
reads them directly rather than re-deriving subsets.

---

## Common source data and shared subsets

### Ranks (`stats`)
- Source: `fgsea/inst/extdata/naive.vs.th1.rnk`, TSV with header `ID<TAB>t`.
- Sort all rows by `t` DECREASING (ties: R `order` stable order).
- SUBSET: take the FIRST 3000 rows of the sorted list.
- `validation/inputs/stats.tsv` — header `stats`, then 3000 doubles (the `t`
  values, in sorted order, index 0..2999).
- `validation/inputs/stats_ids.tsv` — the 3000 gene IDs in the same order (no
  header), for reference.

### Pathways
- Source: `fgsea/inst/extdata/mouse.reactome.gmt`, GMT `name<TAB>id<TAB>gene...`.
- For each pathway, map its gene IDs to 1-based indices into the 3000-gene
  `stats` ordering via `match(genes, stats_ids)`; drop NA and duplicates.
- Eligible = pathways whose mapped size is in [15, 300].
- SUBSET: order eligible by mapped size ascending; pick 20 at evenly spaced
  positions `round(seq(1, n_eligible, length.out=20))`, then dedup. Result is
  20 pathways spanning sizes 15..263.
- `validation/inputs/pathways.tsv` — 20 lines, each a TAB-separated list of
  SORTED 1-based indices for that pathway.
- `validation/inputs/pathway_names.tsv` — 20 pathway names (order matches).
- `validation/inputs/pathway_sizes.tsv` — 20 sizes (order matches).

---

## 1. `calcGseaStatCumulative(stats, selectedStats, gseaParam, scoreType)`

- `stats`: the 3000 doubles from `inputs/stats.tsv`.
- `selectedStats`: 1-based sorted indices of ONE representative pathway — the
  pick whose size is closest to the median pathway size. Dumped at
  `inputs/cumulative_selectedStats.tsv` (one index per line). Its position in
  the 20-pathway list is in `inputs/cumulative_pathway_index0.txt` (0-based).
- `gseaParam` = 1.0.
- `scoreType` runs over {"std","pos","neg"}.
- Output `refs/ref_calcGseaStatCumulative.tsv`, header
  `scoreType<TAB>prefix_index0<TAB>value`. One row per prefix (0-based prefix
  index, 0..k-1) per scoreType. `value` = ES of that prefix.

## 2. `calcGseaStatBatchCpp(stats, selectedGenes, geneRanks)`

- `stats`: `abs(stats)^1` (gseaParam=1) of the 3000 ranks. Dumped at
  `inputs/batchcpp_stats.tsv` (header `statsAbs`).
- `geneRanks`: identity `1..3000`. Dumped at `inputs/batchcpp_geneRanks.tsv`
  (one per line).
- `selectedGenes`: the 20 pathway index lists from `inputs/pathways.tsv`
  (1-based; under identity geneRanks these equal rank positions).
- Output `refs/ref_calcGseaStatBatchCpp.tsv`, header
  `pathway_index0<TAB>pathway_name<TAB>size<TAB>ES`. One row per pathway
  (0-based), ES = single enrichment score.

## 3. `calcGseaStatCumulativeBatch(stats, gseaParam, pathwayScores, pathwaysSizes, iterations, seed, scoreType)`

- `stats`: the 3000 doubles (`inputs/stats.tsv`).
- `gseaParam` = 1.0.
- `pathwayScores`: per-pathway ES = last-prefix value of
  `calcGseaStatCumulative(stats, pathway, 1.0, "std")`, for the 20 pathways.
  Dumped at `inputs/cumbatch_pathwayScores.tsv` (header `pathwayScores`).
- `pathwaysSizes`: the 20 sizes. Dumped at `inputs/cumbatch_pathwaysSizes.tsv`.
- `iterations` = 1000, `seed` = 42, `scoreType` = "std". (Also in
  `inputs/cumbatch_params.tsv`.)
- RNG: `boost::mt19937` seeded with `static_cast<uint32_t>(seed)`; sampling via
  `combination(1, n, k, rng)` (see `util.cppinc`). The Rust port MUST reproduce
  this Mersenne-Twister stream and the `uid_wrapper` rejection sampling exactly.
- Output `refs/ref_calcGseaStatCumulativeBatch.tsv`, header
  `pathway_index0<TAB>leEs<TAB>geEs<TAB>leZero<TAB>geZero<TAB>leZeroSum<TAB>geZeroSum`.
  One row per pathway (0-based); the six accumulator vectors from the returned
  list.

## 4. `fgseaMultilevelCpp(enrichmentScores, ranks, pathwaySize, sampleSize, seed, eps, sign, moveScale, logStatus)`

- `ranks`: INTEGER vector `1..3000` (must be INTSXP). Dumped at
  `inputs/multilevel_ranks.tsv` (one per line).
- `enrichmentScores`: fixed probe set `{0.3, 0.5, 0.7, -0.4, -0.6}`. Dumped at
  `inputs/multilevel_enrichmentScores.tsv` (header `ES`).
- `pathwaySize` runs over {25, 100} (`inputs/multilevel_pathwaySizes.tsv`).
- `sampleSize` = 101, `seed` = 42, `eps` = 0.0, `sign` = FALSE,
  `moveScale` = 1.0, `logStatus` = FALSE. (`inputs/multilevel_params.tsv`.)
- RNG: `boost::mt19937` seeded with `seed`; multilevel resampling in
  `fgseaMultilevelSupplement.cppinc`. Rust must match the MT stream and the
  `score_t` exact-integer comparison logic.
- Output `refs/ref_fgseaMultilevelCpp.tsv`, header
  `pathwaySize<TAB>es_index0<TAB>ES<TAB>cppMPval<TAB>cppIsCpGeHalf<TAB>log2err`.
  One row per (pathwaySize, ES) pair; `cppIsCpGeHalf` is 0/1.

## 5. `gesecaCpp(E, inpScores, genesetSize, sampleSize, seed, eps)`

- `E`: expression matrix from `fgsea/inst/extdata/gse14308_expression_matrix.tsv`
  (10000 genes x 12 samples; first column = gene id used as row name).
  SUBSET: first 200 rows (genes) x all 12 columns, then ROW-CENTER each row
  (subtract row mean; no scaling). Dumped at `inputs/geseca_E.tsv` — 200 rows,
  12 TAB-separated doubles per row, NO header (row-major, matches `E` passed to
  C++ exactly via `%.17g`).
- `inpScores`: fixed `{5.0, 10.0, 20.0}` (`inputs/geseca_inpScores.tsv`, header
  `score`). These are the geneset scores at which p-values are looked up; the
  ruler is extended to `max(inpScores)=20`.
- `genesetSize` = 20, `sampleSize` = 101, `seed` = 42, `eps` = 0.0.
  (`inputs/geseca_params.tsv`.)
- Score function: `getScore(profile) = sum(profile[j]^2)` over the 12 columns,
  where `profile[j] = sum over geneset rows of E[row][j]` (see
  `ScoreCalculation.cppinc`). RNG = `boost::mt19937(seed)`.
- Output `refs/ref_gesecaCpp.tsv`, header
  `score_index0<TAB>score<TAB>pval<TAB>log2err`. One row per inpScore (0-based).

---

## Files

```
validation/
  harness.cpp          unity-include harness + ref_ wrappers
  gen_refs.R           build + generate everything
  srccopy/*.cppinc,*.h verbatim copy of fgsea/src (attributes neutralized)
  SPEC.md              this file
  inputs/              exact input arrays fed to each entry point
  refs/                full-precision reference outputs (5 TSVs)
```

## Reproduce

```
cd /home/fangzq/github/fgsea-rs
Rscript validation/gen_refs.R
```

Requires R 4.5 with Rcpp + BH, Boost headers at /usr/include/boost, g++.
Deterministic: all seeds fixed (42). Outputs overwrite under `refs/`/`inputs/`.
