# GSEA 1.2 -- Gene Set Enrichment Analysis / Broad Institute Executable R
# script to run GSEA Analysis
library("utils")
library("tools")
library("dplyr")
library("GSEA")
cat("\n")
cat("Starting...\n")
cat("\n")
inputds <- readline(prompt = ("Input path to GCT or RNK formatted gene expression dataset (or drop file into R window): "))
if (file_ext(inputds) != "rnk") {
 inputcls <- readline(prompt = ("Input path to experiment CLS file (or drop file into R window): "))
}
gsdb <- readline(prompt = ("Input path to GMT formatted gene set database (or drop file into R window): "))

cat("\n")
if (file_ext(inputds) != "rnk") {
 collapsedataset <- askYesNo("Collapse data set to Gene Symbols? ")
 cat("\n")
} else {
 collapsedataset <- FALSE
}

if (collapsedataset == TRUE & !is.na(collapsedataset)) {
 inputchip <- readline(prompt = ("Input path to CHIP platform file (or drop file into R window): "))
 collapsetype <- menu(c("Max_probe", "Median_of_probes", "Mean_of_probes", "Sum_of_Probes"), 
  graphics = FALSE, title = "Collapsing mode for probe sets => 1 gene")
 if (collapsetype == 1) {
  collapsemode <- "max"
 } else if (collapsetype == 2) {
  collapsemode <- "median"
 } else if (collapsetype == 3) {
  collapsemode <- "mean"
 } else if (collapsetype == 4) {
  collapsemode <- "sum"
 }
} else {
 inputchip <- "NOCHIP"
 collapsemode <- "NOCOLLAPSE"
}

if (file_ext(inputds) != "rnk") {
 reshuffetype <- menu(c("Phenotype", "gene_set"), graphics = FALSE, title = "Select GSEA Permutation Type (recommended: Phenotype)")
 cat("\n")
} else {
 reshuffetype <- 2
}

npermsdefault <- askYesNo("Use default number of permutations for significance testing? (default: 1000) ")
if (npermsdefault == FALSE & !is.na(npermsdefault)) {
 nperms <- readline(prompt = ("Number of permutations: "))
} else {
 nperms <- 1000
}
cat("\n")

maxdefault <- askYesNo("Use default maximum gene set size filter? (default: 500 genes) ")
if (maxdefault == FALSE & !is.na(maxdefault)) {
 maxsize <- readline(prompt = ("Max size: "))
} else {
 maxsize <- 500
}
cat("\n")

mindefault <- askYesNo("Use default minimum gene set size filter? (default: 15 genes) ")
if (mindefault == FALSE & !is.na(mindefault)) {
 minsize <- readline(prompt = ("Min size: "))
} else {
 minsize <- 15
}

if (file_ext(inputds) != "rnk") {
 cat("\n")
 use.s2n <- askYesNo("Use default Signal2Noise metric for ranking genes? ")
 if (use.s2n == TRUE & !is.na(use.s2n)) {
  rankmetric <- "S2N"
 } else if (use.s2n == FALSE) {
  use.ttest <- askYesNo("Use T-Test to rank genes? ")
  if (use.ttest == TRUE & !is.na(use.ttest)) {
   rankmetric <- "ttest"
  } else {
   cat("No other ranking metrics implemented. Defaulting to S2N.\n")
   rankmetric <- "S2N"
  }
 }
}

cat("\n")

outdir <- readline(prompt = ("Drop a directory into R window to use as the output folder or enter directory path: "))
cat("\n")
outname <- readline(prompt = ("Enter a prefix to label output files: "))
cat("\n")

if (reshuffetype == 1) {
 permutation <- "sample.labels"
} else if (reshuffetype == 2) {
 permutation <- "gene.labels"
}

if (file_ext(inputds) != "rnk") {
 rankmethod <- "GSEA"
} else {
 rankmethod <- "preranked"
}

GSEA(
# Input/Output Files :-------------------------------------------------------------------------------
 input.ds = inputds,                    # Input gene expression dataset file in GCT format
 input.cls = inputcls,                  # Input class vector (phenotype) file in CLS format
 gs.db = gsdb,                          # Gene set database in GMT format
 input.chip = inputchip,               # CHIP File
 output.directory      = outdir,        # Directory where to store output and results (default: "")
#  Program parameters :-------------------------------------------------------------------------------
 doc.string            = outname,         # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
 reshuffling.type      = permutation,     # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
 nperm                 = as.integer(nperms),            # Number of random permutations (default: 1000)
 weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
 nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
 fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
 fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
 topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
 adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 gs.size.threshold.min = as.integer(minsize),         # Minimum size (in genes) for database gene sets to be considered (default: 25)
 gs.size.threshold.max = as.integer(maxsize),         # Maximum size (in genes) for database gene sets to be considered (default: 500)
 reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
 preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
 random.seed           = as.integer(as.POSIXct(Sys.time())),            # Random number generator seed. (default: 123456)
 perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
 fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
 replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
 collapse.dataset      = collapsedataset, # Collapse dataset to gene symbols using a user provided chip file (default: F)
 collapse.mode         = collapsemode,
 save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
 use.fast.enrichment.routine = T,         # Use faster routine to compute enrichment for random permutations (default: T)
 gsea.type = rankmethod,                     # Select Standard GSEA (default) or preranked
 rank.metric = rankmetric
)
#----------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
   directory           = outdir,        # Directory where to store output and results (default: "")
   topgs = 20,                                                           # number of top scoring gene sets used for analysis
   height = 16,
   width = 16,
   gsea.type = rankmethod,
   doc.string = outname
)
