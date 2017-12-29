#' @param X matrix. Rows are genes. Columns are samples. Row names are symbols.
#' @param gene_sets list. Each element is a string vector with gene symbols.
#' @param alpha numeric. Parameter for ssGSEA, the default is 0.25
#' @param scale logical. If True, normalize the scores by number of genes in the gene sets.
#' @param norm logical. If True, normalize the scores by the absolute difference between max and min values.
#' @param single logical. If True, use ssGSEA algorithm, otherwise use GSEA.
#'
#' @return matrix containing enrichment scroes. Rows are gene sets, columns are samples.
#'
#' @examples
#' # Create a fake matrix
# m = 100
# n = 100
# set.seed(1)
# X = matrix(rnorm(m*n), m, n)
# # # Assign 'gene symbols' to row names
# rownames(X) = 1:m
# # Create 3 gene sets
# gene_sets = list(a = sample(m, 5), b = sample(m, 5), c = sample(m, 5))
# system.time(assign('a', GSVA::gsva(X, gene_sets, method = 'ssgsea')))
# system.time(assign('b', ssgsea(X, gene_sets, scale = F, norm = T)))
# identical(a, b)
set.seed(1)
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)

    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos

      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

      step_cdf_diff = step_cdf_pos - step_cdf_neg

      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes

      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })

  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))

  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

# setwd("~/github/GSEApy/tests/data")
# comparison with gseapy.ssgsea
X2 = read.table("./data/testSet_rand1200.gct", row.names = 1, header = T,
                comment='#', sep="\t", stringsAsFactors = F)
X3 = as.matrix.data.frame(X2[,-1])
gene_sets = fgsea::gmtPathways("./data/randomSets.gmt")

# test
system.time(assign('a', GSVA::gsva(X3, gene_sets, method = 'ssgsea')))
system.time(assign('b', ssgsea(X3, gene_sets, scale = F, norm = T)))
system.time(assign('b2', ssgsea(X3, gene_sets, scale = T, norm = F)))

identical(a, b)
print("TEST with GSVA::gsva(method='ssgsea')")
print(a)

print("\n\n\n")
print("TEST with ssgsea, NES")
print(b)

print("\n\n\n")
print("TEST with ssgsea, Scaled ES")
print(b2)
