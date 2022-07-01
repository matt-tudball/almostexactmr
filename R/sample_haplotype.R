#' @name sample_haplotype
#' @title Haplotype sampler
#' @description Samples a counterfactual offspring haplotype
#'
#' @param prob Sampling probabilities from \code{prop_score} function
#'
#' @return Vector of offspring haplotypes
#'
sample_haplotype <- function(prob) {
  # Number of combinations of SNPs
  num.snps <- ncol(prob)

  # Each combination of 0 and 1 to describe possible alleles
  vals_alleles <- sample_space(c(0,1), num.snps/2)

  which.row <- sapply(X=1:nrow(prob), FUN=function(x) {
    sample(1:num.snps, size=1, prob=prob[x,])
  })

  out <- vals_alleles[which.row,]

  return(out)
}

