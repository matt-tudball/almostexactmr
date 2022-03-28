#' @name run_test
#' @title Compute almost exact p-values and confidence intervals
#' @description For a set of null hypotheses, computes the observed test statistic,
#' null distribution of test statistics, p-values and confidence intervals
#'
#' @param reps Number of resamples
#' @param beta Vector or scalar of null hypotheses (default is 0)
#' @param pheno List of phenoa which must be of the form \code{list(exp=, out=, cov=)}.
#' \code{exp} is a vector containing the exposure variable, \code{out} is a vector
#' containing the outcome variable and \code{cov} is a matrix containing the covariates
#' to be included in the test statistic.
#' @param prob Sampling probabilities from the \code{prop_score} function.
#' @param ins Matrix of instruments.
#' @param sig Significance level (default is 0.05).
#' @param cores Number of cores to use during parallel computing (default is 1).
#' @param out Vector of strings indicating which objects to return from the function.
#' Options include: "pvalues" for p-values corresponding to \code{beta}; "ci" for a \code{sig}-level
#' confidence interval obtained by inverting the test; "tobs" for the observed test statistic from
#' \code{test_stat}; and "tnull" for the full vector or matrix of counterfactual
#' test statistics (default is to return \code{c("pvalues","ci","tobs")}).
#'
#' @return Named list of objects requested in \code{out}.
#'
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#'
#' @examples
#' @export
run_test <- function(reps, beta=0, CHap, PHap, pheno, prob, snps, sig=0.05, cores=1, out=c("pvalues","tobs")) {
  n <- length(prob)

  ins <- CHap[,colnames(CHap) %in% snps]

  tobs <- test_stat(beta=beta, pheno=pheno, ins=ins)

  cluster <- make_clusters(cores)

  null_stat <- function(j) {
    new.ins <- sample_haplotype(prob)
    return(test_stat(beta, pheno, new.ins))
  }

  tnull <- matrix(t(pbapply::pbsapply(X=1:reps, cl = cluster, simplify=T, FUN=null_stat)), ncol=length(beta))

  on.exit(parallel::stopCluster(cluster))

  pvalues <- sapply(1:length(beta), function(x) 1-mean(tobs[x] >= tnull[,x]))

  main <- sapply(X=out, USE.NAMES=T, simplify=F, FUN=function(x) assign(x, eval(parse(text=x))))
  return(main)
}



