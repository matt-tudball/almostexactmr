#' @name run_test
#' @title Compute almost exact p-values and confidence intervals
#' @description For a set of null hypotheses, computes the observed test statistic,
#' null distribution of test statistics, p-values and confidence intervals
#'
#' @param reps Number of resamples
#' @param beta Vector or scalar of null hypotheses (default is 0)
#' @param dat List of data which must be of the form \code{list(exp=, out=, cov=)}.
#' \code{exp} is a vector containing the exposure variable, \code{out} is a vector
#' containing the outcome variable and \code{cov} is a matrix containing the covariates
#' to be included in the test statistic.
#' @param prob Sampling probabilities from the \code{prop_score} function.
#' @param ins Matrix of instruments.
#' @param sig Significance level (default is 0.05).
#' @param nnodes Number of cores to use during parallel computing (default is 1).
#' @param out Vector of strings indicating which objects to return from the function.
#' Options include: "pvalues" for p-values corresponding to \code{beta}; "ci" for a \code{sig}-level
#' confidence interval obtained by inverting the test; "tobs" for the observed test statistic from
#' \code{test_stat}; and "tnull" for the full vector or matrix of counterfactual
#' test statistics (default is to return \code{c("pvalues","ci","tobs")}).
#' @param verbose Choose whether to show the progress bar (default is \code{TRUE})
#'
#' @return Named list of objects requested in \code{out}.
#'
#' @importFrom stats lm
#' @importFrom parallel clusterExport makeCluster
#' @importFrom pbapply parSapply pblapply
#'
#' @examples
#'
run_test <- function(reps, beta=0, dat, prob, ins, sig=0.05, nnodes=1, out=c("pvalues","ci","tobs"), verbose=T) {
  q <- ncol(prob$m)
  N <- nrow(prob$m)

  tobs <- test_stat(beta=beta,dat=dat,ins=ins)

  if (nnodes == 1) {
    cl = NULL
  }

  if (nnodes > 1) {
    if (get_os() == "windows") {
      cl <- makeCluster(min(reps, nnodes), outfile = "")
      clusterExport(cl, c("sampler","test_stat","dat","prob","beta","N","q","reps"),
                    envir = environment())
    } else {
      cl <- nnodes
    }
  }

  null_dist <- function(j) {
    Gres <- sapply(1:q, function(k) sampler(N,prob$m[,k]) + sampler(N,prob$f[,k]))
    return(test_stat(beta,dat,Gres))
  }

  if(verbose) {
    tnull <- matrix(t(pbsapply(X=1:reps, cl = cl, simplify=T, FUN=null_dist)), ncol=length(beta))
  } else {
    tnull <- matrix(t(parSapply(X=1:reps, cl = cl, simplify=T, FUN=null_dist)), ncol=length(beta))
  }

  if (nnodes > 1) {
    if (get_os() == "windows") {
      stopCluster(cl)
    }
  }

  pvalues <- sapply(1:length(beta), function(x) 1-mean(tobs[x] >= tnull[,x]))

  if("ci" %in% out) {
    id <- which(pvalues > sig)
    ci <- c(-Inf,Inf)
    if(pvalues[1] < sig) {
      ci[1] <- beta[min(id)]
    } else {
      warning('Lower bound of CI is ill-defined.')
    }

    if(pvalues[length(pvalues)] < sig) {
      ci[2] <- beta[max(id)]
    } else {
      warning('Upper bound of CI is ill-defined.')
    }
  }

  main <- sapply(X=out, USE.NAMES=T, simplify=F, FUN=function(x) assign(x,eval(parse(text=x))))
  return(main)
}



