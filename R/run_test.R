#' @name run_test
#' @title Main randomisation inference test
#' @description For a set of null hypotheses, computes the observed test statistic,
#' null distribution of test statistics, p-values and confidence intervals
#'
#' @param reps Number of resamples
#' @param beta Vector or scalar of null hypotheses (default is 0)
#' @param dat List of data
#' @param prob Sampling probabilities from \code{prop_score} function
#' @param ins Matrix of instruments
#' @param sig Significance level (default is 0.05)
#' @param nnodes Number of cores to use during parallel computing (default is 1)
#' @param out List of objects to return from the function (default is main statistics)
#' @param verbose Choose whether to show the progress bar (default is TRUE)
#'
#' @return F-statistic for the instruments in the linear model
#'
#' @importFrom stats lm
#' @importFrom parallel makeCluster clusterExport
#' @importFrom pbapply pblapply
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



