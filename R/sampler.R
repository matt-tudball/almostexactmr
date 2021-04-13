#' @name sampler
#' @description Samples a counterfactual offspring haplotype
#'
#' @param N Sample size
#' @param prob Sampling probabilities from \code{prop_score} function
#'
#' @return Vector of offspring haplotypes
#'
#' @importFrom stats rbinom
#'
#' @examples
#'
sampler <- function(N,prob) {
  Z <- rbinom(N,1,prob)
  return(Z)
}
