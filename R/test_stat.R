#' @name test_stat
#' @title Test statistic
#' @description Computes the test statistic as the F-statistic for the instruments
#' from a linear regression of the adjusted outcome on the instruments and a user-specified
#' set of covariates.
#'
#' @param beta Vector or scalar of null hypotheses to be tested.
#' @param dat List of data which must be of the form \code{list(exp=, out=, cov=)}.
#' \code{exp} is a vector containing the exposure variable, \code{out} is a vector
#' containing the outcome variable and \code{cov} is a matrix containing the covariates
#' to be included in the test statistic.
#' @param ins Matrix of instruments.
#'
#' @return F-statistic for the instruments in the linear model
#'
#' @importFrom stats lm
#'
#' @examples
#'
test_stat <- function(beta, dat, ins) {
  out <- dat$out; exp <- dat$exp; cov <- dat$cov
  # Matrix of outcome variables
  y <- sapply(beta, function(x) out - x*exp)
  m <- nrow(y)

  # Compute residual sum of squares
  con <- matrix(rep(1,m),ncol=1)
  x1 <- cbind(con,cov)
  p1 <- ncol(x1)
  rss1 <- colSums(as.matrix(lm.fit(x=x1,y=y)$residuals^2))

  x2 <- cbind(con,ins,cov)
  p2 <- ncol(x2)

  rss2 <- colSums(as.matrix(lm.fit(x=x2,y=y)$residuals^2))

  # F statistic
  Fstat <- ((rss1 - rss2)/(p2 - p1))/(rss2/(m - p2))

  return(Fstat)
}
