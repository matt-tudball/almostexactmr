#' @name test_stat
#' Compute test statistic
#'
#' @param beta Vector or scalar of null hypotheses
#' @param dat List of data
#' @param ins Matrix of instruments
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
  n <- nrow(y)

  # Compute residual sum of squares
  con <- matrix(rep(1,n),ncol=1)
  x1 <- cbind(con,cov)
  p1 <- ncol(x1)
  rss1 <- colSums(as.matrix(lm.fit(x=x1,y=y)$residuals^2))

  x2 <- cbind(con,ins,cov)
  p2 <- ncol(x2)

  rss2 <- colSums(as.matrix(lm.fit(x=x2,y=y)$residuals^2))

  # F statistic
  Fstat <- ((rss1 - rss2)/(p2 - p1))/(rss2/(n - p2))

  return(Fstat)
}
