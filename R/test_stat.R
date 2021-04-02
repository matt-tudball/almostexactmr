#' Compute test statistic (Should use lm.fit and allow a vector of beta, see example.)
#'
#' @param beta null hypothesis
#' @param X Data frame (or maybe matrix?) of instruments
#' @param W Covariates to adjust in linear model
#'
#' @return F-statistic for the instruments in the linear model
#'
#' @importFrom stats lm
#' @importFrom car linearHypothesis
#'
#' @examples
#' n <- 15000
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- matrix(rnorm(n * 100), n, 100) # no intercept!
#' system.time(lm.fit(x = X, y = y[, 1]))
#' system.time(replicate(100, lm.fit(x = X, y = y[, 1])))
#'
test_stat <- function(beta, X, W = NULL) {
  if (is.null(W)) {
    Tdat <- data.frame(Yadj=Y-beta*D, X)
  } else {
    Tdat <- data.frame(Yadj=Y-beta*D, X, W)
  }
  model <- lm(Yadj~., data=Tdat)
  testvar <- colnames(Tdat)[startsWith(colnames(Tdat),'G')]
  hypo <- sapply(testvar, function(x) {paste(x,'=0',sep='')})
  return(linearHypothesis(model, hypo)$F[2])
}
