# Create function for test statistic
test_stat <- function(beta0, X, W = NULL) {
  if (is.null(W)) { 
    Tdat <- data.frame(Yadj=Y-beta0*D, X) 
  } else {
    Tdat <- data.frame(Yadj=Y-beta0*D, X, W)
  }
  model <- lm(Yadj~., data=Tdat)
  testvar <- colnames(Tdat)[startsWith(colnames(Tdat),'G')]
  hypo <- sapply(testvar, function(x) {paste(x,'=0',sep='')})
  return(linearHypothesis(model, hypo)$F[2])
}
