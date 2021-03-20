# Create function for test statistic
tstat_calc <- function(beta0, X, W = NULL) {
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

# Create function for calculating null distribution
tdist_calc <- function(R, beta0, Prob, W = NULL, Jz, verbose=T) {
  q <- length(Jz)
  TStat_dist <- data.frame(sapply(X=1:R, function(x) {
    if (verbose & x%%500==0) print(paste('Iteration:',x))
    Gres <- matrix(nrow=N,ncol=q)
    # Construct resampled genetic instruments
    for(k in 1:q) {
      Gsample <- conditional_sampler(Prob[[k]])
      Gres[,k] <- Gsample$Zm + Gsample$Zf
    }
    Gres <- namer(Gres,'G',Jz)
    # Compute test statistic
    TStat <- tstat_calc(beta0,X=Gres,W=W)
  }))
  colnames(TStat_dist) <- "beta0"
  return(TStat_dist)
}
