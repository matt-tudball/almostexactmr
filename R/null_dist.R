# Create function for calculating null distribution
null_dist <- function(R, beta, Prob, W = NULL, Jz, verbose=T) {
  q <- length(Jz)
  TStat_dist <- data.frame(sapply(X=1:R, function(x) {
    if (verbose) {
      if(x%%100==0) cat('.')
      if(x%%1000==0) cat(paste(x/1000,'K',sep=''))
      if(x%%5000==0) cat('\n')
    }
    Gres <- matrix(nrow=N,ncol=q)
    # Construct resampled genetic instruments
    for(k in 1:q) {
      Gres[,k] <- conditional_sampler(Prob[[k]]$m) + conditional_sampler(Prob[[k]]$f)
    }
    Gres <- namer(Gres,'G',Jz)
    # Compute test statistic
    TStat <- test_stat(beta,X=Gres,W=W)
  }))
  colnames(TStat_dist) <- "beta"
  return(TStat_dist)
}
