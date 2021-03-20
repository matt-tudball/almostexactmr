# Offspring haplotype (Markov chain model)
conditional_sampler <- function(Prob) {
  # Sampling probabilities
  PrM1 <- Prob$PrM1
  PrF1 <- Prob$PrF1

  # Generate offspring haplotype
  Zm <- rbinom(N,1,PrM1)
  Zf <- rbinom(N,1,PrF1)
  
  return(list(Zm=Zm,Zf=Zf))
}
