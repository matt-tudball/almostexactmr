# Offspring haplotype (Markov chain model)
conditional_sampler <- function(Prob) {
  # Generate offspring haplotype
  Z <- rbinom(N,1,Prob)
  
  return(Z)
}
