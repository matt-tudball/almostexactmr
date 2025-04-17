rm(list=ls())
require(ivmodel); require(car); require(ggplot2); require(mvtnorm); require(devtools); require(pbapply); require(bindata); require(parallel)

setwd("/home/matt-tudball/code/almostexactmr")

# Function to generate correlated Bernoulli draws
generate_correlated_binom <- function(N, p, Sigma) {
  # N     : number of samples
  # p     : vector of marginal probabilities (length k)
  # Sigma : k x k target correlation matrix
  
  k <- length(p)
  if (!all(dim(Sigma) == c(k,k)))
    stop("Sigma must be a square matrix with dimension = length(p)")
  
  X <- rmvbin(n = N,
              margprob = p,
              bincorr  = Sigma)
  
  colnames(X) <- paste0("V", seq_len(k))
  return(X)
}

# ---- Load almostexactmr package ---- #
load_all()

# Overwrite functions with old ones
source('simulation/old-functions/unconditional_sampler.R')
source('simulation/old-functions/prop_score_sim.R')
source('simulation/old-functions/run_test_sim.R')
source('simulation/old-functions/sampler.R')

# ---- SIMULATION PARAMETERS ---- #
# Choose variants to condition on
# Bias from pleiotropy by linkage
#Jset <- list(c(22,25,28),c(47,50,53),c(72,75,78),c(97,100,103),c(122,125,128))
# Correct
Jset <- list(c(23,25,27),c(48,50,52),c(73,75,77),c(98,100,102),c(123,125,127))
# No power
#Jset <- list(c(24,25,26),c(49,50,51),c(74,75,76),c(99,100,101),c(124,125,126))

# Null hypotheses to test
nullvec <- seq(-1.5,2,0.05)

# Number of counterfactuals
lcf <- 1

# ---- Generate the genetic data ----
# Sample size
N <- 1.5e4

# Parental haplotypes
p <- 150 # Number of sites

# Set of instruments
Jz <- c(25,50,75,100,125)
q <- length(Jz)
az <- rep(0,p); az[Jz-1] <- sqrt(0.5/q)

# Pleiotropic instruments
Jy <- c(23,27,48,52,73,77,98,102,123,127)
bz <- rep(0,p); bz[Jy] <- sqrt(0.5/length(Jy))

# Morgan distance between two SNPs
d <- runif(p-1,0,0.75)
d[c(37,62,86,112)] <- Inf # These SNPs are independent (e.g. different chromosome)

# ---- Choose variants to condition on ---- #
# Pleiotropy
#Jset <- list(c(22,25,28),c(47,50,53),c(72,75,78),c(97,100,103),c(122,125,128))
# Correct
Jset <- list(c(23,25,27),c(48,50,52),c(73,75,77),c(98,100,102),c(123,125,127))
# No power
#Jset <- list(c(24,25,26),c(49,50,51),c(74,75,76),c(99,100,101),c(124,125,126))

# ---- Generate the simulation data ----
# Generate the parental haplotypes
# Instructions for using 1000 Genomes in place of simulated data
# 1. Replace Sigma with an LD matrix from 1000 Genomes
# 2. Replace maf with a vector of corresponding effect allele frequencies from 1000 Genomes. The ref/alt alleles 
# must match those in the LD matrix
# 3. Run generate_correlated_binom with the LD matrix and maf vector
maf <- runif(p, 0.45, 0.5)
Sigma <- 0.7^abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
for(type in c('Mm','Mf','Fm','Ff')) {
  out <- generate_correlated_binom(N, maf, Sigma)
  assign(paste(type,'mat',sep=''), out)
}
MHap <- list(m = Mmmat, f = Mfmat)
FHap <- list(m = Fmmat, f = Ffmat)
rm(Sigma, out, Mmmat, Mfmat, Fmmat, Ffmat)

# Genetic data
OHap <- unconditional_sampler(MHap, FHap, p, d, epsilon=1e-8)
Z <- OHap$m+OHap$f

# Genetic instruments
MFHapG <- list(m = MHap$m[,Jz] + MHap$f[,Jz], f = FHap$m[,Jz] + FHap$f[,Jz])

# Unobserved offspring confounder
C <- rnorm(N,0,1)
ac <- sqrt(0.075); bc <- sqrt(0.075)

# Unobserved dynastic confounder
Cm <- rnorm(N,(rowSums(MHap$m+MHap$f)-p * 0.475)/p,1)
Cf <- rnorm(N,(rowSums(FHap$m+FHap$f)-p * 0.475)/p,1)
am <- sqrt(0.075); af <- sqrt(0.075); bm <- sqrt(0.075); bf <- sqrt(0.075)

# Exposure
D0 <- am*Cm + af*Cf + ac*C + rnorm(N,0,sqrt(0.61))
bd <- 0 # No effect of the exposure on the outcome

# Outcome
Y0 <- bm*Cm + bf*Cf + bc*C + rnorm(N,0,sqrt(0.61))

# Junk clean up
rm(ac,af,am,bc,bf,bm,C,Cf,Cm,Jy,type)

# ---- SIMULATION BEGINS HERE ---- #
out <- t(pbsapply(X=1:lcf, cl=NULL, simplify=T, FUN=function(k) {
  # ---- Resample offspring genotype ---- #
  # Genetic data
  OHap <- unconditional_sampler(MHap, FHap, p, d, epsilon=1e-8)
  Z <- OHap$m+OHap$f

  # Genetic instruments
  G <- Z[,Jz]

  # ---- Observed exposure and outcome ---- #
  D <- D0 + Z%*%az
  Y <- Y0 + bd*D + Z%*%bz

  # ---- Conditional sampling probabilities ---- #
  Prob <- list(m = prop_score(MHap, OHap$m, Jset, d),
               f = prop_score(FHap, OHap$f, Jset, d))

  # ---- Choose adjustment set ---- #
  #W <- NULL
  W <- cbind(Prob$m, Prob$f)
  #W <- cbind(MFHapG$m, MFHapG$f)
  #W <- cbind(H, MFHapH$m, MFHapH$f)

  # ---- Compute p-value ----
  results <- run_test(reps=2e3, beta=nullvec, dat=list(out=Y,exp=D,cov=W), prob=Prob,
                      ins=G, nnodes=4, out=c("pvalues"), verbose=TRUE)
  return(results$pvalues)
}))

# ---- Save simulation results ---- #
saveRDS(out, file='FAMMR_FILES/DATA/power_curve_2.rds')

# ---- Create plots ---- #
out <- data.frame(out)
for(j in 1:ncol(out)) {
  name <- colnames(out)[j]
  plot <- ggplot(out, aes_string(x=name)) +
    geom_histogram(color="darkblue", fill="lightblue",bins=20,binwidth=0.05,center=0.025) +
    xlab("p-value") + ylab("Count") + xlim(0,1)
  print(plot)
  ggsave(filename=paste("simulation/pvalue_correct_n",j,"a3",type=".pdf",sep=""),plot=plot,width=4,height=3)
  Sys.sleep(3)
}


