rm(list=ls())
require(ivmodel); require(car); require(ggplot2); require(mvtnorm)

setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents")

# ---- SIMULATION PARAMETERS ---- #
# Choose variants to condition on
# Bias from pleiotropy by linkage
#Jset <- list(c(22,25,28),c(47,50,53),c(72,75,78),c(97,100,103),c(122,125,128))
# Correct
Jset <- list(c(23,25,27),c(48,50,52),c(73,75,77),c(98,100,102),c(123,125,127))
# No power
#Jset <- list(c(24,25,26),c(49,50,51),c(74,75,76),c(99,100,101),c(124,125,126))

# Null hypotheses to test
nullvec <- 0 #seq(-2,2,0.5)
lnull <- length(nullvec)

# Number of counterfactuals
lcf <- 2e3

# ---- Load package ---- #
load_all(path='FAMMR_FILES/CODE/almostexactmr')

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
a <- qnorm(1-0.4); b <- qnorm(1-0.05) # Lower and upper values of threshold
for(type in c('Mm','Mf','Fm','Ff')) {
  sigma <- 0.75^abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
  latent <- rmvnorm(N, mean = rep(0, p), sigma = sigma)
  threshold <- runif(p, a, b)
  out <- matrix(unlist(lapply(1:p, function(x) {1*(latent[,x] > threshold[x])})),ncol=p,byrow=F)
  assign(paste(type,'mat',sep=''), out)
}
MHap <- list(m = Mmmat, f = Mfmat)
FHap <- list(m = Fmmat, f = Ffmat)
rm(sigma, latent, threshold, out, Mmmat, Mfmat, Fmmat, Ffmat)

# Genetic data
OHap <- unconditional_sampler(MHap, FHap, p, d, epsilon=1e-8)
Z <- OHap$m+OHap$f

# Genetic instruments
MFHapG <- list(m = MHap$m[,Jz] + MHap$f[,Jz], f = FHap$m[,Jz] + FHap$f[,Jz])

# Unobserved offspring confounder
C <- rnorm(N,0,1)
ac <- sqrt(0.075); bc <- sqrt(0.075)

# Unobserved dynastic confounder
meansum <- 2*p*(1-integrate(pnorm,a,b)$value/(b-a)) # Mean
Cm <- rnorm(N,(rowSums(MHap$m+MHap$f)-meansum)/p,1)
Cf <- rnorm(N,(rowSums(FHap$m+FHap$f)-meansum)/p,1)
am <- sqrt(0.03); af <- sqrt(0.03); bm <- sqrt(0.03); bf <- sqrt(0.03)

# Exposure
D0 <- am*Cm + af*Cf + ac*C + rnorm(N,0,sqrt(0.7))
bd <- 0 # No effect of the exposure on the outcome

# Outcome
Y0 <- bm*Cm + bf*Cf + bc*C + rnorm(N,0,sqrt(0.7))

# Junk clean up
rm(a,ac,af,am,b,bc,bf,bm,C,Cf,Cm,Jy,meansum,type)

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
  W <- NULL
  #W <- cbind(Prob$m, Prob$f)
  #W <- cbind(MFHapG$m, MFHapG$f)
  #W <- cbind(H, MFHapH$m, MFHapH$f)

  # ---- Compute p-value ----
  results <- run_test(reps=2e3, beta=nullvec, dat=list(out=Y,exp=D,cov=W), prob=Prob,
                      ins=G, nnodes=4, out=c("pvalues"), verbose=F)
  return(results$pvalues)
}))

hist(out)

# ---- Save simulation results ---- #
#saveRDS(out, file='FAMMR_FILES/DATA/test_statistic_sim_correct.rds')

# ---- Create plots ---- #
out <- data.frame(out)
val <- expand.grid(1:lnull, adjset)
for(k in 1:nrow(val)) {
  v <- paste('n',val[k,1],'a',val[k,2],sep='')
  plot <- ggplot(out, aes_string(x=v)) +
    geom_histogram(color="darkblue", fill="lightblue",bins=20,binwidth=0.05,center=0.025) +
    xlab("p-value") + ylab("Count")
  plot
  ggsave(filename=paste("FAMMR_FILES/FIGURES/pvalue_",v,type=".pdf",sep=""),plot=plot,width=4,height=3)
}
