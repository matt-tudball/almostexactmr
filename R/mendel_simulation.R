rm(list=ls())
require(ivmodel); require(car); require(ggplot2); require(mvtnorm)

setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents")

# ---- SIMULATION PARAMETERS ---- #
# Choose adjustment set
# 1: no W
# 2: W = Probdat
# 3: W = MFHapG
adjset <- c(1,2,3); ladj <- length(adjset)

# Choose variants to condition on
# Bias from pleiotropy by linkage
#Jset <- list(c(22,25,28),c(47,50,53),c(72,75,78),c(97,100,103),c(122,125,128))
# Correct
Jset <- list(c(23,25,27),c(48,50,52),c(73,75,77),c(98,100,102),c(123,125,127))
# No power
#Jset <- list(c(24,25,26),c(49,50,51),c(74,75,76),c(99,100,101),c(124,125,126))

# Null hypotheses to test
nullset <- c(-0.3,0,0.3,0.6); lnull <- length(nullset)

# Number of permutations per test
lperm <- 2000

# Number of counterfactuals
lcf <- 100

# ---- Set seed ---- #
set.seed(1215)

# ---- Function for naming data frames ---- #
# Function for naming data frames
namer <- function(D, namelist, indexlist) {
  names <- NULL
  k <- length(namelist)
  p <- ncol(D)/k
  for(i in 1:p) {
    for(j in 1:k) {
      names <- c(names,paste(namelist[j],'_',indexlist[i],'_',sep=''))
    }
  }
  colnames(D) <- names
  return(D)
}

# ---- Function for computing forward-backward weights  ----
source('FAMMR_FILES/CODE/almost-exact-mr/src/forward_backward_weights.R')

# ---- Function for computing sampling probabilities  ----
source('FAMMR_FILES/CODE/almost-exact-mr/src/sampling_probability.R')

# ---- Function for conditionally resampling offspring haplotypes  ----
source('FAMMR_FILES/CODE/almost-exact-mr/src/conditional_sampler.R')

# ---- Function for unconditionally resampling offspring haplotypes  ----
source('FAMMR_FILES/CODE/almost-exact-mr/src/unconditional_sampler.R')

# ---- Function for computing the test statistic ----
source('FAMMR_FILES/CODE/almost-exact-mr/src/test_statistic.R')


# ---- Generate the genetic data ----
# Sample size
N <- 15000

# Probability of a de novo mutation
epsilon <- 1e-8

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
Mmmat <- namer(Mmmat,'Mm',1:p)
Mfmat <- namer(Mfmat,'Mf',1:p)
Fmmat <- namer(Fmmat,'Fm',1:p)
Ffmat <- namer(Ffmat,'Ff',1:p)
MFHap <- data.frame(Mmmat, Mfmat, Fmmat, Ffmat)
rm(sigma, latent, threshold, out, Mmmat, Mfmat, Fmmat, Ffmat)

# Unobserved offspring confounder
C <- rnorm(N,0,1)
ac <- sqrt(0.075); bc <- sqrt(0.075)

# Unobserved dynastic confounder
meansum <- 2*p*(1-integrate(pnorm,a,b)$value/(b-a)) # Mean
Cm <- rnorm(N,(rowSums(MFHap[,startsWith(colnames(MFHap),'M')])-meansum)/p,1)
Cf <- rnorm(N,(rowSums(MFHap[,startsWith(colnames(MFHap),'F')])-meansum)/p,1)
am <- sqrt(0.03); af <- sqrt(0.03); bm <- sqrt(0.03); bf <- sqrt(0.03)

# Exposure
D0 <- am*Cm + af*Cf + ac*C + rnorm(N,0,sqrt(0.7))
bd <- 0 # No effect of the exposure on the outcome

# Outcome
Y0 <- bm*Cm + bf*Cf + bc*C + rnorm(N,0,sqrt(0.7))

# Junk clean up
rm(a,ac,af,am,b,bc,bf,bm,C,Cf,Cm,Jy,meansum,type)

# ---- SIMULATION BEGINS HERE ---- #
out <- lapply(1:ladj, function(X) matrix(ncol=lnull, nrow=lcf))
names(out) <- as.vector(paste('adj',adjset,sep=''))

for(ja in 1:ladj) {
  for(jn in 1:lnull) {
    for(jc in 1:lcf) {
      # ---- Resample offspring genotype ---- #
      # Genetic data
      OHap <- unconditional_sampler(MFHap, p, d, epsilon)
      OHap$m <- namer(OHap$m,'Zm',1:p); OHap$f <- namer(OHap$f,'Zf',1:p)
      Z <- OHap$m+OHap$f; Z <- namer(Z,'Z',1:p)
      
      # Genetic instruments
      G <- Z[,Jz]; G <- namer(G,'G',Jz)
      MFHapG <- MFHap[,grepl(paste('_',Jz,'_', sep='', collapse = "|"), colnames(MFHap))]
      
      # ---- Observed exposure and outcome ---- #
      D <- D0 + Z%*%az
      Y <- Y0 + bd*D + Z%*%bz
      
      # ---- Forward and backward weights ---- #
      W <- fbweights(MFHap, OHap, p, d, epsilon)
      
      # ---- Conditional sampling probabilities ---- #
      Prob <- vector(mode="list", length=q)
      for(k in 1:q) {
        Prob[[k]] <- sampling_probability(MFHap, W, Jset[[k]], d, epsilon)
      }
      Probdat <- data.frame(matrix(unlist(Prob), nrow=N, byrow=FALSE))
      Probdat <- namer(Probdat,c('Pm','Pf'),Jz)
      
      # ---- Generate test statistic distribution for one hypothesis ----
      # Null hypothesis
      beta0 <- nullset[jn]
      
      # Test statistic is fixed effect regression of Y-beta0*D on Z
      if(adjset[ja] == 1) W <- NULL
      if(adjset[ja] == 2) W <- Probdat
      if(adjset[ja] == 3) W <- MFHapG
      
      TStat_obs <- tstat_calc(beta0,X=G,W=W)
      
      TStat_null <- tdist_calc(R=lperm,beta0,Prob,W,Jz,verbose=T)
      
      # p-value
      print(paste('On counterfactual',jc,'testing null of',nullset[jn],
                  'using adjustment set',adjset[ja]))
      out[[ja]][jc,jn] <- 1-mean(TStat_obs >= TStat_null$beta0)
    }
  }
}

# ---- Save simulation results ---- #
saveRDS(out, file='FAMMR_FILES/DATA/test_statistic_sim.rds')