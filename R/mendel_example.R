rm(list=ls())
require(ivmodel); require(car); require(ggplot2); require(mvtnorm)

setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents")

# ---- Set seed ---- #
set.seed(1210)

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
Mmmat <- namer(Mmmat,'Mm',1:p)
Mfmat <- namer(Mfmat,'Mf',1:p)
Fmmat <- namer(Fmmat,'Fm',1:p)
Ffmat <- namer(Ffmat,'Ff',1:p)
MFHap <- data.frame(Mmmat, Mfmat, Fmmat, Ffmat)
rm(sigma, latent, threshold, out)

# Genetic data
OHap <- unconditional_sampler(MFHap, p, d, epsilon)
OHap$m <- namer(OHap$m,'Zm',1:p); OHap$f <- namer(OHap$f,'Zf',1:p)
Z <- OHap$m+OHap$f; Z <- namer(Z,'Z',1:p)

# Genetic instruments
G <- Z[,Jz]; G <- namer(G,'G',Jz)
MFHapG <- MFHap[,grepl(paste('_',Jz,'_', sep='', collapse = "|"), colnames(MFHap))]

# Other genetic variants
Hset <- 1:p
for(k in 1:q) {
Hset <- setdiff(Hset, c((Jset[[k]][1]+1):(Jset[[k]][3]-1)))
}
H <- Z[,Hset]

MFHapH <- MFHap[,grepl(paste('_',Hset,'_', sep='', collapse = "|"), colnames(MFHap))]

# Unobserved offspring confounder
C <- rnorm(N,0,1)
ac <- sqrt(0.075); bc <- sqrt(0.075)

# Unobserved dynastic confounder
meansum <- 2*p*(1-integrate(pnorm,a,b)$value/(b-a)) # Mean
Cm <- rnorm(N,(rowSums(MFHap[,startsWith(colnames(MFHap),'M')])-meansum)/p,1)
Cf <- rnorm(N,(rowSums(MFHap[,startsWith(colnames(MFHap),'F')])-meansum)/p,1)
am <- sqrt(0.03); af <- sqrt(0.03); bm <- sqrt(0.03); bf <- sqrt(0.03)

# Exposure
D <- Z%*%az + am*Cm + af*Cf + ac*C + rnorm(N,0,sqrt(0.7))
bd <- 0 # No effect of the exposure on the outcome

# Outcome
Y <- bd*D +  Z%*%bz + bm*Cm + bf*Cf + bc*C + rnorm(N,0,sqrt(0.7))

# Junk clean up
rm(a,ac,af,am,az,b,bc,bf,bm,bz,C,Cf,Cm,Hset,Jy,k,meansum,type)

# ---- Forward and backward weights ---- #
MWgt <- fbweights(Mmmat, Mfmat, OHap$m, p, d, epsilon)
FWgt <- fbweights(Fmmat, Ffmat, OHap$f, p, d, epsilon)

# ---- Visualise the data ---- #

# Print the first 6 observations
round(data.frame(G[,'G_25_'], D, Y, MFHapG[,grepl(paste('_',25,'_', sep=''), colnames(MFHapG))])[1:6,],2)

# Resample once from the randomisation distribution
Gres <- matrix(nrow=N,ncol=q)
Prob <- vector(mode="list", length=q)
for(k in 1:q) {
  Prob[[k]] <- list(m = sampling_probability(Mmmat, Mfmat, MWgt, Jset[[k]], d, epsilon),
                f = sampling_probability(Fmmat, Ffmat, FWgt, Jset[[k]], d, epsilon))
  Gres[,k] <- conditional_sampler(Prob[[k]]$m) + conditional_sampler(Prob[[k]]$f)
}
Gres <- namer(Gres,'G',Jz)
Probdat <- data.frame(matrix(unlist(Prob), nrow=N, byrow=FALSE))
Probdat <- namer(Probdat,c('Pm','Pf'),Jz)

round(data.frame(Gres[,'G_25_'], D, Y, MFHapG[,grepl(paste('_',25,'_', sep=''), colnames(MFHapG))])[1:6,],2)

# ---- Super-population inference ----
# Unrelated
iv_unrelated <- ivmodel(Y=Y,D=D,Z=G)

print(paste('First stage F-statistic:', round(iv_unrelated$AR$Fstat,2)))

print(paste('TSLS 95% CI: [', round(iv_unrelated$kClass$ci['1',1],2),', ',
            round(iv_unrelated$kClass$ci['1',2],2),']', sep=''))

print(paste('AR 95% CI: [', round(iv_unrelated$AR$ci[1],2),', ',
            round(iv_unrelated$AR$ci[2],2),']', sep=''))

# Only conditioning on parental haplotypes
iv_parental <- ivmodel(Y=Y,D=D,Z=G,X=MFHapG)

print(paste('First stage F-statistic:', round(iv_parental$AR$Fstat,2)))

print(paste('TSLS 95% CI: [', round(iv_parental$kClass$ci['1',1],2),', ',
            round(iv_parental$kClass$ci['1',2],2),']', sep=''))

print(paste('AR 95% CI: [', round(iv_parental$AR$ci[1],2),', ',
            round(iv_parental$AR$ci[2],2),']', sep=''))

# Conditioning on full propensity score
iv_full <- ivmodel(Y=Y,D=D,Z=G,X=Probdat)

print(paste('First stage F-statistic:', round(iv_full$AR$Fstat,2)))

print(paste('TSLS 95% CI: [', round(iv_full$kClass$ci['1',1],2),', ',
            round(iv_full$kClass$ci['1',2],2),']', sep=''))

print(paste('AR 95% CI: [', round(iv_full$AR$ci[1],2),', ',
            round(iv_full$AR$ci[2],2),']', sep=''))

# ---- Generate test statistic distribution for one hypothesis ----
# Null hypothesis
beta0 <- -0.12

# Test statistic is fixed effect regression of Y-beta0*D on Z
#W <- NULL
W <- Probdat
#W <- H,MFHapH,MFHapG
#W <- MFHapG
TStat_obs <- tstat_calc(beta0,X=G,W=W)

TStat_null <- tdist_calc(R=5e4,beta0,Prob,W,Jz,verbose=T)

# p-value
print(1-mean(TStat_obs >= TStat_null$beta0))

# Plot null distribution
cond <- (TStat_obs >= TStat_null$beta0)
plot <- ggplot(TStat_null, aes(x=beta0)) +
  geom_histogram(data=TStat_null, color="darkred",
                 fill="#FF9999", bins=50) +
  geom_histogram(data=subset(TStat_null,cond==TRUE), color="darkblue",
                 fill="lightblue", bins=50) +
  xlab("Test statistic") + ylab("Count")
plot

ggsave(filename=paste("FAMMR_FILES/FIGURES/null_dist",type=".pdf",sep=""),plot=plot)
