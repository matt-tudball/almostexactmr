rm(list=ls())
require(ivmodel); require(car); require(ggplot2); require(mvtnorm); require(devtools); require(pbapply); require(sim1000G)

setwd()

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF(vcf_file, maxNumberOfVariants = 567 , min_maf = 0.01 , max_maf = NA)
readGeneticMap(chromosome = 4)

genetic_map_of_region = system.file("examples",
"chr4-geneticmap.txt",
package = "sim1000G")
readGeneticMapFromFile(genetic_map_of_region)

n_parents <- 200 # Must be even
n_children <- n_parents / 2 # One child per couple

startSimulation(vcf, totalNumberOfIndividuals = n_parents + n_children)

map <- data.frame(
  rsid = do.call(paste, c(SIM$varinfo[,c("#CHROM","POS","REF","ALT")], sep='-')),
  dist = c(NA, diff(SIM$cm))
) # Genetic distances between SNPs

id <- c()
for(i in 1:n_parents) { 
  id[i] <- SIM$addUnrelatedIndividual()
}

pairs <- split(id, rep(1:(length(id) / 2), each = 2)) # Pair up the parents
names(pairs) <- NULL

mother_id <- sapply(X = pairs, FUN = function(x) x[1])
father_id <- sapply(X = pairs, FUN = function(x) x[2])

MHap <- list(m = SIM$gt1[mother_id, ], f = SIM$gt2[mother_id, ])
FHap <- list(m = SIM$gt1[father_id, ], f = SIM$gt2[father_id, ])
colnames(MHap$m) <- colnames(MHap$f) <- colnames(FHap$m) <- colnames(FHap$f) <- map$rsid

for (i in 1:n_children) {
  id[n_parents + i] <- SIM$mate(pairs[[i]][1], pairs[[i]][2]) # Add children
  pairs[[i]] <- c(pairs[[i]], n_parents + i)
}
child_id <- sapply(X = pairs, FUN = function(x) x[3])

OHap <- list(m = SIM$gt1[child_id, ], f = SIM$gt2[child_id, ])
colnames(OHap$m) <- colnames(OHap$f) <- map$rsid

# ---- SIMULATION PARAMETERS ---- #
# Null hypotheses to test
nullvec <- seq(-1.5,2,0.05)

# Number of counterfactuals
lcf <- 1e3

# ---- Load package ---- #
load_all(path='/home/matt-tudball/code/almostexactmr')

# ---- Generate the genetic data ----
# Sample size
N <- n_children

# Parental haplotypes
p <- ncol(OHap$m) # Number of sites

# Set of instruments
Jz <- c(300)
q <- length(Jz)
az <- rep(0,p); az[Jz-1] <- sqrt(0.5/q)

# Pleiotropic instruments
Jy <- c(100)
bz <- rep(0,p); bz[Jy] <- sqrt(0.5/length(Jy))

# ---- Choose variants to condition on ---- #
Jset <- list(c(150, 300, 450))

# ---- Generate the simulation data ----
# Genetic data
Z <- OHap$m+OHap$f

# Genetic instruments
MFHapG <- list(m = MHap$m[,Jz] + MHap$f[,Jz], f = FHap$m[,Jz] + FHap$f[,Jz])

# Unobserved offspring confounder
C <- rnorm(N,0,1)
ac <- sqrt(0.075); bc <- sqrt(0.075)

# Unobserved dynastic confounder
meansum <- mean(rowSums(MHap$m + MHap$f)) # Mean
Cm <- rnorm(N, (rowSums(MHap$m + MHap$f) - meansum) / p, 1)
Cf <- rnorm(N, (rowSums(FHap$m + FHap$f) - meansum) / p, 1)
am <- sqrt(0.075); af <- sqrt(0.075); bm <- sqrt(0.075); bf <- sqrt(0.075)

# Exposure
D0 <- am*Cm + af*Cf + ac*C + rnorm(N,0,sqrt(1 - am^2 - af^2 - ac^2))
bd <- 0 # No effect of the exposure on the outcome

# Outcome
Y0 <- bm*Cm + bf*Cf + bc*C + rnorm(N,0,sqrt(1 - bm^2 - bf^2 - bc^2))

# Junk clean up
rm(ac,af,am,bc,bf,bm,C,Cf,Cm,Jy,meansum)

region <- data.frame(lower=map$rsid[150], snps=map$rsid[300], upper=map$rsid[450])

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
  Prob <- list(m = prop_score(MHap, OHap$m, map, region),
               f = prop_score(FHap, OHap$f, map, region))

  # ---- Choose adjustment set ---- #
  #W <- NULL
  W <- cbind(Prob$m, Prob$f)
  #W <- cbind(MFHapG$m, MFHapG$f)
  #W <- cbind(H, MFHapH$m, MFHapH$f)

  # ---- Compute p-value ----
  results <- run_test(reps = 2e3, beta = nullvec, OHap$m, MHap, pheno=list(out=Y,exp=D,cov=W), prob=Prob,
                      snps = region$snps, cores=4, out=c("pvalues"))
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
  ggsave(filename=paste("FAMMR_FILES/FIGURES/pvalue_correct_n",j,"a3",type=".pdf",sep=""),plot=plot,width=4,height=3)
  Sys.sleep(3)
}


