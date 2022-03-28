rm(list=ls())
setwd("Z:/projects/ieu2/p1/096/working/data")

library(devtools); library(ivmodel)

# Phased offspring haplotypes
OHap <- readRDS("haps/child_haps_chr21.rds")

# Phased maternal haplotypes
MHap <- readRDS("haps/mother_haps_chr21.rds")

# Phenotype data
pheno <- readRDS("pheno/pheno_formatted.rds")

# Genetic map
map <- readRDS("haps/map_chr21.rds")

# Load almostexactmr package
load_all(path='C:/Users/ow18301/OneDrive - University of Bristol/Documents/FAMMR_FILES/CODE/almostexactmr')

# Format phenotype data
pheno <- pheno[,colnames(pheno) %in% c("cidB3744","qlet","f7ms026a","dw042","geneticid")]

bool <- (!(is.na(pheno$f7ms026a) | pheno$f7ms026a < 0)) &
        (!(is.na(pheno$dw042) | pheno$dw042 < 0))

subset_hap <- function(dat, bool) {
  if (nrow(dat) != length(bool)) stop("Lengths do not match...")

  dat$id <- dat$id[bool]
  dat$m <- dat$m[bool,]
  dat$f <- dat$f[bool,]

  return(dat)
}

pheno <- pheno[bool, ]
OHap$id <- OHap$id[bool]
MHap$id <- MHap$id[bool]
OHap$m <- OHap$m[bool,]
OHap$f <- OHap$f[bool,]
MHap$m <- MHap$m[bool,]
MHap$f <- MHap$f[bool,]

# Print data
rsids <- c("rs571312")
bool <- colnames(OHap$m) %in% rsids
print(data.frame(OHap$m[1:3,bool], MHap$m[1:3,bool], MHap$f[1:3,bool], pheno$f7ms026a[1:3], pheno$dw042[1:3]))

# MR in unrelated individuals
ols_model <- lm(dw042 ~ f7ms026a, data=pheno)
summary(ols_model)

iv <- OHap$m[,bool] + OHap$f[,bool]
iv_model <- ivmodel(Y = pheno$dw042, D=pheno$f7ms026a, Z=iv)
summary(iv_model)

# Propensity scores
prob <- prop_score(
  PHap=MHap,
  CHap=OHap$m,
  map=map,
  rsid=rsids,
  mb=1)

# Obtain results
results <- run_test(
  reps=1e5,
  beta=0,
  dat=list(out=pheno$dw042, exp=pheno$f7ms026a, cov=prob),
  prob=prob,
  rsid=rsids,
  nnodes=4,
  out=c("pvalues", "tobs", "tnull"),
  verbose=T)

print(results$pvalues)

# Plot null distribution
TStat_obs <- results$tobs
TStat_null <- data.frame(x = results$tnull[results$tnull<7.5])
cond <- (TStat_obs >= TStat_null$x)
plot <- ggplot(TStat_null, aes(x=x)) +
  geom_histogram(data=TStat_null, color="darkred",
                 fill="#FF9999", bins=50) +
  geom_histogram(data=subset(TStat_null,cond==TRUE), color="darkblue",
                 fill="lightblue", bins=50) +
  xlab("Test statistic") + ylab("Count")
plot

ggsave(filename="~/FAMMR_FILES/FIGURES/alspac_pvalue.pdf", plot=plot, width=5, height=3)

# Obtain p-value curve
beta <- seq(-5,10,0.01)
conf <- run_test(reps=1e4,
                 beta=beta,
                 dat=list(out=pheno$dw042,exp=pheno$f7ms026a,cov=prob),
                 prob=prob,
                 rsid=rsids,
                 nnodes=4,
                 out=c("pvalues"),
                 verbose=T)

data <- data.frame(x=beta, y=conf$pvalues)
plot <- ggplot(data) +
  geom_line(aes(x=x, y=y), color="darkblue") +
  geom_hline(yintercept=0.05, color="#FF9999") +
  xlab("Test statistic") + ylab("Count")
plot
