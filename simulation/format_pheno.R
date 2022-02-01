rm(list=ls())
setwd("Z:/projects/ieu2/p1/096/working/data")

# Phased offspring haplotypes
OHap <- readRDS("haps/child_haps.rds")

# Phenotype data
raw <- read.csv("pheno/B3744_18Jun2021.csv", header=T)

# Linker files (links phenotype ID with genetic ID)
link <- read.csv("pheno/B3744_datasetids.csv", header=T)

famid <- as.integer(substr(OHap$id,1,nchar(OHap$id)-1))
indid <- as.character(substr(OHap$id,nchar(OHap$id),nchar(OHap$id)))

keepid <- sapply(X=1:nrow(raw), FUN=function(i) {
  cid <- link$gi_1000g_g0m_g1[link$cidB3744 == raw$cidB3744[i]]
  bool <- (cid %in% famid) & (raw$qlet[i] %in% c("A","B","C","D"))
  if(bool) return(cid) else return(NA)
})

keep <- !is.na(keepid)
pheno <- raw[keep,]
pheno$geneticid <- keepid[keep]

rm(raw, keepid, keep)

pheno <- data.frame(pheno[match(famid, pheno$geneticid),])

saveRDS(pheno, file='pheno/pheno_formatted.rds')
