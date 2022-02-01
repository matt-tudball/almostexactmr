#' @name prop_score
#' @title Propensity scores
#' @description Computes propensity scores for the instruments
#'
#' @param PHap Named list of parental haplotypes which must be of the form
#' \code{list(m=, f=)}.\code{m} is a (0,1)-matrix of maternally-inherited haplotypes and
#' \code{f} is a (0,1)-matrix of paternally-inherited haplotypes. Each row corresponds
#' to an offspring and each column is a SNP.
#' @param CHap (0,1)-matrix of offspring haplotypes. Each row corresponds
#' to an offspring and each column is a SNP.
#' @param Jset List of indices for boundary and instrument SNPs.
#' @param d Vector of genetic distances corresponding to haplotypes in PHap and CHap
#' @param epsilon De novo mutation rate (default is 1e-8)
#'
#' @return Matrix of propensity scores
#'
#' @examples
#'
prop_score <- function(PHap, CHap, map, rsid, mb, epsilon=1e-8) {
  # Genetic distances
  d <- map$dist[-1]

  # Forward-backward weights
  W <- fbweights(PHap, CHap, d, epsilon)

  # Find indices for instrument and upper/lower edge of window
  j <- which(colnames(CHap) == rsid)
  pos <- map$snppos[map$rsid == rsid]

  lrsid <- map$rsid[max(which(map$snppos <= pos-mb*1e6))]
  l <- which(colnames(CHap) == lrsid)

  hrsid <- map$rsid[min(which(map$snppos >= pos+mb*1e6))]
  h <- which(colnames(CHap) == hrsid)

  # Probabilities for the meiosis indicators
  r_probl <- 0.5*(1+exp(-2*sum(d[l:(j-1)])))
  if(h == j+1) { r_probh = 1 } else { r_probh <- 0.5*(1+exp(-2*sum(d[j:(h-2)]))) }

  # Backward weights
  beta <- list(m = W$b$m[,h-1], f = W$b$f[,h-1])

  # Forward weights
  alpha <- list(m = W$a$m[,l+1], f = W$a$f[,l+1])

  # Meiosis indicator distribution for parental haplotype
  propm <- (r_probh*beta$m + (1-r_probh)*beta$f)*(r_probl*alpha$m + (1-r_probl)*alpha$f)
  propf <- ((1-r_probh)*beta$m + r_probh*beta$f)*((1-r_probl)*alpha$m + r_probl*alpha$f)
  PrPm <- propm/(propm+propf)

  # Probabilities for the offspring haplotypes
  Pm <- PHap$m; Pf <- PHap$f
  PrP1 <- ifelse(Pm[,j] == 1, 1-epsilon, epsilon)*PrPm +
          ifelse(Pf[,j] == 1, 1-epsilon, epsilon)*(1-PrPm)

  return(PrP1)
}
