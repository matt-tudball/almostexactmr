#' @name prop_score
#' @title Propensity scores
#' @description Computes propensity scores for the instruments
#'
#' @param PHap Named list of mother or father haplotypes
#' @param CHap Matrix of offspring haplotypes
#' @param Jset List of indices for boundary and instrument SNPs
#' @param d Vector of genetic distances corresponding to haplotypes in PHap and CHap
#' @param epsilon De novo mutation rate (default is 1e-8)
#'
#' @return Matrix of propensity scores
#'
#' @examples
#'
prop_score <- function(PHap, CHap, Jset, d, epsilon=1e-8) {

  # Forward-backward weights
  W <- fbweights(PHap, CHap, d, epsilon)

  # Indices for lower bound, marker and upper bound respectively
  out <- sapply(Jset, function(x) {
    l <- x[1]; j <- x[2]; h <- x[3]

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
  })

  return(out)
}
