#' @name fbweights
#' @title Forward-backward algorithm
#' @description Computes forward-backward weights via the forward-backward algorithm
#'
#' @param PHap Named list of parental haplotypes which must be of the form
#' \code{list(m=, f=)}.\code{m} is a (0,1)-matrix of maternally-inherited haplotypes and
#' \code{f} is a (0,1)-matrix of paternally-inherited haplotypes. Each row corresponds
#' to an offspring and each column is a SNP.
#' @param CHap (0,1)-matrix of offspring haplotypes. Each row corresponds
#' to an offspring and each column is a SNP.
#' @param d Vector of genetic distances corresponding to haplotypes in PHap and CHap
#' @param epsilon De novo mutation rate (default is 1e-8)
#'
#' @return Named list of forward-backward weights
#'
#' @examples
#'
fbweights <- function(PHap, CHap, d, epsilon=1e-8) {
  # Function for the probability of inheriting a SNP given meiosis indicator
  prsnp <- function(P,Q) {
    return(ifelse(P == Q, 1-epsilon, epsilon))
  }

  # Check equality of rows and columns
  if(!all(sapply(list(PHap$m, PHap$f, CHap), function(x) nrow(x) == nrow(MHap$m))))
    stop('Rows are not equal')
  if(!all(sapply(list(PHap$m, PHap$f, CHap), function(x) ncol(x) == ncol(MHap$m))))
    stop('Columns are not equal')

  # Parental haplotypes
  Pm <- PHap$m; Pf <- PHap$f

  # Sample size
  N <- nrow(CHap)

  # Number of sites
  p <- ncol(CHap)

  # Construct forward weights
  Fweight <- list(m = matrix(nrow=N, ncol=(p+1)),
                  f = matrix(nrow=N, ncol=(p+1)))

  # Zero forward weight
  Fweight$m[,1] <- 0.5; Fweight$f[,1] <- 0.5

  # Define forward weights recursively
  for (j in 2:(p+1)) {
    dat <- list(m = Pm[,j-1], f = Pf[,j-1], o = CHap[,j-1])

    # Initial forward weight
    if (j == 2) {
      Fweight$m[,j] <- 0.5*prsnp(dat$m,dat$o)
      Fweight$f[,j] <- 0.5*prsnp(dat$f,dat$o)
    }

    # Define remaining forward weights recursively
    if (j > 2) {
      # Forward weights one place behind
      fback <- list(m = Fweight$m[,j-1], f = Fweight$f[,j-1])

      # Recombination probability
      r_prob <- 0.5*(1+exp(-2*d[j-2]))

      # Additional forward weights
      Fweight$m[,j] <- prsnp(dat$m,dat$o)*(r_prob*fback$m + (1-r_prob)*fback$f)
      Fweight$f[,j] <- prsnp(dat$f,dat$o)*((1-r_prob)*fback$m + r_prob*fback$f)
    }
  }

  # Construct backward weights
  Bweight <- list(m = matrix(nrow=N, ncol=p),
                  f = matrix(nrow=N, ncol=p))

  # Define backward weights recursively
  for(j in p:1) {
    # Initial backward weight
    if(j == p) {
      Bweight$m[,j] <- 1
      Bweight$f[,j] <- 1
    }

    # Define remaining backward weights recursively
    if(j < p) {
      dat <- list(m = Pm[,j+1], f = Pf[,j+1], o = CHap[,j+1])

      # Backward weights one place ahead
      bfor <- list(m = Bweight$m[,j+1], f = Bweight$f[,j+1])

      # Recombination probability
      r_prob <- 0.5*(1+exp(-2*d[j]))

      # Additional backward weights
      Bweight$m[,j] <- bfor$m*r_prob*prsnp(dat$m,dat$o) + bfor$f*(1-r_prob)*prsnp(dat$f,dat$o)
      Bweight$f[,j] <- bfor$m*(1-r_prob)*prsnp(dat$m,dat$o) + bfor$f*r_prob*prsnp(dat$f,dat$o)
    }
  }
  return(list(b=Bweight, a=Fweight))
}
