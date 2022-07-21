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
#' @param mb Distance in megabases from each instrument to start conditioning on chromosome.
#' @param map Map file of genetic distances between each adjacent SNP. Must contain a
#' variable called dist.
#' @param rsids rsID(s) of the instrument(s).
#' @param epsilon De novo mutation rate (default is 1e-8)
#'
#' @return Matrix of propensity scores
#'
#' @export

# --- Main function ---
prop_score <- function(PHap, CHap, map, region, epsilon=1e-8) {
  # Sanity check
  if (!identical(colnames(CHap), map$rsid)) stop("SNP mismatch between map and input data...")

  # Genetic distances
  dist <- map$dist[-1]

  # Forward-backward weights
  Weights <- forward_backward_weights(PHap, CHap, dist, epsilon)

  # Number of instruments
  num_snps <- length(region$snps)

  # Positions of lower and upper instruments
  snps_col <- sapply(X=region$snps, FUN=function(x) which(colnames(CHap) == x))
  snps_pos <- sapply(X=region$snps, FUN=function(x) map$pos[map$rsid == x])

  # Positions of lower and upper bounds
  l <- which(colnames(CHap) == region$lower)
  h <- which(colnames(CHap) == region$upper)

  # Data frame containing name, base position and column position of SNPs
  dat_snps <- data.frame(name=region$snps, pos=snps_pos, col=snps_col)

  # Each combination of m and f to describe possible inheritance patterns
  vals_meiosis <- sample_space(c("m","f"), num_snps)

  # Probability of each inheritance pattern
  probmf <- apply(X=vals_meiosis, MARGIN=1, FUN=function(x) {
    x.length <- length(x)

    prob <- full_prob(PHap=PHap, CHap=CHap, Weights=Weights, dist=dist, col=dat_snps$col[1], l=l, h=h)

    if(x[1] == "f") prob <- 1-prob

    if(length(x) > 1) {
      for(k in 2:x.length) {
        partial <- partial_prob(PHap=PHap, CHap=CHap, Weights=Weights, dist=dist, u=x[k-1], dat=dat_snps, col=dat_snps$col[k], h=h)

        if(x[k] == "f") partial <- 1-partial

        prob <- prob*partial
      }
    }

    return(prob)
  })

  # Each combination of 0 and 1 to describe possible alleles
  vals_alleles <- sample_space(c(0,1), num_snps)

  # Probability of each allele (this is the propensity score)
  score <- apply(X=vals_alleles, MARGIN=1, FUN=function(alleles) {

    # Calculate emission probability
    emit <- apply(X=vals_meiosis, MARGIN=1, FUN=function(meiosis) {

      prob <- sapply(X=1:length(meiosis), FUN=function(x) {
        bool <- as.numeric(meiosis[x]=="m")
        P <- bool*PHap$m + (1-bool)*PHap$f

        return(ifelse(P[,dat_snps$col[x]] == alleles[x], 1-epsilon, epsilon))
      })

      prob <- apply(X=prob, MARGIN=1, FUN=function(x) prod(x))
      return(prob)
    })

    out <- sapply(X=1:nrow(vals_meiosis), FUN=function(x) probmf[,x]*emit[,x])
    return(rowSums(out))
  })

  return(score)
}

# --- Helper functions ---
full_prob <- function(PHap, CHap, Weights, dist, col, l, h, epsilon=1e-8) {
  # Column position of given SNP
  j <- col

  # Probabilities for the meiosis indicators
  r_probl <- 0.5*(1+exp(-2*sum(dist[l:(j-1)])))
  if(h == j+1) { r_probh = 1 } else { r_probh <- 0.5*(1+exp(-2*sum(dist[j:(h-2)]))) }

  # Backward weights
  beta <- list(m = Weights$backward$m[,h-1], f = Weights$backward$f[,h-1])

  # Forward weights
  alpha <- list(m = Weights$forward$m[,l+1], f = Weights$forward$f[,l+1])

  # Meiosis indicator distribution for parental haplotype
  propm <- (r_probh*beta$m + (1-r_probh)*beta$f)*(r_probl*alpha$m + (1-r_probl)*alpha$f)
  propf <- ((1-r_probh)*beta$m + r_probh*beta$f)*((1-r_probl)*alpha$m + r_probl*alpha$f)
  PrPm <- propm/(propm+propf)

  return(PrPm)
}


partial_prob <- function(PHap, CHap, Weights, dist, u="m", dat, col, h, epsilon=1e-8) {
  # Column position of given SNP
  j <- col

  # Column position of previous SNP
  j.prev <- dat$col[which(dat$col==j)-1]

  # Probabilities for the meiosis indicators
  r_probl <- 0.5*(1+exp(-2*sum(dist[j.prev:(j-1)])))
  if(h == j+1) { r_probh = 1 } else { r_probh <- 0.5*(1+exp(-2*sum(dist[j:(h-2)]))) }

  # Backward weights
  beta <- list(m = Weights$backward$m[,h-1], f = Weights$backward$f[,h-1])

  # Meiosis indicator distribution for parental haplotype
  bool <- as.integer(u=="m")
  propm <- (r_probh*beta$m + (1-r_probh)*beta$f)*(r_probl*bool + (1-r_probl)*(1-bool))
  propf <- ((1-r_probh)*beta$m + r_probh*beta$f)*((1-r_probl)*bool+ r_probl*(1-bool))
  PrPm <- propm/(propm+propf)

  return(PrPm)
}

sample_space <- function(vals, num) {
  list.vals <- rep(list(vals), num)
  combos <- as.matrix(expand.grid(list.vals))
  return(combos)
}

