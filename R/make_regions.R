#' @name make_regions
#' @title Makes conditioning regions around SNPs
#' @description This function takes a genetic map, list of instruments (rsIDs) and
#' distance (in megabases) and returns a list of conditioning regions consisting of:
#' 1) the instruments included and 2) the lower and upper boundaries in base positions.
#'
#' @param map A data frame containing a genetic map.
#' @param snps A vector of rsIDs for the instruments.
#' @param mb Distance in megabases from each instrument to start conditioning on chromosome.
#'
#' @return A list of conditioning regions.
#' @export

make_regions <- function(map, snps, mb) {
  # Check that all supplied SNPs are found in the map file
  bool <- sapply(X=snps, FUN=function(x) x %in% map$rsid)
  if (!all(bool)) {
    bad_snps <- paste(snps[!bool], collapse=", ")
    stop(paste(bad_snps, "not found in map file"))
  }

  # Create a data frame containing rsids, base position and lower and upper positions
  # for regions around each SNP.
  dat <- data.frame(snps=snps)

  dat$pos <- sapply(X=snps, FUN=function(x) map$pos[map$rsid==x])

  dat$lower <- dat$pos-mb*1e6

  dat$upper <- dat$pos+mb*1e6

  # Cluster SNPs into regions, allowing for overlap.
  dat$region <- 1
  if(nrow(dat) > 1) {
    for(i in 2:nrow(dat)) {
      if ((dat$pos[i] >= dat$lower[i-1]) & (dat$pos[i] <= dat$upper[i-1])) dat$region[i] <- dat$region[i-1]
      else dat$region[i] <- dat$region[i-1]+1
    }
  }

  # Create a list summarising each region.
  regions <- lapply(X = 1:max(dat$region), FUN=function(x) {
    bool <- dat$region==x
    out <- list(
      snps = dat$snps[bool],
      pos = dat$pos[bool],
      lower = map$rsid[map$pos == max(map$pos[map$pos <= min(dat$lower[bool])])],
      upper = map$rsid[map$pos == min(map$pos[map$pos >= max(dat$lower[bool])])])
    return(out)
  })

  return(regions)
}

