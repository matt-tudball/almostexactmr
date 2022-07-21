#' @name make_clusters
#' @title Make clusters for parallel processing
#' @description This function takes a number of cores and creates the corresponding
#' number of clusters.
#'
#' @param cores Number of threads to be used for parallel processing.
#'
#' @return A cluster object.
#'
#' @importFrom parallel makeCluster clusterExport

make_clusters <- function(cores=1) {
  cluster <- parallel::makeCluster(cores, outfile = "")

  parallel::clusterExport(
    cl=cluster,
    varlist=c(ls(), as.vector(lsf.str(env=environment(run_test)))),
    envir = environment())

  return(cluster)
}
