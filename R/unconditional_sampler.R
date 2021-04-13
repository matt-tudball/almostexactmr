# Offspring haplotype (Markov chain model)
unconditional_sampler <- function(MHap, FHap, p, d, epsilon) {

  # Check equality of rows and columns
  if(!all(sapply(list(MHap$m, MHap$f, FHap$m, FHap$f), function(x) nrow(x) == nrow(MHap$m))))
    stop('Rows are not equal')
  if(!all(sapply(list(MHap$m, MHap$f, FHap$m, FHap$f), function(x) ncol(x) == ncol(MHap$m))))
    stop('Columns are not equal')

  # Sample size
  N <- nrow(MHap$m)

  # Initialise the hidden state
  Umsite <- matrix(nrow=N, ncol=p)
  Ufsite <- matrix(nrow=N, ncol=p)

  # Initial position: a = 0, b = 1
  Umsite[,1] <- rbinom(N,1,0.5)
  Ufsite[,1] <- rbinom(N,1,0.5)

  # Position inherited from parental haplotype based on Morgan distance
  for(j in 2:p) {
    which_site_m <- rbinom(N, 1, 0.5*(1+exp(-2*d[j-1])))
    which_site_f <- rbinom(N, 1, 0.5*(1+exp(-2*d[j-1])))

    prev_m <- Umsite[,j-1]; prev_f <- Ufsite[,j-1]
    Umsite[,j] <- which_site_m*prev_m + (1-which_site_m)*(1-prev_m)
    Ufsite[,j] <- which_site_f*prev_f + (1-which_site_f)*(1-prev_f)
  }

  # Generate the offspring haplotype
  Zm <- matrix(nrow=N, ncol=p)
  Zf <- matrix(nrow=N, ncol=p)

  for (j in 1:p) {
    # Inheritance from mother/father haplotype
    Mm <- MHap$m[,j]
    Mf <- MHap$f[,j]
    Fm <- FHap$m[,j]
    Ff <- FHap$f[,j]

    um_site <- Umsite[,j]; uf_site <- Ufsite[,j]
    Zm[,j] <- um_site*Mf + (1-um_site)*Mm
    Zf[,j] <- uf_site*Ff + (1-uf_site)*Fm

    # De novo mutation
    mutate_m <- rbinom(N,1,epsilon)
    mutate_f <- rbinom(N,1,epsilon)

    um <- Zm[,j]; uf <- Zf[,j]
    Zm[,j] <- mutate_m*(1-um) + (1-mutate_m)*um
    Zf[,j] <- mutate_f*(1-uf) + (1-mutate_f)*uf
  }

  return(list(m=Zm, f=Zf))
}
