# Forward and backward weights
fbweights <- function(Pm, Pf, Z, p, d, epsilon) {
  # Function for the probability of inheriting a SNP given meiosis indicator
  prsnp <- function(P,Q) {
    return(ifelse(P == Q, 1-epsilon, epsilon))
  }
  
  # Check equality of rows and columns
  if(!((ncol(Pm) == ncol(Pf)) & (ncol(Pm) == ncol(Z)))) stop('Columns are not equal')
  if(!((nrow(Pm) == nrow(Pf)) & (nrow(Pm) == nrow(Z)))) stop('Rows are not equal')
  
  # Sample size
  N <- nrow(Pm)
  
  # Construct forward weights
  Fweight <- namer(matrix(nrow=N, ncol=2*(p+1)), c('wm','wf'), 0:p)
  
  # Zero forward weight
  Fweight[,c('wm_0_','wf_0_')] <- c(0.5,0.5) 
  
  # Define forward weights recursively
  for (j in 1:p) {
    dat <- data.frame(Pmj = Pm[,grepl(paste('_',j,'_', sep=''), colnames(Pm))],
                      Pfj = Pf[,grepl(paste('_',j,'_', sep=''), colnames(Pf))],
                      Zj = Z[,grepl(paste('_',j,'_', sep=''), colnames(Z))])
    
    # Initial forward weight
    if (j == 1) {
      Fweight[,paste('wm_',j,'_',sep='')] <- 0.5*prsnp(dat$Pmj,dat$Zj)
      Fweight[,paste('wf_',j,'_',sep='')] <- 0.5*prsnp(dat$Pfj,dat$Zj)
    } 
    
    # Define remaining forward weights recursively
    if (j > 1) {
      # Forward weights one place behind
      fback <- data.frame(alpham = Fweight[,paste('wm_',j-1,'_',sep='')],
                          alphaf = Fweight[,paste('wf_',j-1,'_',sep='')])
                         
      # Recombination probability
      r_prob <- 0.5*(1+exp(-2*d[j-1])) 
      
      # Additional forward weights
      Fweight[,paste('wm_',j,'_',sep='')] <- prsnp(dat$Pmj,dat$Zj)*(r_prob*fback$alpham + (1-r_prob)*fback$alphaf)
      Fweight[,paste('wf_',j,'_',sep='')] <- prsnp(dat$Pfj,dat$Zj)*((1-r_prob)*fback$alpham + r_prob*fback$alphaf)
    }
  }
  
  # Construct backward weights
  Bweight <- namer(matrix(nrow=N, ncol=2*p), c('wm','wf'), 1:p)
  
  # Define backward weights recursively
  for(j in p:1) {
    # Initial backward weight
    if(j == p) {
      Bweight[,paste('wm_',j,'_',sep='')] <- 1
      Bweight[,paste('wf_',j,'_',sep='')] <- 1
    } 
    
    # Define remaining backward weights recursively
    if(j < p) {
      dat <- data.frame(Pmj = Pm[,grepl(paste('_',j+1,'_', sep=''), colnames(Pm))],
                        Pfj = Pf[,grepl(paste('_',j+1,'_', sep=''), colnames(Pf))],
                        Zj = Z[,grepl(paste('_',j+1,'_', sep=''), colnames(Z))])
      
      # Backward weights one place ahead
      bfor <- data.frame(betam = Bweight[,paste('wm_',j+1,'_',sep='')],
                         betaf = Bweight[,paste('wf_',j+1,'_',sep='')])
      
      # Recombination probability
      r_prob <- 0.5*(1+exp(-2*d[j])) 
      
      # Additional backward weights
      Bweight[,paste('wm_',j,'_',sep='')] <- bfor$betam*r_prob*prsnp(dat$Pmj,dat$Zj) + bfor$betaf*(1-r_prob)*prsnp(dat$Pfj,dat$Zj)
      Bweight[,paste('wf_',j,'_',sep='')] <- bfor$betam*(1-r_prob)*prsnp(dat$Pmj,dat$Zj) + bfor$betaf*r_prob*prsnp(dat$Pfj,dat$Zj)
    }
  }
  return(list(Bw=Bweight, Fw=Fweight))
}
