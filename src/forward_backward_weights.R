# Forward and backward weights
fbweights <- function(PHap, Z, p, d, epsilon) {
  # Function for the probability of inheriting a SNP given meiosis indicator
  prsnp <- function(P,Q) {
    return(ifelse(P == Q, 1-epsilon, epsilon))
  }
  
  # Sample size
  N <- nrow(PHap)
  
  # Construct forward weights
  Fmweight <- namer(matrix(nrow=N, ncol=2*(p+1)), c('wm','wf'), 0:p)
  Ffweight <- namer(matrix(nrow=N, ncol=2*(p+1)), c('wm','wf'), 0:p)
  
  # Zero forward weight
  Fmweight[,c('wm_0_','wf_0_')] <- c(0.5,0.5) 
  Ffweight[,c('wm_0_','wf_0_')] <- c(0.5,0.5) 
  
  # Define forward weights recursively
  for (j in 1:p) {
    dat <- data.frame(Mm = PHap[,paste('Mm_',j,'_',sep='')],
                      Mf = PHap[,paste('Mf_',j,'_',sep='')],
                      Fm = PHap[,paste('Fm_',j,'_',sep='')],
                      Ff = PHap[,paste('Ff_',j,'_',sep='')],
                      Zm = Z$m[,paste('Zm_',j,'_',sep='')],
                      Zf = Z$f[,paste('Zf_',j,'_',sep='')])
    
    # Initial forward weight
    if (j == 1) {
      Fmweight[,paste('wm_',j,'_',sep='')] <- 0.5*prsnp(dat$Mm,dat$Zm)
      Fmweight[,paste('wf_',j,'_',sep='')] <- 0.5*prsnp(dat$Mf,dat$Zm)
      Ffweight[,paste('wm_',j,'_',sep='')] <- 0.5*prsnp(dat$Fm,dat$Zf)
      Ffweight[,paste('wf_',j,'_',sep='')] <- 0.5*prsnp(dat$Ff,dat$Zf)
    } 
    
    # Define remaining forward weights recursively
    if (j > 1) {
      # Forward weights one place behind
      fback <- data.frame(alphamm = Fmweight[,paste('wm_',j-1,'_',sep='')],
                          alphamf = Fmweight[,paste('wf_',j-1,'_',sep='')],
                          alphafm = Ffweight[,paste('wm_',j-1,'_',sep='')],
                          alphaff = Ffweight[,paste('wf_',j-1,'_',sep='')])
      
      # Recombination probability
      r_prob <- 0.5*(1+exp(-2*d[j-1])) 
      
      # Additional forward weights
      Fmweight[,paste('wm_',j,'_',sep='')] <- prsnp(dat$Mm,dat$Zm)*(r_prob*fback$alphamm + (1-r_prob)*fback$alphamf)
      Fmweight[,paste('wf_',j,'_',sep='')] <- prsnp(dat$Mf,dat$Zm)*((1-r_prob)*fback$alphamm + r_prob*fback$alphamf)
      Ffweight[,paste('wm_',j,'_',sep='')] <- prsnp(dat$Fm,dat$Zf)*(r_prob*fback$alphafm + (1-r_prob)*fback$alphaff)
      Ffweight[,paste('wf_',j,'_',sep='')] <- prsnp(dat$Ff,dat$Zf)*((1-r_prob)*fback$alphafm + r_prob*fback$alphaff)
    }
  }
  
  # Construct backward weights
  Bmweight <- namer(matrix(nrow=N, ncol=2*p), c('wm','wf'), 1:p)
  Bfweight <- namer(matrix(nrow=N, ncol=2*p), c('wm','wf'), 1:p)
  
  # Define backward weights recursively
  for(j in p:1) {
    # Initial backward weight
    if(j == p) {
      Bmweight[,paste('wm_',j,'_',sep='')] <- 1
      Bmweight[,paste('wf_',j,'_',sep='')] <- 1
      Bfweight[,paste('wm_',j,'_',sep='')] <- 1
      Bfweight[,paste('wf_',j,'_',sep='')] <- 1
    } 
    
    # Define remaining backward weights recursively
    if(j < p) {
      dat <- data.frame(Mm = PHap[,paste('Mm_',j+1,'_',sep='')],
                        Mf = PHap[,paste('Mf_',j+1,'_',sep='')],
                        Fm = PHap[,paste('Fm_',j+1,'_',sep='')],
                        Ff = PHap[,paste('Ff_',j+1,'_',sep='')],
                        Zm = Z$m[,paste('Zm_',j+1,'_',sep='')],
                        Zf = Z$f[,paste('Zf_',j+1,'_',sep='')])
      
      # Backward weights one place ahead
      bfor <- data.frame(betamm = Bmweight[,paste('wm_',j+1,'_',sep='')],
                         betamf = Bmweight[,paste('wf_',j+1,'_',sep='')],
                         betafm = Bfweight[,paste('wm_',j+1,'_',sep='')],
                         betaff = Bfweight[,paste('wf_',j+1,'_',sep='')])
      
      # Recombination probability
      r_prob <- 0.5*(1+exp(-2*d[j])) 
      
      # Additional backward weights
      Bmweight[,paste('wm_',j,'_',sep='')] <- bfor$betamm*r_prob*prsnp(dat$Mm,dat$Zm) + bfor$betamf*(1-r_prob)*prsnp(dat$Mf,dat$Zm)
      Bmweight[,paste('wf_',j,'_',sep='')] <- bfor$betamm*(1-r_prob)*prsnp(dat$Mm,dat$Zm) + bfor$betamf*r_prob*prsnp(dat$Mf,dat$Zm)
      Bfweight[,paste('wm_',j,'_',sep='')] <- bfor$betafm*r_prob*prsnp(dat$Fm,dat$Zf) + bfor$betaff*(1-r_prob)*prsnp(dat$Ff,dat$Zf)
      Bfweight[,paste('wf_',j,'_',sep='')] <- bfor$betafm*(1-r_prob)*prsnp(dat$Fm,dat$Zf) + bfor$betaff*r_prob*prsnp(dat$Ff,dat$Zf)
    }
  }
  return(list(Bm=Bmweight, Bf=Bfweight, Fm=Fmweight, Ff=Ffweight))
}
