# Offspring haplotype (Markov chain model)
sampling_probability <- function(PHap, W, Jset, d, epsilon) {
  # Indices for lower bound, marker and upper bound respectively
  l <- Jset[1]; j <- Jset[2]; h <- Jset[3]
  
  # Probabilities for the meiosis indicators
  r_probl <- 0.5*(1+exp(-2*sum(d[l:(j-1)])))
  if(h == j+1) { r_probh = 1 } else { r_probh <- 0.5*(1+exp(-2*sum(d[j:(h-1)]))) }
  
  # Backward weights
  beta <- data.frame(mm = W$Bm[,paste('wm_',h-1,'_',sep='')],
                     mf = W$Bm[,paste('wf_',h-1,'_',sep='')],
                     fm = W$Bf[,paste('wm_',h-1,'_',sep='')],
                     ff = W$Bf[,paste('wf_',h-1,'_',sep='')])
  
  # Forward weights
  alpha <- data.frame(mm = W$Fm[,paste('wm_',l,'_',sep='')],
                      mf = W$Fm[,paste('wf_',l,'_',sep='')],
                      fm = W$Ff[,paste('wm_',l,'_',sep='')],
                      ff = W$Ff[,paste('wf_',l,'_',sep='')])
  
  # Meiosis indicator distribution for maternal haplotype
  propm <- (r_probh*beta$mm + (1-r_probh)*beta$mf)*(r_probl*alpha$mm + (1-r_probl)*alpha$mf)
  propf <- ((1-r_probh)*beta$mm + r_probh*beta$mf)*((1-r_probl)*alpha$mm + r_probl*alpha$mf)
  PrMm <- propm/(propm+propf)
  
  # Meiosis indicator distribution for paternal haplotype
  propm <- (r_probh*beta$fm + (1-r_probh)*beta$ff)*(r_probl*alpha$fm + (1-r_probl)*alpha$ff)
  propf <- ((1-r_probh)*beta$fm + r_probh*beta$ff)*((1-r_probl)*alpha$fm + r_probl*alpha$ff)
  PrFm <- propm/(propm+propf)
  
  # Probabilities for the offspring haplotypes
  PrM1 <- ifelse(PHap[,paste('Mm_',j,'_',sep='')] == 1, 1-epsilon, epsilon)*PrMm + ifelse(PHap[,paste('Mf_',j,'_',sep='')] == 1, 1-epsilon, epsilon)*(1-PrMm)
  PrF1 <- ifelse(PHap[,paste('Fm_',j,'_',sep='')] == 1, 1-epsilon, epsilon)*PrFm + ifelse(PHap[,paste('Ff_',j,'_',sep='')] == 1, 1-epsilon, epsilon)*(1-PrFm)
  
  return(list(PrM1=PrM1,PrF1=PrF1))
}
