# Offspring haplotype (Markov chain model)
prop_score <- function(Pm, Pf, W, Jset, d, epsilon) {
  # Indices for lower bound, marker and upper bound respectively
  l <- Jset[1]; j <- Jset[2]; h <- Jset[3]
  
  # Probabilities for the meiosis indicators
  r_probl <- 0.5*(1+exp(-2*sum(d[l:(j-1)])))
  if(h == j+1) { r_probh = 1 } else { r_probh <- 0.5*(1+exp(-2*sum(d[j:(h-2)]))) }
  
  # Backward weights
  beta <- data.frame(m = W$Bw[,paste('wm_',h-1,'_',sep='')],
                     f = W$Bw[,paste('wf_',h-1,'_',sep='')])
  
  # Forward weights
  alpha <- data.frame(m = W$Fw[,paste('wm_',l,'_',sep='')],
                      f = W$Fw[,paste('wf_',l,'_',sep='')])
  
  # Meiosis indicator distribution for parental haplotype
  propm <- (r_probh*beta$m + (1-r_probh)*beta$f)*(r_probl*alpha$m + (1-r_probl)*alpha$f)
  propf <- ((1-r_probh)*beta$m + r_probh*beta$f)*((1-r_probl)*alpha$m + r_probl*alpha$f)
  PrPm <- propm/(propm+propf)
  
  # Probabilities for the offspring haplotypes
  PrP1 <- ifelse(Pm[,grepl(paste('_',j,'_', sep=''), colnames(Pm))] == 1, 1-epsilon, epsilon)*PrPm + 
          ifelse(Pf[,grepl(paste('_',j,'_', sep=''), colnames(Pf))] == 1, 1-epsilon, epsilon)*(1-PrPm)
  
  return(PrP1)
}
