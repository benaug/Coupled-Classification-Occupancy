NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  for(i in 1:n.species){
    Beta0.lam[i] ~ dnorm(1,sd=5)
    Beta0.theta[i] ~ dlogis(0,1)
    Beta0.psi[i] ~ dlogis(0,1)
    for(i2 in 1:n.species){
      log.alpha[i,i2] ~ dnorm(0,2) # priors for log(alpha) parameters
      G.alpha[i,i2] <- exp(log.alpha[i,i2])
    }
  }
  
  #--------------------------------------------------------------
  # occupancy model
  #--------------------------------------------------------------
  for(j in 1:J){
    for(i in 1:n.species){
      logit(psi[i,j]) <- Beta0.psi[i]
      z[i,j] ~ dbern(psi[i,j])
      for(k in 1:K){
        logit(theta[i,j,k]) <- Beta0.theta[i]
        w[i,j,k] ~ dbern(theta[i,j,k]*z[i,j]*K2D.indicator[j,k])
        log(lambda[i,j,k]) <- Beta0.lam[i]
      }
    }
    for(k in 1:K){
      bigLam[j,k] <- sum(lambda[1:n.species,j,k]*w[1:n.species,j,k])
      y2D[j,k] ~ dpois(bigLam[j,k]*K2D[j,k]) #likelihood for counts summed over species
      site.prob[1:n.species,j,k] <- (lambda[1:n.species,j,k]*w[1:n.species,j,k])/max(bigLam[j,k],2.220446e-16)
      #max part to avoid dividing by zero due to numerical underflow
      #and warning from Nimble for sites with all w=0 and no samples.
    }
  }
  for(i in 1:n.species){
    z.counts[i] <- sum(z[i,1:J]) #total number of occupied sites per species
    w.counts[i]<-sum(w[i,1:J,1:K])
  }
  #--------------------------------------------------------------
  # classification model
  #--------------------------------------------------------------
  for(l in 1:n.samples){
    IDtrue[l] ~ dcat(site.prob[1:n.species,G.site[l],G.occ[l]]) #"ecological prior" for species ID
    G.obs[l,1:n.species] ~ ddirch(G.alpha[IDtrue[l],1:n.species]) #likelihood of species ID
  }
})# end model
