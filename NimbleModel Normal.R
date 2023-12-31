NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  for(i in 1:n.species){
    Beta0.lam[i] ~ dnorm(1,sd=5)
    Beta0.psi[i] ~ dlogis(0,1)
    G.mu[i] ~ dnorm(0,sd=10) #using normal prior on mu for conjugate sampler
    G.tau[i] ~ dgamma(0.001,0.001) #using gamma prior on tau for conjugate sampler
    G.sigma[i] <- sqrt(1/G.tau[i])
  }
  #--------------------------------------------------------------
  # occupancy model
  #--------------------------------------------------------------
  for(j in 1:J){
    for(i in 1:n.species){
      logit(psi[i,j]) <- Beta0.psi[i]
      z[i,j] ~ dbern(psi[i,j])
      for(k in 1:K){
        log(lambda[i,j,k]) <- Beta0.lam[i]
      }
    }
    for(k in 1:K){
      bigLam[j,k] <- sum(lambda[1:n.species,j,k]*z[1:n.species,j])
      y2D[j,k] ~ dpois(bigLam[j,k]*K2D[j,k]) #likelihood for counts summed over species
      site.prob[1:n.species,j,k] <- (lambda[1:n.species,j,k]*z[1:n.species,j])/max(bigLam[j,k],2.220446e-16)
      #max part to avoid dividing by zero due to numerical underflow
    }
  }
  for(i in 1:n.species){
    z.counts[i] <- sum(z[i,1:J]) #total number of occupied sites per species
  }
  #--------------------------------------------------------------
  # classification model
  #--------------------------------------------------------------
  for(l in 1:n.samples){
    IDtrue[l] ~ dcat(site.prob[1:n.species,G.site[l],G.occ[l]]) #"ecological prior" for species ID
    G.obs[l] ~ dnorm(G.mu[IDtrue[l]],tau=G.tau[IDtrue[l]]) #likelihood of species ID
  }
})# end model

