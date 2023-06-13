##############known ID - w (site use) and z (occupancy) marginalized out ##################
NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  for(i in 1:n.species){
    Beta0.lam[i] ~ dnorm(1,sd=5)
    Beta0.psi[i] ~ dlogis(0,1)
    for(i2 in 1:n.levels){
      alpha[i,i2] <- 1 #classification prob prior parameters
      #all 1's here means prior is each classification equally likely. Can relax
      #by using a prior for each alpha.
    }
    G.theta[i,1:n.levels] ~ ddirch(alpha[i,1:n.species])
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
    for(k in 1:K){ #likelihood for counts summed over species
      bigLam[j,k] <- sum(lambda[1:n.species,j,k]*z[1:n.species,j])
      y2D[j,k] ~ dpois(bigLam[j,k]*K2D[j,k])
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
    G.obs[l] ~ dcat(G.theta[IDtrue[l],1:n.levels]) #likelihood of species ID
  }
})# end model

