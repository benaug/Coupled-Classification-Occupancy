## sampler to update z, w, and ID
LatentSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    G.site <- control$G.site
    G.occ <- control$G.occ
    n.samples <- control$n.samples
    n.species <- control$n.species
    n.levels <- control$n.levels
    J <- control$J
    K <- control$K
    zpossible <- control$zpossible
    n.community <- control$n.community
    plausible.j <- control$plausible.j
    plausible.jk <- control$plausible.jk
    maxdet <- control$maxdet
    G.obs3D <- control$G.obs3D
    IDknown3D <- control$IDknown3D
    IDtrue3D <- control$IDtrue3D
    IDknown <- control$IDknown
    y.known<- control$y.known
    K2D <- control$K2D
    calcNodes <- model$getDependencies(c("z","w","IDtrue")) #Everything we need to recalculate loglik for after update
  },
  run = function() {
    wpossible=zpossible
    psi <- model$psi
    theta <- model$theta
    lambda <- model$lambda
    G.theta <- model$G.theta
    y2D <- model$y2D
    sp.pred1D <- rep(0,sum(y2D))
    z.pred <- matrix(0,nrow=n.species,ncol=J)
    w.pred <- array(0,dim=c(n.species,J,K))
    sp.pred <- array(0,dim=c(J,K,maxdet))
    sp.idx <- 1
    
    #do one site at a time
    for(j in 1:J){
      lp.comm.z <- rep(log(0),n.community)
      lp.comm.w <- array(log(0),dim=c(n.community,K,n.community))
      #Not storing lp.spec because it can max out memory. Recalculating part necessary for 
      #sampling from joint distribution later.
      # lp.spec <- array(log(0),dim=c(n.community,K,n.community,maxdet,n.species))
      #get log probs for each z combo
      z.prob <- matrix(0,nrow=n.community,ncol=n.species)
      for(i in 1:n.community){
        z.prob[i,] <- dbinom(zpossible[i,],1,psi[,j],log=TRUE)
      }
      lp.z <- rep(log(0),n.community)
      logprobzz <- rep(0,n.community)
      for(zz in 1:n.community){ #for each possible z vector state at this site
        #ll for each plausible z vector
        if(plausible.j[j,zz]){
          lp.z[zz] <- sum(z.prob[zz,])
        }else{
          lp.z[zz] <- log(0)
        }
        #now do ll for each plausible w vector and then detections
        if(zz>1 & plausible.j[j,zz]){ #if detections at site j and plausible.j
          ll.sum.k <- rep(0,K)
          for(k in 1:K){ #for each occasion
            if(K2D[j,k]>0){ #if site-occasion operational
              #get log probs for each w combo
              w.prob <- matrix(0,nrow=n.community,ncol=n.species)
              for(i in 1:n.community){
                w.prob[i,] <- dbinom(wpossible[i,],1,theta[,j,k],log=TRUE)
              }
              
              logprobww <- rep(0,n.community)
              lp.w <- rep(log(0),n.community)
              for(ww in 1:n.community){ #for each possible w vector state at this site
                #need to skip this ww if zpossible rules out wpossible
                #cannot have w state = 1 if z state = 0
                zwconflicts <- sum(wpossible[ww,]==1&zpossible[zz,]==0)
                if(zwconflicts==0){
                  #ll for each plausible w vector
                  if(plausible.jk[j,k,ww]){
                    lp.w[ww] <- sum(w.prob[ww,]*zpossible[zz,]) #don't count w=0 if z=0
                  }else{
                    lp.w[ww] <- log(0)
                  }
                  if(ww>1 & plausible.jk[j,k,ww]){ #if detections on occasion k at site j and plausible.jk
                    pois.sum.terms <- rep(0,n.species)
                    for(i in 1:n.species){ #for each species
                      if(wpossible[ww,i]==1){ #site j used by species i on occasion k
                        pois.sum.terms[i] <- lambda[i,j,k]
                      }
                    }
                    logprobww[ww] <- dpois(y2D[j,k],sum(pois.sum.terms)*K2D[j,k],log=TRUE) #poisson encounter ll
                    # ecological prior part of encounter (add feature score part later) assuming known IDs, relax later
                    if(y2D[j,k]>0){ #if samples here
                      logmixprobs <- log(pois.sum.terms/sum(pois.sum.terms))
                      for(l in 1:y2D[j,k]){
                        if(IDknown3D[j,k,l]){
                          logprobww[ww] <- logprobww[ww]+logmixprobs[IDtrue3D[j,k,l]]+ #add ecological prior categorical likelihood
                            dcat(G.obs3D[j,k,l],G.theta[IDtrue3D[j,k,l],1:n.levels],log=TRUE)  #add feature score likelihood
                        }else{
                          logprob.tmp <- rep(0,n.species) #store feature score logprob over all possible species
                          for(i in 1:n.species){
                            logprob.tmp[i] <- dcat(G.obs3D[j,k,l],G.theta[i,1:n.levels],log=TRUE)
                          }
                          #log-sum-exp trick to avoid numerical overflow
                          lp.spec <- logprob.tmp+logmixprobs
                          maxlp <- max(lp.spec)#deal with overflow
                          logprobww[ww] <- logprobww[ww]+maxlp+log(sum(exp(lp.spec-maxlp))) #marginal ID likelihood
                        }
                      }
                    }
                  }
                }
              }
              #log-sum-exp trick to avoid numerical overflow
              lp.comm.w[zz,k,] <- lp.w+logprobww
              maxlp <- max(lp.comm.w[zz,k,])
              ll.sum.k[k] <- maxlp+log(sum(exp(lp.comm.w[zz,k,]-maxlp)))
            }
          }
          logprobzz[zz] <- sum(ll.sum.k)
        }
      }
      lp.comm.z <- lp.z+logprobzz
      
      #sample from joint distributions
      maxlp <- max(lp.comm.z)#deal with overflow
      z.probs <- exp(lp.comm.z-maxlp)
      z.probs <- z.probs/sum(z.probs)
      pick.z <- rcat(1,z.probs)
      z.pred[,j] <- zpossible[pick.z,]
      if(pick.z>1){
        for(k in 1:K){
          if(K2D[j,k]>0){
            maxlp <- max(lp.comm.w[pick.z,k,])#deal with overflow
            w.probs <- exp(lp.comm.w[pick.z,k,]-maxlp)
            w.probs <- w.probs/sum(w.probs)
            pick.w <- rcat(1,w.probs)
            w.pred[,j,k] <- wpossible[pick.w,]
            if(pick.w>1){
              if(y2D[j,k]>0){
                for(l in 1:y2D[j,k]){
                  if(IDknown3D[j,k,l]==0){
                    #calculate lp.sp here instead of storing it aboveso we don't need
                    #to store an array that can max out memory with more species
                    pois.sum.terms <- rep(0,n.species)
                    for(i in 1:n.species){ #for each species
                      if(w.pred[i,j,k]==1){ #site j used by species i on occasion k
                        pois.sum.terms[i] <- lambda[i,j,k]
                      }
                    }
                    logmixprobs <- log(pois.sum.terms/sum(pois.sum.terms))
                    logprob.tmp <- rep(0,n.species) #store feature score logprob over all possible species
                    for(i in 1:n.species){
                      logprob.tmp[i] <- dcat(G.obs3D[j,k,l],G.theta[i,1:n.levels],log=TRUE)
                    }
                    lp.spec <- logprob.tmp+logmixprobs
                    maxlp <- max(lp.spec)#deal with overflow
                    l.probs <- exp(lp.spec-maxlp)
                    l.probs <- l.probs/sum(l.probs)
                    sp.pred[j,k,l] <- pick.l <- rcat(1,l.probs)
                  }else{
                    sp.pred[j,k,l] <- IDtrue3D[j,k,l]
                  }
                  sp.pred1D[sp.idx] <- sp.pred[j,k,l]
                  sp.idx <- sp.idx+1
                }
              }
            }
          }
        }
      }
    }
    model$z <<- z.pred
    model$w <<- w.pred
    model$IDtrue <<- sp.pred1D
    
    # #update lp
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)