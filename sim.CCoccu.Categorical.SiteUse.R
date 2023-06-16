sim.CCoccu.Categorical.SiteUse <-
  function(n.species=NA,psi=NA,theta=NA,lambda=NA,lambda.sd=NA,theta.sd=NA,
           K=NA,K2D=NA,J=NA,n.cov=NA,G.theta=NA,
           pObs=1,pKnown=0){
    if(nrow(G.theta)!=n.species)stop("G.theta must be of length n.species")
    if(!missing(K2D)){
      print("Using K2D to simulate instead of K")
      if(dim(K2D)[1]!=J||dim(K2D)[2]!=K)stop("K2D must be of dimension J x K")
    }else{
      K2D <- array(1,dim=c(J,K))
    }
    #get lambdas, either fixed or random effects
    lambda.use <- array(NA,dim=c(n.species,J,K))
    if(!any(is.na(lambda.sd))){
      print("Simulating species by site by occasion encounter random effects on log scale.")
      log.lambda <- log(lambda)
      for(i in 1:n.species){
        lambda.use[i,,] <- matrix(exp(log.lambda[i]+rnorm(J*K,0,lambda.sd[i])),nrow=J,ncol=K)
      }
    }else{
      for(i in 1:n.species){
        lambda.use[i,,] <- lambda[i]
      }
    }
    #get theta, either fixed or random effects
    theta.use <- matrix(NA,n.species,J)
    if(!any(is.na(theta.sd))){
      print("Simulating species by site theta random effects on logit scale.")
      logit.theta <- qlogis(theta)
      for(i in 1:n.species){
        theta.use[i,] <- plogis(logit.theta[i]+rnorm(J,0,theta.sd[i]))
      }
    }else{
      for(i in 1:n.species){
        theta.use[i,] <- theta[i]
      }
    }
    
    #simulate ecological model data
    z <- matrix(0,nrow=n.species,ncol=J)
    w <- y.true <- array(0,dim=c(n.species,J,K))
    for(i in 1:n.species){
      for(j in 1:J){
        z[i,j] <- rbinom(1,1,psi[i]) #site use
        for(k in 1:K){
          if(K2D[j,k]>0){ #site-occasion operation
            w[i,j,k] <- rbinom(1,1,theta.use[i,j]*z[i,j]) #site visitation|site use
            if(w[i,j,k]==1){
              y.true[i,j,k] <- rpois(1,lambda.use[i,j,k]*K2D[j,k]) #detection|site visitation
            }
          }
        }
      }
    }
    y2D <- apply(y.true,c(2,3),sum)
    
    #disaggregate counts - order site -> occasion -> sample
    n.samples <- sum(y.true)
    ID <- G.site <- G.occ <- rep(NA,n.samples)
    idx <- 1
    for(j in 1:J){ #loop through sites
      for(k in 1:K){ #then occasions
        for(i in 1:n.species){ #then species
          if(y.true[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y.true[i,j,k]){ #then samples
              ID[idx] <- i
              G.site[idx] <- j
              G.occ[idx] <- k
              idx <- idx+1
            }
          }
        }
      }
    }
    
    #Simulate observed G with measurement error
    G.obs <- rep(NA,n.samples)
    for(i in 1:n.species){
      G.obs[ID==i]=rcat(sum(ID==i),G.theta[i,])
    }
    
    
    #Simulate missing G.obs values
    rem <- which(rbinom(n.samples,1,pObs)==0)
    G.obs[rem] <- NA
    
    #Simulate known ID samples
    IDknown <- rbinom(n.samples,1,pKnown)
    parm.vals <- list(psi=psi,theta=theta,
                      lambda=lambda,lambda.sd=lambda.sd,theta.sd=theta.sd,G.theta=G.theta,pObs=pObs,pKnown=pKnown)
    return(list(y2D=y2D,G.obs=G.obs,G.site=G.site,G.occ=G.occ,IDtrue=ID,IDknown=IDknown, #observed data
                y.true=y.true,z=z,w=w, #true data
                n.species=n.species,parm.vals=parm.vals, #sim params
                K2D=K2D))
  }
