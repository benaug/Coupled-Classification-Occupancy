buildNimData.Dirichlet <- function(data=data){
  n.species=data$n.species
  #Build zpossible and wpossible
  if(n.species==2){
    zpossible=rbind(c(1,1),c(1,0),c(0,1))
  }else{
    zpossible=cbind(c(1,1,0),c(1,0,1))
    if(n.species>2){
      for (i in (3:n.species)) {
        zpossible=cbind(c(rep(1,2^(i-1)),rep(0,((2^(i-1))-1))), rbind(zpossible, rep(0, (i-1)), zpossible))
      }
    }
  }
  zpossible=rbind(0,zpossible)
  wpossible=zpossible
  n.community=nrow(zpossible)

  #pull out observed data
  y2D=data$y2D
  G.obs=data$G.obs
  G.site=data$G.site
  G.occ=data$G.occ
  IDtrue=data$IDtrue
  IDknown=data$IDknown
  n.samples=nrow(G.obs)
  J=dim(y2D)[1]
  K=dim(y2D)[2]
    
  #restructure feature score structures for more efficient MCMC
  maxdet=max(y2D)
  IDtrue3D=IDknown3D=array(0,dim=c(J,K,maxdet))
  G.obs4D=array(0,dim=c(J,K,maxdet,n.species))
  for(j in 1:J){
    for(k in 1:K){
      idx=which(G.site==j&G.occ==k)
      if(length(idx)==0)next
      IDtrue3D[j,k,1:length(idx)]=IDtrue[idx]
      IDknown3D[j,k,1:length(idx)]=IDknown[idx]
      G.obs4D[j,k,1:length(idx),]=G.obs[idx,]
    }
  }
  
  #build y.true.known for validated samples
  y.true.known=array(0,dim=c(n.species,J,K))
  for(l in 1:n.samples){
    if(IDknown[l]){
      y.true.known[IDtrue[l],G.site[l],G.occ[l]]=y.true.known[IDtrue[l],G.site[l],G.occ[l]]+1
    }
  }
  y2D.known=apply(y.true.known,c(2,3),sum)
  
  #build structures to only marginalize over latent states consistent with observed data
  
  #which sites have detections?
  # any.detections.j=1*(rowSums(y2D.known)>0)
  any.detections.j=1*(rowSums(y2D)>0)
  plausible.j=matrix(1,J,n.community)
  plausible.j[which(any.detections.j==1),1]=0
  #more to remove based on known ID detections
  y2Dknownij=apply(y.true.known,c(1,2),sum)
  for(j in 1:J){
    if(sum(y2Dknownij[,j]>0)){
      for(i in 1:n.community){
        plausible.j[j,i]=1*(all(which(y2Dknownij[,j]>0)%in%which(zpossible[i,]==1)))
      }
    }
  }
  #which site occasions have detections?
  # any.detections.jk=1*(apply(y.true.known,c(2,3),sum)>0)
  any.detections.jk=1*(y2D>0)
  plausible.jk=array(1,dim=c(J,K,n.community))
  for(j in 1:J){
    plausible.jk[j,which(any.detections.jk[j,]==1),1]=0 #if detection at j k, all 0 is not plausible w state
  }
  #more to remove based on known ID detections
  #assuming all known for now
  for(j in 1:J){
    for(k in 1:K){
      if(sum(y.true.known[,j,k])>0){
        for(i in 1:n.community){
          plausible.jk[j,k,i]=1*(all(which(y.true.known[,j,k]>0)%in%which(wpossible[i,]==1)))
        }
      }
    }
  }
  
  z.init=matrix(1*(rowSums(y2D)>0),nrow=n.species,ncol=J,byrow=TRUE)
  w.init=array(0,dim=c(n.species,J,K))
  for(i in 1:n.species){
    w.init[i,,]=1*(y2D>0)
  }

  return(list(zpossible=zpossible,n.community=n.community,
              G.obs4D=G.obs4D,IDknown3D=IDknown3D,IDtrue3D=IDtrue3D,
              plausible.j=plausible.j,plausible.jk=plausible.jk,maxdet=maxdet,y.known=y.true.known,
              z.init=z.init,w.init=w.init))
}

