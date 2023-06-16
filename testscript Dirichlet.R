#The model in this script is a regular occupancy model where instead of observing "species ID", you observe
#a Dirichlet random variable. To improve estimation, you may also have validation data with Dirichlet observations
#and a true species IDas might be available from humans. 
#This script assumes the validation data is randomly chosen from the focal
#survey, but independent validation data may also be used with some modification. An advantage of selecting
#validation data from the focal survey is that it fixes some z states (and w states for SiteUse versions).

#This test script uses a sampler for z and sample species ID that samples from their marginal distributions.
#You can not use the provided sampler and use nimble-assigned samplers, but this approach will sometimes
#not converge/fully explore the posterior. While the provided sampler will converge/fully explore the posterior,
#it is increasingly slower as you add more species and can't handle more than 8 or so with decent runtimes.

#This model is slower with Dirichlet RVs

library(nimble)
library(coda)

source("sim.CCoccu.Dirichlet.R")
source("buildNimData.Dirichlet.R")
source("NimbleModel Dirichlet.R")
source("Latent Sampler Dirichlet.R")

n.species <- 3 #species
J <- 50 #sites
K <- 10 #occasions
K2D <- matrix(1,nrow=J,ncol=K) #optional trap operation/effort

psi <- sample(c(0.2,0.3,0.4),n.species,replace=TRUE) #occupancy probs
lambda <- sample(c(5,10,15),n.species,replace=TRUE) #detection rates|occupancy
 
#feature score parameters - n.species x n.species
G.alpha=matrix(0,nrow=n.species,ncol=n.species)
for(i in 1:n.species){
  G.alpha[i,]=1
  G.alpha[i,i]=5 #higher values on diagonals means better discrimination
}

G.alpha

#visualize overlap in feature score distributions
#plot here looks decent for 3 species, probably won't be legible for more than 4 or so
par(mfrow=c(n.species,n.species),ask=FALSE)
par(mar = c(5.1, 4.1, 1.1, 1.1))
n.sims=10000
sims=array(NA,dim=c(n.species,n.species,n.sims))
for(i in 1:n.species){
  for(j in 1:n.sims){
    sims[i,,j]=nimble::rdirch(1,alpha=G.alpha[i,])
  }
  for(i2 in 1:n.species){
    hist(sims[i,i2,],freq=FALSE,main="",xlab=paste("Observed Value for Species",i2),ylab=paste("True Species",i))
  }
}
par(mar = c(5.1, 4.1, 4.1, 2.1))

pObs <- 1 #We might not observe the partial ID covariates for all samples. Only considering all observed in MCMC.
pKnown <- 0.10 #We might know the true IDs for some samples, e.g. validation samples. If not, set to 0.
#Data simulator assumes these are selected at random across sites, but can choose them in any way you want.

#simulate data
data <- sim.CCoccu.Dirichlet(n.species=n.species,psi=psi,lambda=lambda,K=K,J=J,
                                     G.alpha=G.alpha,pObs=pObs,pKnown=pKnown,K2D=K2D)

#what is the observed data?
#1) we know how many detections there were for each trap-occasion, just not which species they are
str(data$y2D)
#2) we might have known ID samples. Here, they are randomly selected from focal survey, but can be from elsewhere
head(data$IDtrue[data$IDknown==1],10)
#we can link these to a site of detection
head(data$G.site[data$IDknown==1],10)
#and the occasion of detection
head(data$G.occ[data$IDknown==1],10)
#3) then we have unknown ID samples, where we also know the site and occasion of capture
head(data$G.site[data$IDknown==0],10)
head(data$G.occ[data$IDknown==0],10)
#4) finally, we observe a Dirichlet random variable, which is a vector 
#we observe this for validated and unvalidated samples
head(data$G.obs,10)

#format data for nimble
nimbuild <- buildNimData.Dirichlet(data)
n.samples <- nrow(data$G.obs)
n.levels <- ncol(G.alpha)

#supply ID inits and data for validated samples
IDtrue.init <- IDtrue.data <- rep(NA,n.samples)
known.idx <- which(data$IDknown==1) #validated samples
IDtrue.init[known.idx] <- data$IDtrue[known.idx]
IDtrue.data[known.idx] <- data$IDtrue[known.idx]
unknown.idx <- which(data$IDknown==0) #unvalidated samples

#initializing ID to species with highest value - helps prevent label switching, speeds convergence
IDtrue.init[unknown.idx]=apply(data$G.obs[unknown.idx,],1,function(x){which(x==max(x))})

#fit model
constants <- list(n.species=n.species,J=J,K=K,K2D=data$K2D,
                G.site=data$G.site,G.occ=data$G.occ,n.samples=n.samples)
Niminits <- list(z=nimbuild$z.init,IDtrue=IDtrue.init)
Nimdata <- list(y2D=data$y2D,G.obs=data$G.obs,IDtrue=IDtrue.data)

#set parameters to monitor
parameters <- c('Beta0.lam','Beta0.psi','G.alpha','z.counts')
parameters2 <- c('IDtrue')
#thinning rates
nt <- 1
nt2 <- 25 #record fewer iters

# Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,monitors2=parameters2, thin=nt, thin2=nt2, useConjugacy=FALSE)

#remove z, w, and IDtrue samplers and replace
#NOTE: If you skip this replacement, you can use the nimble-assigned samplers and see how poorly they perform
conf$removeSampler("z")
conf$removeSampler("IDtrue")
conf$addSampler(target = paste("z[1:",n.species,",1:",J,"]", sep=""),
                type = 'LatentSampler',
                control = list(G.site=data$G.site,G.occ=data$G.occ,n.samples=n.samples,n.species=n.species,
                               J=J,K=K,zpossible=nimbuild$zpossible,n.community=nimbuild$n.community,
                               n.levels=n.levels,plausible.j=nimbuild$plausible.j,maxdet=nimbuild$maxdet,
                               G.obs4D=nimbuild$G.obs4D,IDknown3D=nimbuild$IDknown3D,
                               IDtrue3D=nimbuild$IDtrue3D,K2D=K2D,IDknown=data$IDknown),
                silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for easier debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

burnin1 <- 250
mvSamples  <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:burnin1),]))

#Note! If no validation data is used (or informative priors), you may see label switching
#The parameter estimates are correct, but it might be hard to tell which species is which
#using real data (as opposed to simulated data where you can easily determine this).

#Truth
rowSums(data$z) #number of occupied sites per species
log(lambda) #beta0.lam
qlogis(psi) #beta0.psi
c(G.alpha) #classification probability matrix parameters

#Species ID posteriors
burnin2 <- 5 #burnin for mvSamples 2 with different thinning rate
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
mvSamples2 <- mvSamples2[-c(1:burnin2),]

#Can check posterior ID for latent ID samples
idx <- which(data$IDknown==0) #list of latent ID samples
check <- 25 #which one to look at
table(mvSamples2[,idx[check]]) #number of iterations sample is assigned to species X
#posterior probability sample is of species X (need good effective sample size for reliable estimates)
table(mvSamples2[,idx[check]])/sum(table(mvSamples2[,idx[check]]))

