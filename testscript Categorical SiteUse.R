#This test script expands the categorical model to allow further zero inflation
#through a site use indicator. Warning, updating latent states in this version
#is slower.

library(nimble)
library(coda)

source("sim.CCoccu.Categorical.SiteUse.R")
source("buildNimData.R")
source("NimbleModel Categorical SiteUse.R")
source("Latent Sampler Categorical SiteUse.R")

n.species <- 3 #species
J <- 50 #sites
K <- 10 #occasions
K2D <- matrix(1,nrow=J,ncol=K) #optional trap operation/effort

psi <- sample(c(0.2,0.3,0.4),n.species,replace=TRUE) #occupancy probs
theta <- sample(c(0.5,0.6,0.7),n.species,replace=TRUE) #site use probs
lambda <- sample(c(5,10,15),n.species,replace=TRUE) #detection rates|occupancy
 
#feature score parameters - n.species x n.levels. Data simulator considers n.levels=n.species
G.theta <- matrix(c(0.8,0.1,0.1,0.05,0.9,0.05,0.1,0.1,0.8),nrow=n.species,byrow=TRUE)

pObs <- 1 #We might not observe the partial ID covariates for all samples. Only considering all observed in MCMC.
pKnown <- 0.25 #We might know the true IDs for some samples, e.g. validation samples.

#simulate data
data <- sim.CCoccu.Categorical.siteUse(n.species=n.species,psi=psi,theta=theta,lambda=lambda,K=K,J=J,
                                     G.theta=G.theta,pObs=pObs,pKnown=pKnown,K2D=K2D)

#format data for nimble
nimbuild <- buildNimData(data)
n.samples <- length(data$G.obs)
n.levels <- ncol(G.theta)

#supply ID inits and data for validated samples
IDtrue.init <- IDtrue.data <- rep(NA,n.samples)
known.idx <- which(data$IDknown==1) #validated samples
IDtrue.init[known.idx] <- data$IDtrue[known.idx]
IDtrue.data[known.idx] <- data$IDtrue[known.idx]
unknown.idx <- which(data$IDknown==0) #unvalidated samples

#initializing ID to observed category. Assuming categories are species ID, but could be something else
IDtrue.init[unknown.idx] <- data$G.obs[unknown.idx]

#fit model
constants <- list(n.species=n.species,J=J,K=K,K2D=data$K2D,K2D.indicator=1*(K2D>0),
                G.site=data$G.site,G.occ=data$G.occ,n.samples=n.samples,n.levels=n.levels)
Niminits <- list(z=nimbuild$z.init,w=nimbuild$w.init,IDtrue=IDtrue.init)
Nimdata <- list(y2D=data$y2D,G.obs=data$G.obs,IDtrue=IDtrue.data)

#set parameters to monitor
parameters <- c('Beta0.lam','Beta0.psi','Beta0.theta','G.theta','z.counts','w.counts')
parameters2 <- c('IDtrue')
#thinning rates
nt <- 1
nt2 <- 25 #record fewer iters

# Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel,constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,monitors2=parameters2,thin=nt,thin2=nt2,useConjugacy=FALSE)

#remove z, w, and IDtrue samplers and replace
conf$removeSampler("z")
conf$removeSampler("w")
conf$removeSampler("IDtrue")
conf$addSampler(target = paste("z[1:",n.species,",1:",J,"]", sep=""),
                type = 'LatentSampler',control = list(G.site=data$G.site,G.occ=data$G.occ,n.samples=n.samples,n.species=n.species,
                                                  J=J,K=K,zpossible=nimbuild$zpossible,n.community=nimbuild$n.community,n.levels=n.levels,
                                                  plausible.j=nimbuild$plausible.j,plausible.jk=nimbuild$plausible.jk,
                                                  maxdet=nimbuild$maxdet,G.obs3D=nimbuild$G.obs3D,IDknown3D=nimbuild$IDknown3D,
                                                  IDtrue3D=nimbuild$IDtrue3D,K2D=K2D,IDknown=data$IDknown,y.known=nimbuild$y.known),
                silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for easier debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

burnin1 <- 250
mvSamples  <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:burnin1),]))

#Truth
rowSums(data$z) #number of occupied sites per species
rowSums(data$w) #number of used site-occasions per species
log(lambda) #beta0.lam
qlogis(psi) #beta0.psi
qlogis(theta) #beta0.theta
c(G.theta) #classification probability matrix parameters

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
