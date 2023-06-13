#This test script uses a sampler for z and sample species ID that samples from their marginal distributions.
#You can not use the provided sampler and use nimble-assigned samplers, but this approach will often
#not converge/fully explore the posterior. While the provided sampler will converge/fully explore the posterior,
#it is increasingly slower as you add more species and can't handle more than 8 or so with decent runtimes.

#This model works with general categorical covariates aside from "species ID", but is set up here to just handle
#"species ID". You'll see n.levels below. This is set to n.species here because the categorical covaraites are species ID.

library(nimble)
library(coda)

source("sim.CCoccu.Categorical.R")
source("buildNimData.R")
source("NimbleModel Categorical.R")
source("Latent Sampler Categorical.R")

n.species <- 3 #species
J <- 50 #sites
K <- 10 #occasions
K2D <- matrix(1,nrow=J,ncol=K) #optional trap operation/effort

psi <- sample(c(0.2,0.3,0.4),n.species,replace=TRUE) #occupancy probs
lambda <- sample(c(5,10,15),n.species,replace=TRUE) #detection rates|occupancy
 
#feature score parameters - n.species x n.levels. Data simulator considers n.levels=n.species
G.theta <- matrix(c(0.8,0.1,0.1,0.05,0.9,0.05,0.1,0.1,0.8),nrow=n.species,byrow=TRUE)

pObs <- 1 #We might not observe the partial ID covariates for all samples. Only considering all observed in MCMC.
pKnown <- 0.25 #We might know the true IDs for some samples, e.g. validation samples.
#Data simulator assumes these are selected at random across sites, but can choose them in any way you want.

#simulate data
data <- sim.CCoccu.Categorical(n.species=n.species,psi=psi,lambda=lambda,K=K,J=J,
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

#initializing ID to observed species ID
IDtrue.init[unknown.idx] <- data$G.obs[unknown.idx]

#fit model
constants <- list(n.species=n.species,J=J,K=K,K2D=data$K2D,
                G.site=data$G.site,G.occ=data$G.occ,n.samples=n.samples,n.levels=n.levels)
Niminits <- list(z=nimbuild$z.init,IDtrue=IDtrue.init)
Nimdata <- list(y2D=data$y2D,G.obs=data$G.obs,IDtrue=IDtrue.data)

#set parameters to monitor
parameters <- c('Beta0.lam','Beta0.psi','G.theta','z.counts')
parameters2 <- c('IDtrue')
#thinning rates
nt <- 1
nt2 <- 25 #record fewer iters

# Build the model, configure the mcmc, and compileConfigure
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,monitors2=parameters2, thin=nt, thin2=nt2, useConjugacy=False)

#remove z, w, and IDtrue samplers and replace
#NOTE: If you skip this replacement, you can use the nimble-assigned samplers and see how poorly they perform
conf$removeSampler("z")
conf$removeSampler("IDtrue")
conf$addSampler(target = paste("z[1:",n.species,",1:",J,"]", sep=""),
                type = 'LatentSampler',
                control = list(G.site=data$G.site,G.occ=data$G.occ,n.samples=n.samples,n.species=n.species,
                               J=J,K=K,zpossible=nimbuild$zpossible,n.community=nimbuild$n.community,
                               n.levels=n.levels,plausible.j=nimbuild$plausible.j,maxdet=nimbuild$maxdet,
                               G.obs3D=nimbuild$G.obs3D,IDknown3D=nimbuild$IDknown3D,
                               IDtrue3D=nimbuild$IDtrue3D,K2D=K2D,IDknown=data$IDknown),
                silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for easier debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
Cmcmc$run(1000,reset=FALSE) #Can keep extending the run by rerunning this line
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

burnin1 <- 250
mvSamples  <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:burnin1),]))

#Truth
rowSums(data$z) #number of occupied sites per species
log(lambda) #beta0.lam
qlogis(psi) #beta0.psi
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

