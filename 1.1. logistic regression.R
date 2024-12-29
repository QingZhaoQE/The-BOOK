#=================================================================================================
# Logistic regression
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 1_Background/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(1)

# Basic values
nsite <- 360 # number of sites
ncovs <- 3   # number of environmental covariates
beta <- c(0, 0.8, -0.4, 0.2) # intercept and slopes for occurrence

# Define function for data simulation
data_sim <- function(nsite, ncovs, beta) {
  x <- matrix(rnorm(nsite*ncovs, 0, 1), nsite, ncovs) # environmental covariates
  psi <- inv.logit(cbind(1, x) %*% beta) # occurrence probability
  y <- rbinom(nsite, 1, psi) # observations of occurrence

  list(x=x, y=y)
} # data_sim

#=======================
# Define MCMC algorithm
#=======================
logistic_mcmc <- function(y, x, nmcmc) {

  # Setup variables
  nsite <- dim(x)[1]
  ncovs <- dim(x)[2]

  beta_save <- matrix(0, ncovs+1, nmcmc)

  # Priors
  beta_mean <- 0
  beta_sd <- 2

  # Starting values
  beta <- rep(0, ncovs + 1)

  # Tuning factor
  beta_tune <- 0.15

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample beta
    beta_star <- rnorm(ncovs+1, beta, beta_tune)
    mh1 <- sum(dbinom(y, 1, inv.logit(cbind(1,x) %*% beta_star), log=T)) + 
           sum(dnorm(beta_star, beta_mean, beta_sd, log=T))
    mh2 <- sum(dbinom(y, 1, inv.logit(cbind(1,x) %*% beta     ), log=T)) + 
           sum(dnorm(beta     , beta_mean, beta_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta <- beta_star
    }

    beta_save[,k] <- beta
  } # k

  # Write output
  list(beta_save=beta_save)
} # logistic_mcmc

#==========
# Run MCMC
#==========
# Check convergence
library(foreach)    # for parallel computing
library(doParallel) # for parallel computing

numCores <- round(detectCores() / 2) # use half of the cores for parallel computing
registerDoParallel(numCores) # setup parallel computing

data <- data_sim(nsite=nsite, ncovs=ncovs, beta=beta) # simulate data
nmcmc <- 50000 # number of iterations
chain <- 3     # number of chains

start_time <- Sys.time() # start time of computing
out <- foreach (i=1:chain, .packages='boot') %dopar% {
  logistic_mcmc(y=data$y, x=data$x, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

pdf(file='1.1 logistic regression_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta <- c(expression(beta[0]), expression(beta[1]), expression(beta[2]), expression(beta[3]))

par(mfrow=c(2,2))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,0))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta[i] - yint * 8
  ymax <- beta[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  if (i %in% 1:2) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.5)
  }
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.5, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta[i], col='grey16', lwd=1.5)
  title(main=ylab_beta[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=1.2, outer=T)

dev.off()

# Evaluate bias
numCores <- round(detectCores() / 2) # only use half of the cores for parallel computing
registerDoParallel(numCores) # setup parallel computing

nsims <- 100   # number of simulations
nmcmc <- 50000 # number of iterations
chain <- 3     # number of chains

set.seed(NULL) # remove seed to generate a different dataset each time
out_nsims <- list()
for (sim in 1:nsims) {
  data <- data_sim(nsite=nsite, ncovs=ncovs, beta=beta)
  out_nsims[[sim]] <- foreach (i=1:chain, .packages='boot') %dopar% {
    logistic_mcmc(y=data$y, x=data$x, nmcmc=nmcmc)
  } # i
} # sim
stopImplicitCluster() # clean up the cluster

pdf(file='1.2 logistic regression_violin.pdf', width=10, height=8)

library(vioplot) # for making violin plots

par(mfrow=c(2,2))
par(mar=c(1,3,3,1))
par(oma=c(0,3,0,0))

for (i in 1:(ncovs+1)) {
  tt <- array(, dim=c(nmcmc*0.8, chain, nsims))
  for (s in 1:nsims) {
    for (j in 1:chain) {
      tt[,j,s] <- out_nsims[[s]][[j]]$beta_save[i,(nmcmc*0.2+1):nmcmc]
    } # j
  } # s
  plot(1, type='n', ylim=c(beta[i]-1,beta[i]+1), xlab='', ylab='', axes=F)
  vioplot(tt, col='lightcoral', rectCol='grey36', lineCol='grey36', border=NA, wex=0.8, cex=1.8, add=T)
  abline(h=beta[i], col='royalblue', lwd=1.5)
  axis(2, at=seq(beta[i]-1,beta[i]+1,length.out=5), cex.axis=1.5, las=2)
  title(main=ylab_beta[i], cex.main=2, line=1.2)
} # i

title(ylab='Posterior Value', cex.lab=2.5, line=1.2, outer=T)

dev.off()


