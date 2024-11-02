#=================================================================================================
# Static N-mixture model
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

#setwd('')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(5)

# Basic values
nsite <- 360 # number of sites
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 2   # number of observational covariates

beta <- c(1.8, 1, -0.8)    # intercept and slopes for true abundance
alpha <- c(0.3, -0.6, 0.3) # intercept and slopes for detection probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- cbind(x1[,1], x2[,2]) # environmental covariates

# Simulate data
lambda <- exp(cbind(1,x) %*% beta) # expectation of abundance
N <- rpois(nsite, lambda) # true abundance

w <- array(rnorm(nsite * nreps * npcvs, 0, 1), dim=c(nsite, nreps, npcvs)) # observational covariates

p <- matrix(, nsite, nreps) # detection probability
for (i in 1:nsite) {
  p[i,] <- inv.logit(cbind(1,w[i,,]) %*% alpha)
} # i

y <- matrix(, nsite, nreps) # count data
for (i in 1:nsite) {
  y[i,] <- rbinom(nreps, N[i], p[i,])
} # i

#=======================
# Define MCMC algorithm
#=======================
Nmix_mcmc <- function(y, x, w, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nreps <- dim(y)[2]
  ncovs <- dim(x)[2]
  npcvs <- dim(w)[3]
  ymax <- apply(y, 1, max)

  beta_save <- matrix(0, ncovs + 1, nmcmc)
  alpha_save <- matrix(0, npcvs + 1, nmcmc)
  N_save <- matrix(0, nsite, 2000)

  # Priors
  beta_mean <- rep(0, ncovs + 1)
  beta_sd <- 2
  alpha_mean <- rep(0, npcvs + 1)
  alpha_sd <- 2

  # Starting values
  beta <- rep(0, ncovs + 1)
  alpha <- rep(0, npcvs + 1)
  lambda <- exp(cbind(1,x) %*% beta)
  N <- round((ymax + 1) / 0.5)
  p <- matrix(0.5, nsite, nreps)

  # Tuning factors
  beta_tune <- 0.15
  alpha_tune <- c(0.3, 0.1, 0.1)
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- rpois(nsite, N + N_tune)
    mh1 <- apply(dbinom(y, N_star, p, log=T), 1, sum) + 
           dpois(N_star, lambda, log=T) + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , p, log=T), 1, sum) + 
           dpois(N     , lambda, log=T) + 
           dpois(N_star, N+N_tune, log=T)
    mh <- exp(mh1 - mh2)
    Nkeep <- ((mh > runif(nsite)) & (N_star >= ymax))
    N[Nkeep] <- N_star[Nkeep]

    ### Sample beta
    beta_star <- rnorm(ncovs + 1, beta, beta_tune)
    lambda_star <- exp(cbind(1,x) %*% beta_star)
    mh1 <- sum(dpois(N, lambda_star, log=T)) + 
           sum(dnorm(beta_star, beta_mean, beta_sd, log=T))
    mh2 <- sum(dpois(N, lambda     , log=T)) + 
           sum(dnorm(beta     , beta_mean, beta_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta <- beta_star
      lambda <- lambda_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs + 1, alpha, alpha_tune)
    p_star <- matrix(, nsite, nreps)
    for (i in 1:nsite) {
      p_star[i,] <- inv.logit(cbind(1,w[i,,]) %*% alpha_star)
    } # i
    mh1 <- sum(dbinom(y, N, p_star, log=TRUE)) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(dbinom(y, N, p     , log=TRUE)) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Save samples
    beta_save[,k] <- beta
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        N_save[,k-nmcmc+2000] <- N
      }
    }
  } # k

  # Write output
  list(beta_save=beta_save, alpha_save=alpha_save, N_save=N_save)
} # Nmix_mcmc

#==========
# Run MCMC
#==========
library(foreach)    # for parallel computing
library(doParallel) # for parallel computing

numCores <- round(detectCores() / 2) # only use half of the cores for parallel computing
registerDoParallel(numCores) # setup parallel computing

nmcmc <- 50000 # number of iterations of MCMC computing
chain <- 3     # number of chains

start_time <- Sys.time() # start time of computing
out <- foreach (i=1:chain, .packages='boot') %dopar% {
  Nmix_mcmc(y=y, x=x, w=w, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.1 N-mixture model_output.RData')

#==============
# Plot results
#==============
pdf(file='3.1.1 N-mixture model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta <- c(expression(beta[0]), expression(beta["pond"]), expression(beta["temperature"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]))

par(mfrow=c(2,3))
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
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta[i], col='grey16', lwd=1.5)
  title(main=ylab_beta[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i] - yint * 4
  ymax <- alpha[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1)
  text(x=nmcmc*0.6, y=ymax, pos=1, cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=0.4, outer=T)

dev.off()

pdf(file='3.1.2 N-mixture model_N.pdf', width=8.5, height=8)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    N_post <- out[[i]]$N_save
  } else {
    N_post <- cbind(N_post, out[[i]]$N_save)
  }
} # i
N_med <- apply(N_post, 1, median)

nscale <- 100
Nmax <- max(c(log(N+1), log(N_med+1)))
N_scale <- round(log(N+1) / Nmax * (nscale-1)) + 1
N_med_scale <- round(log(N_med+1) / Nmax * (nscale-1)) + 1

colt <- colorRampPalette(c('white','lightseagreen','navy'))(nscale)
cole <- colorRampPalette(c('white','gold','firebrick'))(nscale)

layout(mat=matrix(c(rep(c(1,1,2,2),times=3),rep(3,8)), nrow=4, ncol=5, byrow=FALSE))
par(mar=c(1,2,3,1))
par(oma=c(4,4,0,6))

plot(cov$lon, cov$lat, xlim=c(-10,10), ylim=c(-8,8), axes=F, xlab='', ylab='', 
     pch=21, cex=cov$size[,1]*6, col='grey18', bg=colt[N_scale])
axis(1, at=seq(-10,10,5), labels=rep('',5))
axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
box()
title(main='True Abundance', cex.main=2)

plot(cov$lon, cov$lat, xlim=c(-10,10), ylim=c(-8,8), axes=F, xlab='', ylab='', 
     pch=21, cex=cov$size[,1]*6, col='grey18', bg=cole[N_med_scale])
axis(1, at=seq(-10,10,5), cex.axis=1.6)
axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
box()
title(main='Estimated Abundance', cex.main=2)

axis(1, at=0, labels='Easting', tick=F, line=2.5, cex.axis=3)
title(ylab='Northing', cex.lab=3, line=1.5, outer=T)

par(mar=c(1,1,3,2))
plot(1, ylim=c(1250, 1500), type='n', axes=F, xlab='', ylab='')
vioplot(apply(N_post,2,sum), col='lightcoral', rectCol='grey36', lineCol='grey36', 
        border=NA, cex=3, add=T)
abline(h=sum(N), col='royalblue', lwd=1.8)
axis(4, at=seq(1250, 1500, 50), cex.axis=1.6, las=2)
axis(4, at=1375, label='Total Abundance', line=5.5, cex.axis=3, tick=F)

dev.off()


