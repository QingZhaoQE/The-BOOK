#=================================================================================================
# Dynamic N-mixture model with population growth (non-spatial)
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
nyear <- 20  # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 3   # number of observational covariates

beta0 <- c(1.8, 1, -0.8)          # intercept and slopes for initial abundance
beta_rho <- c(0, -0.8, 0.6, -0.4) # intercept and slopes for population growth
delta <- 0.5                      # expectation of the number of immigrants
alpha <- c(0.3, -0.6, 0.3, -0.3)  # intercept and slopes for detection probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- array(, dim=c(nsite, nyear, ncovs)) # environmental covariates
x[,,1] <- x1
x[,,2] <- x2
### Standardize covariates for each site
for (i in 1:nsite) {
  for (h in 1:ncovs) {
    x[i,2:nyear,h] <- (x[i,2:nyear,h] - mean(x[i,,h])) / sd(x[i,,h])
  } # h
} # i

# Simulate data
lambda <- matrix(, nsite, nyear) # expectation of abundance
N <- matrix(, nsite, nyear) # abundance
lambda[,1] <- exp(cbind(1,x[,1,]) %*% beta0) # expectation of initial abundance
N[,1] <- rpois(nsite, lambda[,1]) # initial abundance

rho <- matrix(0, nsite, nyear-1) # population growth rate
for (t in 2:nyear) {
  rho[,t-1] <- exp(cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[,t,]) %*% beta_rho)
  lambda[,t] <- N[,t-1] * rho[,t-1] + delta
  N[,t] <- rpois(nsite, lambda[,t])
} # t

w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

p <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    p[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha)
  } # t
} # i

y <- array(, dim=c(nsite, nyear, nreps)) # count data
for (i in 1:nsite) {
  for (t in 1:nyear) {
    y[i,t,] <- rbinom(nreps, N[i,t], p[i,t,])
  } # t
} # i

#=======================
# Define MCMC algorithm
#=======================
growth_Nmix_mcmc <- function(y, x, w, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:2, max)

  beta0_save <- matrix(0, ncovs+1, nmcmc) 
  beta_rho_save <- matrix(0, ncovs+2, nmcmc)
  delta_save <- rep(0, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc) 
  Nsum_save <- matrix(0, nyear, nmcmc)

  # Priors
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_rho_mean <- rep(0, ncovs+2)
  beta_rho_sd <- 2
  log_delta_mean <- 0
  log_delta_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_rho <- rep(0, ncovs+2)
  delta <- 1
  alpha <- rep(0, npcvs+1)
  lambda0 <- rep(1, nsite) # expectation of initial abundance
  rho <- matrix(1, nsite, nyear-1)
  N <- round((ymax + 1) / 0.5)
  p <- array(0.5, dim=c(nsite, nyear, nreps))

  # Tuning factors
  beta0_tune <- 0.08
  beta_rho_tune <- 0.15
  delta_tune <- 0.15
  alpha_tune <- 0.01
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- matrix(rpois(nsite*nyear, N + N_tune), nsite, nyear)
    mh1 <- apply(dbinom(y, N_star, p, log=T), 1:2, sum) + 
           dpois(N_star, cbind(lambda0, N[,-nyear]*rho + delta), log=T) + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , p, log=T), 1:2, sum) + 
           dpois(N     , cbind(lambda0, N[,-nyear]*rho + delta), log=T) + 
           dpois(N_star, N+N_tune, log=T)
    mh <- exp(mh1 - mh2)
    Nkeep <- ((mh > runif(nsite*nyear)) & (N_star >= ymax))
    N[Nkeep] <- N_star[Nkeep]

    ### Sample beta0
    beta0_star <- rnorm(ncovs+1, beta0, beta0_tune)
    lambda0_star <- exp(cbind(1,x[,1,]) %*% beta0_star)
    mh1 <- sum(dpois(N[,1], lambda0_star, log=T)) + 
           sum(dnorm(beta0_star, beta0_mean, beta0_sd, log=T))
    mh2 <- sum(dpois(N[,1], lambda0     , log=T)) + 
           sum(dnorm(beta0     , beta0_mean, beta0_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta0 <- beta0_star
      lambda0 <- lambda0_star
    }

    ### Sample beta_rho
    beta_rho_star <- rnorm(ncovs+2, beta_rho, beta_rho_tune)
    rho_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      rho_star[,t-1] <- exp(cbind(1, (N[,t-1]-lambda0)/lambda0, x[,t,]) %*% beta_rho_star)
    } # t
    mh1 <- sum(dpois(N[,-1], N[,-nyear] * rho_star + delta, log=T)) + 
           sum(dnorm(beta_rho_star, beta_rho_mean, beta_rho_sd, log=T))
    mh2 <- sum(dpois(N[,-1], N[,-nyear] * rho      + delta, log=T)) + 
           sum(dnorm(beta_rho     , beta_rho_mean, beta_rho_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_rho <- beta_rho_star
      rho <- rho_star
    }

    ### Sample delta
    delta_star <- exp(rnorm(1, log(delta), delta_tune))
    mh1 <- sum(dpois(N[,-1], N[,-nyear] * rho + delta_star, log=T)) + 
           dnorm(log(delta_star), log_delta_mean, log_delta_sd, log=T)
    mh2 <- sum(dpois(N[,-1], N[,-nyear] * rho + delta     , log=T)) + 
           dnorm(log(delta     ), log_delta_mean, log_delta_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      delta <- delta_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs+1, alpha, alpha_tune)
    p_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        p_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_star)
      } # t
    } # i
    mh1 <- sum(dbinom(y, N, p_star, log=TRUE), na.rm=T) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(dbinom(y, N, p     , log=TRUE), na.rm=T) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Save samples
    beta0_save[,k] <- beta0
    beta_rho_save[,k] <- beta_rho
    delta_save[k] <- delta
    alpha_save[,k] <- alpha
    Nsum_save[,k] <- colSums(N)
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_rho_save=beta_rho_save, delta_save=delta_save, 
       alpha_save=alpha_save, 
       Nsum_save=Nsum_save)

} # growth_Nmix_mcmc

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
  growth_Nmix_mcmc(y=y, x=x, w=w, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.2 population growth N-mixture model_output.RData')

#==============
# Plot results
#==============
pdf(file='3.2.1 population growth N-mixture model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["pond"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_rho <- c(expression(beta[0]^""["["*rho*"]"]), expression(beta["density"]^""["["*rho*"]"]), 
                   expression(beta["pond"]^""["["*rho*"]"]), expression(beta["temperature"]^""["["*rho*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))

par(mfrow=c(4,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,0))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta0_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0[i] - yint * 4
  ymax <- beta0[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i], cex.main=2, line=1.4)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_rho_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_rho[i] - yint * 8
  ymax <- beta_rho[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_rho[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_rho[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$delta_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- delta - yint * 16
ymax <- delta + yint * 16
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*8), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=delta, col='grey16', lwd=1.5)
title(main=expression(delta), cex.main=2, line=1.1)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i] - yint * c(4,2,2,2)[i]
  ymax <- alpha[i] + yint * c(4,2,2,2)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  if (i == 1) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*c(2,1,1,1)[i]), digits=1), cex.axis=1.1, las=2)
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

pdf(file='3.2.2 population growth N-mixture model_Nsum.pdf', width=10, height=8)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    Nsum_post <- out[[i]]$Nsum_save[,(nmcmc/2+1):nmcmc]
  } else {
    Nsum_post <- cbind(Nsum_post, out[[i]]$Nsum_save[,(nmcmc/2+1):nmcmc])
  }
} # i

par(mfrow=c(1,1))
par(mar=c(5,7,1,1))

plot(1, xlim=c(1,nyear), ylim=c(1000,4200), type='n', xlab='', ylab='', axes=F)
vioplot(t(Nsum_post), add=T, col='lightcoral', rectCol='grey36', lineCol='grey36', border=NA)
lines(colSums(N) ~ c(1:nyear), type='o', pch=16, cex=0.8, lwd=1.2, col='royalblue')
axis(1, at=seq(0,nyear,5), cex.axis=1.5)
axis(2, at=seq(1000,4000,1000), cex.axis=1.5, las=2)
title(xlab='Year', cex.lab=2.5, line=3.5)
title(ylab='Total Population Size', cex.lab=2.5, line=4.8)

points(x=5, y=3944, pch=16, cex=1, col='royalblue')
lines(x=c(4.4,5.6), y=c(3944,3944), lwd=1.2, col='royalblue')
vioplot(Nsum_post[8,]-median(Nsum_post[8,])+3944, at=10, add=T, col='lightcoral', rectCol='grey36', lineCol='grey36', border=NA)
text(x= 5.6, y=3944, labels='True', pos=4, cex=1.8)
text(x=10.6, y=3944, labels='Estimated', pos=4, cex=1.8)

dev.off()


