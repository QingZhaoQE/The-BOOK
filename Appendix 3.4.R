#=============================================================================================================
# Spatial dynamic N-mixture model with population growth & distance-, density- and environment-based movement
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=============================================================================================================

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
npcvs <- 2   # number of observational covariates

beta0 <- c(1.8, 1, -0.8)               # intercept and slopes for initial abundance
beta_rho <- c(0, -0.8, 0.6, -0.4, 0.2) # intercept and slopes for population growth
kappa <- 2                             # distance effect on movement
beta_eta <- c(-0.8, 0.6, 0.2)          # slopes for movement
alpha <- c(0.3, -0.6, 0.4)             # intercept and slopes for detection probability

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

lon <- cov$lon # longitude or easting
lat <- cov$lat # latitude or northing
d <- matrix(, nsite, nsite) # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1

rho <- matrix(0, nsite, nyear-1) # population growth rate
eta <- array(, dim=c(nsite, nsite, nyear-1)) # unstandardized colonization rate
theta <- array(, dim=c(nsite, nsite, nyear-1)) # standardized colonization rate
for (t in 2:nyear) {
  rho[,t-1] <- exp(cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[,t,], (N[,t-1]-lambda[,1])/lambda[,1]*x[,t,1]) %*% beta_rho)
  for (i in 1:nsite) {
    eta[i,,t-1] <- exp(-1 * kappa * d[i,] + cbind((N[,t-1]-N[i,t-1])/exp(beta0[1]), x[,t,1]-x[i,t,1], (N[,t-1]-N[i,t-1])/exp(beta0[1])*(x[,t,1]-x[i,t,1])) %*% beta_eta)
  } # i
  theta[,,t-1] <- eta[,,t-1] / rowSums(eta[,,t-1])
  lambda[,t] <- (N[,t-1] * rho[,t-1]) %*% theta[,,t-1]
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
dist_cov_growth_Nmix_mcmc <- function(y, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:2, max)

  beta0_save <- matrix(0, ncovs+1, nmcmc) 
  beta_rho_save <- matrix(0, ncovs+3, nmcmc)
  kappa_save <- rep(0, nmcmc)
  beta_eta_save <- matrix(0, ncovs+1, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc)
  N_save <- array(0, dim=c(nsite, nyear, 2000))

  # Priors
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_rho_mean <- rep(0, ncovs+3)
  beta_rho_sd <- 2
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  beta_eta_mean <- rep(0, ncovs+1)
  beta_eta_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_rho <- rep(0, ncovs+3)
  kappa <- 1
  beta_eta <- rep(0, ncovs+1)
  alpha <- rep(0, npcvs+1)
  lambda0 <- exp(cbind(1,x[,1,]) %*% beta0) # expectation of initial abundance
  rho <- matrix(1, nsite, nyear-1)
  eta <- theta <- array(, dim=c(nsite, nsite, nyear-1))
  for (t in 1:(nyear-1)) {
    for (i in 1:nsite) {
      eta[i,,t] <- exp(-1 * kappa * d[i,])
    } # i
    theta[,,t] <- eta[,,t] / rowSums(eta[,,t])
  } # t
  N <- round((ymax + 1) / 0.5)
  p <- array(0.5, dim=c(nsite, nyear, nreps))

  # Tuning factors
  beta0_tune <- 0.08
  beta_rho_tune <- 0.08
  kappa_tune <- 0.15
  beta_eta_tune <- 0.1
  alpha_tune <- 0.02
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- matrix(rpois(nsite*nyear, N + N_tune), nsite, nyear)
    prob1 <- prob2 <- matrix(, nsite, nyear)
    prob1[,1] <- dpois(N_star[,1], lambda0, log=T)
    prob2[,1] <- dpois(N     [,1], lambda0, log=T)
    for (t in 2:nyear) {
      prob1[,t] <- dpois(N_star[,t], (N[,t-1] * rho[,t-1]) %*% theta[,,t-1], log=T)
      prob2[,t] <- dpois(N     [,t], (N[,t-1] * rho[,t-1]) %*% theta[,,t-1], log=T)
    } # t
    mh1 <- apply(dbinom(y, N_star, p, log=T), 1:2, sum) + 
           prob1 + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , p, log=T), 1:2, sum) + 
           prob2 + 
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
    beta_rho_star <- rnorm(ncovs+3, beta_rho, beta_rho_tune)
    rho_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      rho_star[,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,],(N[,t-1]-lambda0)/lambda0*x[,t,1]) %*% beta_rho_star)
    } # t
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[,t] <- dpois(N[,t+1], (N[,t] * rho_star[,t]) %*% theta[,,t], log=T)
      prob2[,t] <- dpois(N[,t+1], (N[,t] * rho     [,t]) %*% theta[,,t], log=T)
    } # t
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_rho_star, beta_rho_mean, beta_rho_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta_rho     , beta_rho_mean, beta_rho_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_rho <- beta_rho_star
      rho <- rho_star
    }

    ### Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    eta_star <- theta_star <- array(, dim=c(nsite, nsite, nyear-1))
    for (t in 2:nyear) {
      for (i in 1:nsite) {
        eta_star[i,,t-1] <- exp(-1 * kappa_star * d[i,] + cbind((N[,t-1]-N[i,t-1])/exp(beta0[1]),x[,t,1]-x[i,t,1],(N[,t-1]-N[i,t-1])/exp(beta0[1])*(x[,t,1]-x[i,t,1])) %*% beta_eta)
      } # i
      theta_star[,,t-1] <- eta_star[,,t-1] / rowSums(eta_star[,,t-1])
    } # t
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[,t] <- dpois(N[,t+1], (N[,t] * rho[,t]) %*% theta_star[,,t], log=T)
      prob2[,t] <- dpois(N[,t+1], (N[,t] * rho[,t]) %*% theta     [,,t], log=T)
    } # t
    mh1 <- sum(prob1) + 
           dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=T)
    mh2 <- sum(prob2) + 
           dnorm(log(kappa     ), log_kappa_mean, log_kappa_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
      eta <- eta_star
      theta <- theta_star
    }

    ### Sample beta_eta
    beta_eta_star <- rnorm(ncovs+1, beta_eta, beta_eta_tune)
    eta_star <- theta_star <- array(, dim=c(nsite, nsite, nyear-1))
    for (t in 2:nyear) {
      for (i in 1:nsite) {
        eta_star[i,,t-1] <- exp(-1 * kappa * d[i,] + cbind((N[,t-1]-N[i,t-1])/exp(beta0[1]),x[,t,1]-x[i,t,1],(N[,t-1]-N[i,t-1])/exp(beta0[1])*(x[,t,1]-x[i,t,1])) %*% beta_eta_star)
      }
      theta_star[,,t-1] <- eta_star[,,t-1] / rowSums(eta_star[,,t-1])
    } # t
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[,t] <- dpois(N[,t+1], (N[,t] * rho[,t]) %*% theta_star[,,t], log=T)
      prob2[,t] <- dpois(N[,t+1], (N[,t] * rho[,t]) %*% theta     [,,t], log=T)
    } # t
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_eta_star, beta_eta_mean, beta_eta_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta_eta     , beta_eta_mean, beta_eta_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_eta <- beta_eta_star
      eta <- eta_star
      theta <- theta_star
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
    kappa_save[k] <- kappa
    beta_eta_save[,k] <- beta_eta
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        N_save[,,k-nmcmc+2000] <- N
      }
    }
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_rho_save=beta_rho_save, kappa_save=kappa_save, beta_eta_save=beta_eta_save, 
       alpha_save=alpha_save, 
       N_save=N_save)

} # dist_cov_growth_Nmix_mcmc

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
  dist_cov_growth_Nmix_mcmc(y=y, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.3.2 dist cov population growth N-mixture model_output.RData')

#==============
# Plot results
#==============
pdf(file='3.3.2.1 dist cov population growth N-mixture model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["pond"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_rho <- c(expression(beta[0]^""["["*rho*"]"]), expression(beta["density"]^""["["*rho*"]"]), 
                   expression(beta["pond"]^""["["*rho*"]"]), expression(beta["temperature"]^""["["*rho*"]"]), 
                   expression(beta["density "%*%" pond"]^""["["*rho*"]"]))
ylab_beta_eta <- c(expression(beta["density"]^""["["*eta*"]"]), expression(beta["pond"]^""["["*eta*"]"]), 
                   expression(beta["density "%*%" pond"]^""["["*eta*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))

par(mfrow=c(5,3))
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

for (i in 1:(ncovs+3)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_rho_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_rho[i] - yint * c(8,4,8,4,8)[i]
  ymax <- beta_rho[i] + yint * c(8,4,8,4,8)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(4,2,4,2,4)[i]), digits=1), cex.axis=1.1, las=2)
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
  tt[,j] <- out[[j]]$kappa_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa - yint * 8
ymax <- kappa + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa, col='grey16', lwd=1.5)
title(main=expression(kappa), cex.main=2, line=1.1)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_eta_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_eta[i] - yint * c(8,8,4)[i]
  ymax <- beta_eta[i] + yint * c(8,8,4)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(4,4,2)[i]), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_eta[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_eta[i], cex.main=2, line=1.2)
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
  ymin <- alpha[i] - yint * 2
  ymax <- alpha[i] + yint * 2
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
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

pdf(file='3.3.2.2 dist cov population growth N-mixture model_theta.pdf', width=9, height=8)

kappa_post <- c(out[[1]]$kappa[(nmcmc/2+1):nmcmc], out[[2]]$kappa[(nmcmc/2+1):nmcmc], out[[3]]$kappa[(nmcmc/2+1):nmcmc])
beta0_post <- cbind(out[[1]]$beta0[,(nmcmc/2+1):nmcmc], out[[2]]$beta0[,(nmcmc/2+1):nmcmc], out[[3]]$beta0[,(nmcmc/2+1):nmcmc])
beta_eta_post <- cbind(out[[1]]$beta_eta[,(nmcmc/2+1):nmcmc], out[[2]]$beta_eta[,(nmcmc/2+1):nmcmc], out[[3]]$beta_eta[,(nmcmc/2+1):nmcmc])
N_post <- array(, dim=c(nsite, nyear, 6000))
N_post[,,   1:2000] <- out[[1]]$N
N_post[,,2001:4000] <- out[[2]]$N
N_post[,,4001:6000] <- out[[3]]$N
kappa_med <- median(kappa_post)
beta0_med <- apply(beta0_post, 1, median)
beta_eta_med <- apply(beta_eta_post, 1, median)
N_med <- apply(N_post, 1:2, median)

years_plot <- c(7, 14)
eta_true <- theta_true <- eta_pred <- theta_pred <- array(, dim=c(nsite, nsite, length(years_plot)))
for (t in 1:length(years_plot)) {
  yr <- years_plot[t]
  for (i in 1:nsite) {
    eta_true[i,,t] <- exp(-1 * kappa     * d[i,] + cbind((N    [,yr]-N    [i,yr])/exp(beta0    [1]),x[,yr,1]-x[i,yr,1],(N    [,yr]-N    [i,yr])/exp(beta0    [1])*(x[,yr,1]-x[i,yr,1])) %*% beta_eta    )
    eta_pred[i,,t] <- exp(-1 * kappa_med * d[i,] + cbind((N_med[,yr]-N_med[i,yr])/exp(beta0_med[1]),x[,yr,1]-x[i,yr,1],(N_med[,yr]-N_med[i,yr])/exp(beta0_med[1])*(x[,yr,1]-x[i,yr,1])) %*% beta_eta_med)
  } # i
  theta_true[,,t] <- eta_true[,,t] / rowSums(eta_true[,,t])
  theta_pred[,,t] <- eta_pred[,,t] / rowSums(eta_pred[,,t])
} # t

plot_cut <- 0.1

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(4,5,5,2))

for (t in 1:length(years_plot)) {
  plot(cov$lon, cov$lat, xlim=c(-10.2,10.2), ylim=c(-8.2,8.2), axes=F, xlab='', ylab='', type='n')
  axis(1, at=seq(-10,10,5), labels=rep('',5))
  if (t == 1) {
    axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
  } else {
    axis(2, at=seq(-08,08,4), labels=rep('',5))
    axis(4, at=0, labels='True', cex.axis=2.5, tick=F, line=.5)
  }
  axis(3, at=0, labels=paste('Year ', years_plot[t], sep=''), cex.axis=2.5, line=-.8, tick=F)
  box()
  for (i in 1:(nsite-1)) {
    for (j in (i+1):nsite) {
      length <- ifelse(d[i,j] / 8 > 0.1, 0.1, d[i,j] / 8)
      if (theta_true[i,j,t] > theta_true[j,i,t]) {
        if ((theta_true[i,j,t] - theta_true[j,i,t]) > plot_cut) {
          arrows(x0=cov$lon[j], y0=cov$lat[j], x1=cov$lon[i], y1=cov$lat[i], 
                 lwd=(theta_true[i,j,t] - theta_true[j,i,t])*5, length=length, col='royalblue')
        }
      } else {
        if ((theta_true[j,i,t] - theta_true[i,j,t]) > plot_cut) {
          arrows(x0=cov$lon[i], y0=cov$lat[i], x1=cov$lon[j], y1=cov$lat[j], 
                 lwd=(theta_true[j,i,t] - theta_true[i,j,t])*5, length=length, col='royalblue')
        }
      }
    } # j
  } # i
  points(cov$lon, cov$lat, pch=21, cex=rowMeans(cov$size)*4.6, col='grey18', bg=NA)
} # t

for (t in 1:length(years_plot)) {
  plot(cov$lon, cov$lat, xlim=c(-10.2,10.2), ylim=c(-8.2,8.2), axes=F, xlab='', ylab='', type='n')
  axis(1, at=seq(-10,10,5), cex.axis=1.6)
  if (t == 1) {
    axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
  } else {
    axis(2, at=seq(-08,08,4), labels=rep('',5))
    axis(4, at=0, labels='Estimated', cex.axis=2.5, tick=F, line=.5)
  }
  box()
  for (i in 1:(nsite-1)) {
    for (j in (i+1):nsite) {
      length <- ifelse(d[i,j] / 8 > 0.1, 0.1, d[i,j] / 8)
      if (theta_pred[i,j,t] > theta_pred[j,i,t]) {
        if ((theta_pred[i,j,t] - theta_pred[j,i,t]) > plot_cut) {
          arrows(x0=cov$lon[j], y0=cov$lat[j], x1=cov$lon[i], y1=cov$lat[i], 
                 lwd=(theta_pred[i,j,t] - theta_pred[j,i,t])*5, length=length, col='lightcoral')
        }
      } else {
        if ((theta_pred[j,i,t] - theta_pred[i,j,t]) > plot_cut) {
          arrows(x0=cov$lon[i], y0=cov$lat[i], x1=cov$lon[j], y1=cov$lat[j], 
                 lwd=(theta_pred[j,i,t] - theta_pred[i,j,t])*5, length=length, col='lightcoral')
        }
      }
    } # j
  } # i
  points(cov$lon, cov$lat, pch=21, cex=rowMeans(cov$size)*4.6, col='grey18', bg=NA)
} # t

title(xlab='Easting' , line=2.6, cex.lab=2.8, outer=T)
title(ylab='Northing', line=2.4, cex.lab=2.8, outer=T)
title(main=expression(paste('Colonization Rate (', theta, ')', sep='')), line=3.2, cex.main=3, outer=T)

dev.off()


