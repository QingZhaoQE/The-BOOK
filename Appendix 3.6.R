#=================================================================================================
# Dail-Madsen model without age structure
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

#setwd('')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(6)

# Basic values
nsite <- 360 # number of sites
nyear <- 20  # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 3   # number of observational covariates

beta0 <- c(1.8, 1, -0.8)              # intercept and slopes for the initial abundance
beta_phi   <- c(0.8, -0.4, 0.4, -0.2) # intercept and slopes for apparent survival
beta_gamma <- c(0.4, -0.6, 0.6, -0.3) # intercept and slopes for gain (recruitment)
alpha <- c(0.3, -0.6, 0.3, -0.3)      # intercept and slopes for detection

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

phi   <- matrix(, nsite, nyear-1) # apparent survival
gamma <- matrix(, nsite, nyear-1) # recrutiment
for (t in 2:nyear) {
  phi  [,t-1] <- inv.logit(cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[,t,]) %*% beta_phi  )
  gamma[,t-1] <- exp      (cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[,t,]) %*% beta_gamma)
  lambda[,t] <- N[,t-1] * phi[,t-1] + gamma[,t-1]
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
Dail_Madsen_mcmc <- function(y, x, w, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:2, max)

  beta0_save <- matrix(0, ncovs+1, nmcmc) 
  beta_phi_save   <- matrix(0, ncovs+2, nmcmc)
  beta_gamma_save <- matrix(0, ncovs+2, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc) 
  N_save <- matrix(0, nyear, nmcmc)
  phi_save   <- matrix(0, nyear-1, nmcmc)
  gamma_save <- matrix(0, nyear-1, nmcmc)

  # Priors
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_phi_mean <- rep(0, ncovs+2)
  beta_phi_sd <- 2
  beta_gamma_mean <- rep(0, ncovs+2)
  beta_gamma_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_phi   <- rep(0, ncovs+2)
  beta_gamma <- rep(0, ncovs+2)
  alpha <- rep(0, npcvs+1)
  lambda0 <- rep(1, nsite) # expectation of initial abundance
  phi   <- matrix(0.5, nsite, nyear-1)
  gamma <- matrix(1  , nsite, nyear-1)
  N <- round((ymax + 1) / 0.5)
  p <- array(0.5, dim=c(nsite, nyear, nreps))

  # Tuning factors
  beta0_tune <- 0.05
  beta_phi_tune   <- c(0.15, 0.10, 0.10, 0.10)
  beta_gamma_tune <- c(0.08, 0.05, 0.05, 0.05)
  alpha_tune <- 0.01
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- matrix(rpois(nsite*nyear, N + N_tune), nsite, nyear)
    mh1 <- apply(dbinom(y, N_star, p, log=T), 1:2, sum) + 
           dpois(N_star, cbind(lambda0, N[,-nyear] * phi + gamma), log=T) + 
           dpois(N, N_star + N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , p, log=T), 1:2, sum) + 
           dpois(N     , cbind(lambda0, N[,-nyear] * phi + gamma), log=T) + 
           dpois(N_star, N + N_tune, log=T)
    mhN <- exp(mh1 - mh2)
    Nkeep <- ((mhN > runif(nsite*nyear)) & (N_star >= ymax))
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

    ### Sample beta_phi
    beta_phi_star <- rnorm(ncovs+2, beta_phi, beta_phi_tune)
    phi <- phi_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi     )
      phi_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_star)
    } # t
    mh1 <- sum(dpois(N[,-1], N[,-nyear] * phi_star + gamma, log=T)) + 
           sum(dnorm(beta_phi_star, beta_phi_mean, beta_phi_sd, log=T))
    mh2 <- sum(dpois(N[,-1], N[,-nyear] * phi      + gamma, log=T)) + 
           sum(dnorm(beta_phi     , beta_phi_mean, beta_phi_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_phi <- beta_phi_star
    }

    ### Sample beta_gamma
    beta_gamma_star <- rnorm(ncovs+2, beta_gamma, beta_gamma_tune)
    gamma <- gamma_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      gamma     [,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma     )
      gamma_star[,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma_star)
    } # t
    mh1 <- sum(dpois(N[,-1], N[,-nyear] * phi + gamma_star, log=T)) + 
           sum(dnorm(beta_gamma_star, beta_gamma_mean, beta_gamma_sd, log=T))
    mh2 <- sum(dpois(N[,-1], N[,-nyear] * phi + gamma     , log=T)) + 
           sum(dnorm(beta_gamma     , beta_gamma_mean, beta_gamma_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_gamma <- beta_gamma_star
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
    beta_phi_save[,k] <- beta_phi
    beta_gamma_save[,k] <- beta_gamma
    alpha_save[,k] <- alpha
    N_save[,k] <- colMeans(N)
    phi_save[,k] <- colMeans(phi)
    gamma_save[,k] <- colMeans(gamma)
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_phi_save=beta_phi_save, beta_gamma_save=beta_gamma_save, 
       alpha_save=alpha_save, 
       N_save=N_save, phi_save=phi_save, gamma_save=gamma_save)

} # Dail_Madsen_mcmc

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
  Dail_Madsen_mcmc(y=y, x=x, w=w, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.5.1 Dail-Madsen model_output.RData')

#==============
# Plot results
#==============
pdf(file='3.5.1.1 Dail-Madsen model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["pond"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["density"]^""["["*phi*"]"]), 
                   expression(beta["pond"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta["density"]^""["["*gamma*"]"]), 
                     expression(beta["pond"]^""["["*gamma*"]"]), expression(beta["temperature"]^""["["*gamma*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]), 
                expression(alpha[3]), expression(alpha[4]))

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
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi[i] - yint * c(8,4,8,8)[i]
  ymax <- beta_phi[i] + yint * c(8,4,8,8)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(4,2,4,4)[i]), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_gamma_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma[i] - yint * c(4,8,4,4)[i]
  ymax <- beta_gamma[i] + yint * c(4,8,4,4)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(2,4,2,2)[i]), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.6, 
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
  if (i == 1) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1)
  text(x=nmcmc*0.6, y=ymax, pos=1, cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=0.4, outer=T)

dev.off()

pdf(file='3.5.1.2 Dail-Madsen model_N phi gamma.pdf', width=8, height=8)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    N_mean <- out[[1]]$N_save
    phi_mean <- out[[1]]$phi_save
    gamma_mean <- out[[1]]$gamma_save
  } else {
    N_mean <- cbind(N_mean, out[[1]]$N_save)
    phi_mean <- cbind(phi_mean, out[[1]]$phi_save)
    gamma_mean <- cbind(gamma_mean, out[[1]]$gamma_save)
  }
} # i

par(mfrow=c(3,1))
par(mar=c(1.5,1.5,0,0))
par(oma=c(4,5,1,1))

plot(1, xlim=c(1,nyear), ylim=c(0,20), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(N), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), labels=rep('',5))
axis(2, at=seq(0,20,4), cex.axis=1.6, las=2)
axis(2, at=10, labels='Population Size', cex.axis=2.5, line=3.2, tick=F)

points(x=7, y=18, pch=16, cex=1, col='royalblue')
lines(x=c(6.5,7.5), y=c(18,18), lwd=1.2, col='royalblue')
vioplot(rnorm(10000,18,0.5), at=10.5, add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
text(x=7.5, y=18, labels='True', pos=4, cex=2.5)
text(x=11, y=18, labels='Estimated', pos=4, cex=2.5)

plot(1, xlim=c(1,nyear), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(phi), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), labels=rep('',5))
axis(2, at=seq(0,1,0.2), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Apparent Survival', cex.axis=2.5, line=3.2, tick=F)

plot(1, xlim=c(1,nyear), ylim=c(0,5), type='l', axes=F, xlab='', ylab='')
vioplot(t(gamma_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(gamma), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), cex.axis=1.6)
axis(2, at=seq(0,5,1), cex.axis=1.6, las=2)
axis(2, at=2.5, labels='Recruitment', cex.axis=2.5, line=3.2, tick=F)

title(xlab='Year', cex.lab=3, line=2.4, outer=T)

dev.off()

pdf(file='3.5.1.3 Dail-Madsen model_correlation.pdf', width=8, height=8)

for (i in 1:chain) {
  if (i == 1) {
    beta_phi0_post   <- out[[i]]$beta_phi  [1,(nmcmc/2+1):nmcmc]
    beta_gamma0_post <- out[[i]]$beta_gamma[1,(nmcmc/2+1):nmcmc]
  } else {
    beta_phi0_post   <- c(beta_phi0_post  , out[[i]]$beta_phi  [1,(nmcmc/2+1):nmcmc])
    beta_gamma0_post <- c(beta_gamma0_post, out[[i]]$beta_gamma[1,(nmcmc/2+1):nmcmc])
  }
} # i

nplot <- 10000
sel <- sort(sample(1:(nmcmc/2*chain), nplot, replace=F))
a <- beta_phi0_post[sel]
b <- beta_gamma0_post[sel]
a2 <- (a - min(a)) / (max(a) - min(a))
b2 <- (b - min(b)) / (max(b) - min(b))
ab <- matrix(, nplot, nplot)
for (i in 1:nplot) {
  for (j in 1:nplot) {
    ab[i,j] <- sqrt((a2[i] - a2[j]) ^ 2 + (b2[i] - b2[j]) ^ 2)
  } # j
} # i

ab01 <- ifelse(ab < 0.08, 1, 0)
ab_scale <- rowSums(ab01) - min(rowSums(ab01)) + 1
col2 <- c(colorRampPalette(c('navy','seagreen'))(floor(max(ab_scale)*0.6)), 
          colorRampPalette(c('seagreen','yellow'))(ceiling(max(ab_scale)*0.4)))

par(mfrow=c(1,1))
par(mar=c(5,6,1,1))

plot(a, b, col=col2[ab_scale], pch=16, cex=1.4, xlim=c(0.2,1.4), ylim=c(0.1,0.7), axes=F, xlab='', ylab='')
axis(1, at=seq(0.2,1.6,0.4), cex.axis=1.5)
axis(2, at=seq(0.1,0.8,0.2), cex.axis=1.5, las=2)
axis(1, at=0.8, labels=ylab_beta_phi  [1], cex.axis=2, line=2.5, tick=F)
axis(2, at=0.4, labels=ylab_beta_gamma[1], cex.axis=2, line=2.8, las=2, tick=F)
box()

text(x=0.84, y=0.68, labels=paste('correlation = ', round(cor(a, b), digits=2), sep=''), pos=4, cex=1.8)

dev.off()


