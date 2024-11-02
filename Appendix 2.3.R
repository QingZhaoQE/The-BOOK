#=================================================================================================
# Spatial dynamic occupancy model with adjacency-based movement
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

#setwd('')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(3)

# Basic values
nsite <- 360 # number of sites
nyear <- 20  # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 3   # number of observational covariates

beta0 <- c(0.4, 1, -0.4)         # intercept and slopes for initial occupancy
beta_phi <- c(-0.2, 0.6, -0.3)   # intercept and slopes for persistence
zeta <- 0.4                      # probability of neighborhood colonization
delta <- 0.1                     # probability of long-distance colonization
alpha <- c(0.3, -0.6, 0.4, -0.2) # intercept and slopes for detection probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- array(, dim=c(nsite, nyear, ncovs)) # environmental covariates
x[,,1] <- x1
x[,,2] <- x2

# Simulate data
psi <- matrix(, nsite, nyear) # occupancy probability
z <- matrix(, nsite, nyear) # occupancy status
psi[,1] <- inv.logit(cbind(1,x[,1,]) %*% beta0) # initial occupancy probability
z[,1] <- rbinom(nsite, 1, psi[,1]) # initial occupancy status

lon <- cov$lon # longitude or easting
lat <- cov$lat # latitude or northing
d <- matrix(, nsite, nsite) # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1
a <- ifelse(d < 1, 1, 0) # adjacency matrix
diag(a) <- 0 # a site is not a "neighbor" to ifself

phi   <- matrix(, nsite, nyear-1) # probability of persistence
gamma <- matrix(, nsite, nyear-1) # probability of colonization
for (t in 2:nyear) {
  for (i in 1:nsite) {
    phi  [i,t-1] <- inv.logit(c(1,x[i,t,]) %*% beta_phi)
    gamma[i,t-1] <- 1 - (1 - zeta) ^ sum(a[i,] * z[,t-1]) * (1 - delta)
    psi[i,t] <- z[i,t-1] * phi[i,t-1] + (1 - z[i,t-1]) * gamma[i,t-1]
    z[i,t] <- rbinom(1, 1, psi[i,t])
  } # i
} # t

w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

p <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    p[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha)
  } # t
} # i

y <- array(, dim=c(nsite, nyear, nreps)) # detection/non-detection data
for (i in 1:nsite) {
  for (t in 1:nyear) {
    y[i,t,] <- rbinom(nreps, 1, z[i,t] * p[i,t,])
  } # t
} # i

#=======================
# Define MCMC algorithm
#=======================
adj_dynam_occu_mcmc <- function(y, x, w, a, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ysum <- apply(y, 1:2, sum)

  beta0_save <- matrix(0, ncovs + 1, nmcmc)
  beta_phi_save <- matrix(0, ncovs + 1, nmcmc)
  zeta_save <- rep(0, nmcmc)
  delta_save <- rep(0, nmcmc)
  alpha_save <- matrix(0, npcvs + 1, nmcmc)
  z_save <- array(0, dim=c(nsite, nyear, 2000))

  # Priors
  beta0_mean <- rep(0, ncovs + 1)
  beta0_sd <- 2
  beta_phi_mean <- rep(0, ncovs + 1)
  beta_phi_sd <- 2
  logit_zeta_mean <- 0
  logit_zeta_sd <- 2
  logit_delta_mean <- 0
  logit_delta_sd <- 2
  alpha_mean <- rep(0, npcvs + 1)
  alpha_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs + 1)
  beta_phi <- rep(0, ncovs + 1)
  zeta <- 0.5
  delta <- 0.5
  alpha <- rep(0, npcvs + 1)
  z <- ifelse(ysum==0, 0, 1)
  phi   <- matrix(0.5, nsite, nyear-1)
  gamma <- matrix(, nsite, nyear-1)
  for (t in 2:nyear) {
    for (i in 1:nsite) {
      gamma[i,t-1] <- 1 - (1 - zeta) ^ sum(a[i,] * z[,t-1]) * (1 - delta)
    } # i
  } # t
  p <- array(0.5, dim=c(nsite, nyear, nreps))

  # Tuning factors
  beta0_tune <- 0.3
  beta_phi_tune <- 0.05
  zeta_tune <- 0.05
  delta_tune <- 0.05
  alpha_tune <- 0.025

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample beta0
    beta0_star <- rnorm(ncovs + 1, beta0, beta0_tune)
    mh1 <- sum(dbinom(z[,1], 1, inv.logit(cbind(1,x[,1,]) %*% beta0_star), log=TRUE)) + 
           sum(dnorm(beta0_star, beta0_mean, beta0_sd, log=TRUE))
    mh2 <- sum(dbinom(z[,1], 1, inv.logit(cbind(1,x[,1,]) %*% beta0     ), log=TRUE)) + 
           sum(dnorm(beta0     , beta0_mean, beta0_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta0 <- beta0_star
    }

    ### Sample beta_phi
    beta_phi_star <- rnorm(ncovs + 1, beta_phi, beta_phi_tune)
    phi_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi_star[,t-1] <- inv.logit(cbind(1,x[,t,]) %*% beta_phi_star)
    } # t
    mh1 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi_star + (1 - z[,-nyear]) * gamma, log=TRUE)) + 
           sum(dnorm(beta_phi_star, beta_phi_mean, beta_phi_sd, log=TRUE))
    mh2 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi      + (1 - z[,-nyear]) * gamma, log=TRUE)) + 
           sum(dnorm(beta_phi     , beta_phi_mean, beta_phi_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_phi <- beta_phi_star
      phi <- phi_star
    }

    ### Sample zeta
    zeta_star <- inv.logit(rnorm(1, logit(zeta), zeta_tune))
    gamma <- gamma_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      for (i in 1:nsite) {
        gamma     [i,t-1] <- 1 - (1 - zeta     ) ^ sum(a[i,] * z[,t-1]) * (1 - delta)
        gamma_star[i,t-1] <- 1 - (1 - zeta_star) ^ sum(a[i,] * z[,t-1]) * (1 - delta)
      } # i
    } # t
    mh1 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi + (1 - z[,-nyear]) * gamma_star, log=TRUE)) + 
           dnorm(logit(zeta_star), logit_zeta_mean, logit_zeta_sd, log=TRUE)
    mh2 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi + (1 - z[,-nyear]) * gamma     , log=TRUE)) + 
           dnorm(logit(zeta     ), logit_zeta_mean, logit_zeta_sd, log=TRUE)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      zeta <- zeta_star
    }

    ### Sample delta
    delta_star <- inv.logit(rnorm(1, logit(delta), delta_tune))
    gamma <- gamma_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      for (i in 1:nsite) {
        gamma     [i,t-1] <- 1 - (1 - zeta) ^ sum(a[i,] * z[,t-1]) * (1 - delta     )
        gamma_star[i,t-1] <- 1 - (1 - zeta) ^ sum(a[i,] * z[,t-1]) * (1 - delta_star)
      } # i
    } # t
    mh1 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi + (1 - z[,-nyear]) * gamma_star, log=TRUE)) + 
           dnorm(logit(delta_star), logit_delta_mean, logit_delta_sd, log=TRUE)
    mh2 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi + (1 - z[,-nyear]) * gamma     , log=TRUE)) + 
           dnorm(logit(delta     ), logit_delta_mean, logit_delta_sd, log=TRUE)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      delta <- delta_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs + 1, alpha, alpha_tune)
    p_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        p_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_star)
      } # t
    } # i
    mh1t <- mh2t <- numeric(nyear)
    for (t in 1:nyear) {
      mh1t[t] <- sum(dbinom(y[which(z[,t]==1),t,], 1, p_star[which(z[,t]==1),t,], log=TRUE))
      mh2t[t] <- sum(dbinom(y[which(z[,t]==1),t,], 1, p     [which(z[,t]==1),t,], log=TRUE))
    } # t
    mh1 <- sum(mh1t) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(mh2t) + 
           sum(dnorm(alpha,      alpha_mean, alpha_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Sample z
    s_temp <- matrix(, nsite, nyear)
    s_temp[,1] <- inv.logit(cbind(1,x[,1,]) %*% beta0)
    for (t in 2:nyear) {
      s_temp[,t] <- z[,t-1] * phi[,t-1] + (1 - z[,t-1]) * gamma[,t-1] + 1e-10
    } # t
    num_temp <- s_temp * apply((1 - p), 1:2, prod)
    psi_temp <- num_temp / (num_temp + (1 - s_temp))
    z[which(ysum==0)] <- rbinom(length(which(ysum==0)), 1, psi_temp[which(ysum==0)])

    ### Save samples
    beta0_save[,k] <- beta0
    beta_phi_save[,k] <- beta_phi
    zeta_save[k] <- zeta
    delta_save[k] <- delta
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        z_save[,,k-nmcmc+2000] <- z
      }
    }
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_phi_save=beta_phi_save, zeta_save=zeta_save, delta_save=delta_save, 
       alpha_save=alpha_save, 
       z_save=z_save)
} # adj_dynam_occu_mcmc

#==========
# Run MCMC
#==========
library(foreach)    # for parallel computing
library(doParallel) # for parallel computing

numCores <- round(detectCores() / 2) # only use half of the cores for parallel computing
registerDoParallel(numCores) # setup parallel computing

nmcmc <- 50000 # number of iterations
chain <- 3     # number of chains

start_time <- Sys.time() # start time of computing
out <- foreach (i=1:chain, .packages='boot') %dopar% {
  adj_dynam_occu_mcmc(y=y, x=x, w=w, a=a, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='2.3.1 adj dynamic occupancy model_output.RData')

#==============
# Plot results
#==============
pdf(file='2.3.1.1 adj dynamic occupancy model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["meadow"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["meadow"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]))
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
  ymin <- beta0[i] - yint * 8
  ymax <- beta0[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i], cex.main=2, line=1.5)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi[i] - yint * 4
  ymax <- beta_phi[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$zeta_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- zeta - yint * 2
ymax <- zeta + yint * 2
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=zeta, col='grey16', lwd=1.5)
title(main=expression(zeta), cex.main=2, line=1.2)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$delta_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- delta - yint * 2
ymax <- delta + yint * 2
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=delta, col='grey16', lwd=1.5)
title(main=expression(delta), cex.main=2, line=1.2)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

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
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=1, outer=T)

dev.off()

pdf(file='2.3.1.2 adj dynamic occupancy model_colonization.pdf', width=10, height=8)

library(vioplot)    # for making violin plots

nb_size <- matrix(, nsite, nyear)
for (i in 1:nsite) {
  for (t in 1:nyear) {
    nb_size[i,t] <- sum(cov$size[,t] * a[i,])
  } # t
} # i
nb_large <- ifelse(nb_size > median(nb_size), 1, 0)

colon_small <- colon_large <- numeric(nyear-1)
for (t in 1:(nyear-1)) {
  colon_small[t] <- length(which(z[which(nb_large[,t]==0),t]==0 & z[which(nb_large[,t]==0),t+1]==1)) / length(which(nb_large[,t]==0))
  colon_large[t] <- length(which(z[which(nb_large[,t]==1),t]==0 & z[which(nb_large[,t]==1),t+1]==1)) / length(which(nb_large[,t]==1))
} # i

z_post <- array(, dim=c(nsite, nyear, 6000))
for (i in 1:chain) {
  z_post[,,2000*(i-1)+c(1:2000)] <- out[[i]]$z_save
} # i

colon_post_small <- colon_post_large <- matrix(, 6000, nyear-1)
for (t in 1:(nyear-1)) {
  colon_post_small[,t] <- colMeans(ifelse(z_post[which(nb_large[,t]==0),t,]==0 & z_post[which(nb_large[,t]==0),t+1,]==1, 1, 0))
  colon_post_large[,t] <- colMeans(ifelse(z_post[which(nb_large[,t]==1),t,]==0 & z_post[which(nb_large[,t]==1),t+1,]==1, 1, 0))
} # i

par(mfrow=c(1,1))
par(mar=c(5,6.5,1,1))

plot(1, xlim=c(1,nyear-1), ylim=c(0,0.6), type='n', axes=F, xlab='', ylab='')
vioplot(colon_post_small, col='lightcoral', rectCol=NA, lineCol=NA, border=NA, add=T)
lines(colon_small, type='o', pch=16, cex=0.8, lwd=1.2, col='royalblue')
vioplot(colon_post_large, col='firebrick', rectCol=NA, lineCol=NA, border=NA, add=T)
lines(colon_large, type='o', pch=16, cex=0.8, lwd=1.2, col='navy')
axis(1, at=seq(4,nyear,5), labels=seq(4,nyear,5)+1, cex.axis=1.5)
axis(2, at=seq(0, 0.6, 0.1), cex.axis=1.5, las=2)
title(xlab='Year', cex.lab=2.5, line=3.5)
title(ylab='Proportion of Colonized Sites', cex.lab=2.5, line=4)

points(x=1.5, y=0.54, pch=16, cex=1, col='royalblue')
lines(x=c(1, 2), y=c(0.54,0.54), lwd=1.2, col='royalblue')
vioplot(rnorm(1e4,0.54,0.005), at=4.5, col='lightcoral', rectCol=NA, lineCol=NA, border=NA, add=T)
text(x=2.1, y=0.54, labels='True', pos=4, cex=1.6)
text(x=5.1, y=0.54, labels='Estimated', pos=4, cex=1.6)
text(x=0.7, y=0.58, labels='Small neighboring area', pos=4, cex=2)

points(x=11.5, y=.54, pch=16, cex=1, col='navy')
lines(x=c(11,12), y=c(0.54,0.54), lwd=1.2, col='navy')
vioplot(rnorm(1e4,0.54,0.005), at=14.5, col='firebrick', rectCol=NA, lineCol=NA, border=NA, add=T)
text(x=12.1, y=0.54, labels='True', pos=4, cex=1.6)
text(x=15.1, y=0.54, labels='Estimated', pos=4, cex=1.6)
text(x=10.7, y=0.58, labels='Large neighboring area', pos=4, cex=2)

dev.off()


