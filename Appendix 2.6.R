#=================================================================================================
# Spatial dynamic occupancy model for interacting species (with distance-based movement)
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

#setwd('')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(4)

# Basic values
nsite <- 360 # number of sites
nyear <- 20  # number of years
nreps <- 5   # number of within-season replicates
nspec <- 2   # number of species
ncovs <- 2   # number of environmental covariates
npcvs <- 2   # number of observational covariates

### For following parameters, the first column is for the predator species, and the second column is for the prey species
beta0 <- matrix(c(-1, 0.4, -0.2, 1, 1, -0.5), ncovs+1, nspec)            # intercept and slopes for initial occupancy
beta_phi <- matrix(c(-1, 1, 0, 0, 0.3, -0.3, 0.8, -0.4), ncovs+2, nspec) # intercept and slopes for persistence
kappa <- matrix(c(3, 1, 1.5, 2.5), nspec, nspec)                         # decay parameter for colonization
alpha <- matrix(c(0.8, -0.6, 0.3, 0.6, -0.4, 0.2), npcvs+1, nspec)       # intercept and slopes for detection probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- array(, dim=c(nsite, nyear, ncovs)) # environmental covariates
x[,,1] <- x1
x[,,2] <- x2

# Simulate data
psi <- array(, dim=c(nsite, nyear, nspec)) # occupancy probability
z <- array(, dim=c(nsite, nyear, nspec)) # occupancy status
psi[,1,] <- inv.logit(cbind(1,x[,1,]) %*% beta0) # initial occupancy probability
z[,1,] <- rbinom(nsite*nspec, 1, psi[,1,]) # initial occupancy status

lon <- cov$lon # longitude or easting
lat <- cov$lat # latitude or northing
d <- matrix(, nsite, nsite) # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1
eta <- array(, dim=c(nsite, nsite, nspec, nspec)) # unstandardized colonization rate
theta <- array(, dim=c(nsite, nsite, nspec, nspec)) # standardized colonization rate
for (s1 in 1:nspec) {
  for (s2 in 1:nspec) {
    eta[,,s1,s2] <- exp(-1 * kappa[s1,s2] * d)
    theta[,,s1,s2] <- eta[,,s1,s2] / rowSums(eta[,,s1,s2])
  } # s2
} # s1

phi   <- array(, dim=c(nsite, nyear-1, nspec)) # probability of persistence
gamma <- array(, dim=c(nsite, nyear-1, nspec)) # probability of colonization
for (t in 2:nyear) {
  for (i in 1:nsite) {
    for (s in 1:nspec) {
      phi  [i,t-1,s] <- inv.logit(c(1 - z[i,t-1,-s], z[i,t-1,-s], x[i,t,]) %*% beta_phi[,s])
      gamma[i,t-1,s] <- sum(theta[i,,,s] %*% c(1 - z[i,t-1,-s], z[i,t-1,-s]) * z[,t-1,s])
      psi[i,t,s] <- z[i,t-1,s] * phi[i,t-1,s] + (1 - z[i,t-1,s]) * gamma[i,t-1,s]
      z[i,t,s] <- rbinom(1, 1, psi[i,t,s])
    } # s
  } # i
} # t

w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

p <- array(, dim=c(nsite, nyear, nreps, nspec)) # detection probability
for (t in 1:nyear) {
  for (j in 1:nreps) {
    p[,t,j,] <- inv.logit(cbind(1,w[,t,j,]) %*% alpha)
  } # j
} # t

y <- array(, dim=c(nsite, nyear, nreps, nspec)) # detection/non-detection data
for (t in 1:nyear) {
  for (j in 1:nreps) {
    for (s in 1:nspec) {
      y[,t,j,s] <- rbinom(nsite, 1, z[,t,s] * p[,t,j,s])
    } # s
  } # j
} # t

#=======================
# Define MCMC algorithm
#=======================
sdom_ii_mcmc <- function(y, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  nspec <- dim(y)[4]
  ncovs <- dim(x)[3] 
  npcvs <- dim(w)[4] 
  ysum <- apply(y, c(1,2,4), sum)

  beta0_save <- array(0, dim=c(ncovs + 1, nspec, nmcmc))
  beta_phi_save <- array(0, dim=c(ncovs + 2, nspec, nmcmc))
  kappa_save <- array(0, dim=c(nspec, nspec, nmcmc))
  alpha_save <- array(0, dim=c(npcvs + 1, nspec, nmcmc))

  # Priors
  beta0_mean <- matrix(0, ncovs + 1, nspec)
  beta0_sd <- 2
  beta_phi_mean <- matrix(0, ncovs + 2, nspec)
  beta_phi_sd <- 2
  log_kappa_mean <- matrix(0, nspec, nspec)
  log_kappa_sd <- 2
  alpha_mean <- matrix(0, npcvs + 1, nspec)
  alpha_sd <- 2

  # Starting values
  beta0 <- matrix(0, ncovs + 1, nspec)
  beta_phi <- matrix(0, ncovs + 2, nspec)
  kappa <- matrix(1, nspec,nspec)
  alpha <- matrix(0, npcvs + 1, nspec)
  phi   <- array(0.5, dim=c(nsite, nyear-1, nspec))
  gamma <- array(0.5, dim=c(nsite, nyear-1, nspec))
  eta <- theta <- array(, dim=c(nsite, nsite, nspec, nspec))
  for (s1 in 1:nspec) {
    for (s2 in 1:nspec) {
      eta[,,s1,s2] <- exp(-1 * kappa[s1,s2] * d)
      theta[,,s1,s2] <- eta[,,s1,s2] / rowSums(eta[,,s1,s2])
    } # s2
  } # s1
  p <- array(0.5, dim=c(nsite, nyear, nreps, nspec))
  z <- ifelse(ysum==0, 0, 1)

  # Tuning factors
  beta0_tune <- matrix(0.4, ncovs + 1, nspec)
  beta_phi_tune <- matrix(0.05, ncovs + 2, nspec)
  kappa_tune <- matrix(0.1, nspec, nspec)
  alpha_tune <- matrix(0.05, npcvs + 1, nspec)

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample beta0
    beta0_star <- matrix(rnorm((ncovs + 1) * nspec, beta0, beta0_tune), ncovs + 1, nspec)
    for (s in 1:nspec) {
      mh1 <- sum(dbinom(z[,1,s], 1, inv.logit(cbind(1,x[,1,]) %*% beta0_star[,s]), log=TRUE)) + 
             sum(dnorm(beta0_star[,s], beta0_mean[,s], beta0_sd, log=TRUE))
      mh2 <- sum(dbinom(z[,1,s], 1, inv.logit(cbind(1,x[,1,]) %*% beta0     [,s]), log=TRUE)) + 
             sum(dnorm(beta0     [,s], beta0_mean[,s], beta0_sd, log=TRUE))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta0[,s] <- beta0_star[,s]
      }
    } # s

    ### Sample beta_phi
    beta_phi_star <- matrix(rnorm((ncovs + 2) * nspec, beta_phi, beta_phi_tune), ncovs + 2, nspec)
    phi <- phi_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (t in 2:nyear) {
      for (s in 1:nspec) {
        phi     [,t-1,s] <- inv.logit(cbind(1 - z[,t-1,-s], z[,t-1,-s], x[,t,]) %*% beta_phi     [,s])
        phi_star[,t-1,s] <- inv.logit(cbind(1 - z[,t-1,-s], z[,t-1,-s], x[,t,]) %*% beta_phi_star[,s])
      } # k
    } # t
    for (s in 1:nspec) {
      mh1 <- sum(dbinom(z[,-1,s], 1, z[,-nyear,s] * phi_star[,,s] + (1 - z[,-nyear,s]) * gamma[,,s], log=TRUE)) + 
             sum(dnorm(beta_phi_star[,s], beta_phi_mean[,s], beta_phi_sd, log=TRUE))
      mh2 <- sum(dbinom(z[,-1,s], 1, z[,-nyear,s] * phi     [,,s] + (1 - z[,-nyear,s]) * gamma[,,s], log=TRUE)) + 
             sum(dnorm(beta_phi     [,s], beta_phi_mean[,s], beta_phi_sd, log=TRUE))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta_phi[,s] <- beta_phi_star[,s]
      }
    } # s

    ### Sample kappa
    kappa_star <- matrix(exp(rnorm(nspec * nspec, log(kappa), kappa_tune)), nspec, nspec)
    eta_star <- theta_star <- array(, dim=c(nsite, nsite, nspec, nspec))
    for (s1 in 1:nspec) {
      for (s2 in 1:nspec) {
        eta_star[,,s1,s2] <- exp(-1 * kappa_star[s1,s2] * d)
        theta_star[,,s1,s2] <- eta_star[,,s1,s2] / rowSums(eta_star[,,s1,s2])
      } # s2
    } # s1
    gamma <- gamma_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (t in 2:nyear) {
      for (i in 1:nsite) {
        for (s in 1:nspec) {
          gamma     [i,t-1,s] <- sum(theta     [i,,,s] %*% c(1 - z[i,t-1,-s], z[i,t-1,-s]) * z[,t-1,s])
          gamma_star[i,t-1,s] <- sum(theta_star[i,,,s] %*% c(1 - z[i,t-1,-s], z[i,t-1,-s]) * z[,t-1,s])
        } # s
      } # i
    } # t
    for (s in 1:nspec) {
      mh1 <- sum(dbinom(z[,-1,s], 1, z[,-nyear,s] * phi[,,s] + (1 - z[,-nyear,s]) * gamma_star[,,s], log=TRUE)) + 
             sum(dnorm(log(kappa_star[,s]), log_kappa_mean[,s], log_kappa_sd, log=TRUE))
      mh2 <- sum(dbinom(z[,-1,s], 1, z[,-nyear,s] * phi[,,s] + (1 - z[,-nyear,s]) * gamma     [,,s], log=TRUE)) + 
             sum(dnorm(log(kappa     [,s]), log_kappa_mean[,s], log_kappa_sd, log=TRUE))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        kappa[,s] <- kappa_star[,s]
        eta[,,,s] <- eta_star[,,,s]
        theta[,,,s] <- theta_star[,,,s]
      }
    } # s

    ### Sample alpha
    alpha_star <- matrix(rnorm((npcvs + 1) * nspec, alpha, alpha_tune), npcvs + 1, nspec)
    p_star <- array(, dim=c(nsite, nyear, nreps, nspec))
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        p_star[,t,j,] <- inv.logit(cbind(1,w[,t,j,]) %*% alpha_star)
      } # j
    } # t
    for (s in 1:nspec) {
      mh1t <- mh2t <- numeric(nyear)
      for (t in 1:nyear) {
        mh1t[t] <- sum(dbinom(y[which(z[,t,s]==1),t,,s], 1, p_star[which(z[,t,s]==1),t,,s], log=TRUE))
        mh2t[t] <- sum(dbinom(y[which(z[,t,s]==1),t,,s], 1, p     [which(z[,t,s]==1),t,,s], log=TRUE))
      } # t
      mh1 <- sum(mh1t) + 
             sum(dnorm(alpha_star[,s], alpha_mean[,s], alpha_sd, log=TRUE))
      mh2 <- sum(mh2t) + 
             sum(dnorm(alpha     [,s], alpha_mean[,s], alpha_sd, log=TRUE))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        alpha[,s] <- alpha_star[,s]
        p[,,,s] <- p_star[,,,s]
      }
    } # s

    ### Sample z
    s_temp <- array(, dim=c(nsite, nyear, nspec))
    for (s in 1:nspec) {
      s_temp[,1,s] <- inv.logit(cbind(1,x[,1,]) %*% beta0[,s])
      for (t in 2:nyear) {
        s_temp[,t,s] <- z[,t-1,s] * phi[,t-1,s] + (1 - z[,t-1,s]) * gamma[,t-1,s]
      } # t
    } # s
    num_temp <- s_temp * apply((1 - p), c(1,2,4), prod)
    psi_temp <- num_temp / (num_temp + (1 - s_temp))
    z[which(ysum==0)] <- rbinom(length(which(ysum==0)), 1, psi_temp[which(ysum==0)])

    ### Save samples
    beta0_save[,,k] <- beta0
    beta_phi_save[,,k] <- beta_phi
    kappa_save[,,k] <- kappa
    alpha_save[,,k] <- alpha
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_phi_save=beta_phi_save, kappa_save=kappa_save, 
       alpha_save=alpha_save)
} # sdom_ii_mcmc

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
  sdom_ii_mcmc(y=y, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='2.4 sdom_ii_output.RData')

#==============
# Plot results
#==============
pdf(file='2.4.1 sdom_ii_chains.pdf', width=10, height=10)

library(rstan)

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- matrix(
  c(expression(beta["0,ant"]^"[0]"), expression(beta["meadow,ant"]^"[0]"), expression(beta["temperature,ant"]^"[0]"), 
    expression(beta["0,butterfly"]^"[0]"), expression(beta["meadow,butterfly"]^"[0]"), expression(beta["temperature,butterfly"]^"[0]")), 
  ncovs+1, nspec)
ylab_beta_phi <- matrix(
  c(expression(beta["0,ant,z"[butterfly]==0]^""["["*phi*"]"]), expression(beta["0,ant,z"[butterfly]==1]^""["["*phi*"]"]), 
    expression(beta["meadow,ant"]^""["["*phi*"]"]), expression(beta["temperature,ant"]^""["["*phi*"]"]), 
    expression(beta["0,butterfly,z"[ant]==0]^""["["*phi*"]"]), expression(beta["0,butterfly,z"[ant]==1]^""["["*phi*"]"]), 
    expression(beta["meadow,butterfly"]^""["["*phi*"]"]), expression(beta["temperature,butterfly"]^""["["*phi*"]"])), 
  ncovs+2, nspec)
ylab_kappa <- matrix(
  c(expression(kappa["ant,z"[butterfly]==0]), expression(kappa["ant,z"[butterfly]==1]), 
    expression(kappa["butterfly,z"[ant]==0]), expression(kappa["butterfly,z"[ant]==1])), 
  nspec, nspec)
ylab_alpha <- matrix(
  c(expression(alpha["0,ant"]), expression(alpha["1,ant"]), expression(alpha["2,ant"]), 
    expression(alpha["0,butterfly"]), expression(alpha["1,butterfly"]), expression(alpha["2,butterfly"])), 
    npcvs+1, nspec)

par(mfrow=c(8,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,1))

for (s in 1:nspec) {
for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta0_save[i,s,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0[i,s] - yint * 8
  ymax <- beta0[i,s] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0[i,s], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i,s], cex.main=2, line=1.4)
  text(x=nmcmc*0.42, y=beta0[i,s]+yint*6, pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i
} # s

for (s in 1:nspec) {
for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_save[i,s,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ad <- matrix(c(8,8,8,8,4,8,4,4), ncovs+2, nspec)
  ymin <- beta_phi[i,s] - yint * ad[i,s]
  ymax <- beta_phi[i,s] + yint * ad[i,s]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*ad[i,s]/2), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi[i,s], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi[i,s], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_phi[i,s]+yint*ad[i,s]*.75, pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i
} # s

for (s1 in 1:nspec) {
for (s2 in 1:nspec) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$kappa_save[s2,s1,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- kappa[s2,s1] - yint * 8
  ymax <- kappa[s2,s1] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=kappa[s2,s1], col='grey16', lwd=1.5)
  title(main=ylab_kappa[s2,s1], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=kappa[s2,s1]+yint*6, pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # s2
} # s1

for (s in 1:nspec) {
for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_save[i,s,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i,s] - yint * c(4,2)[s]
  ymax <- alpha[i,s] + yint * c(4,2)[s]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  if (s == 1) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*c(2,1)[s]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i,s], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i,s], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=alpha[i,s]+yint*c(3,1.5)[s], pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i
} # s

title(xlab='Iteration', cex.lab=3, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=3, line=0.8, outer=T)

dev.off()

pdf(file='2.4.2 sdom_ii_eta.pdf', width=9, height=8)

kappa_post <- array(, dim=c(nspec, nspec, nmcmc/2*chain))
for (s1 in 1:nspec) {
  for (s2 in 1:nspec) {
    for (i in 1:chain) {
      if (i == 1) {
        tt <- out[[i]]$kappa_save[s1,s2,(nmcmc/2+1):nmcmc]
      } else {
        tt <- c(tt, out[[i]]$kappa_save[s1,s2,(nmcmc/2+1):nmcmc])
      }
    } # i
    kappa_post[s1,s2,] <- tt
  } # s2
} # s1
kappa_med <- apply(kappa_post, 1:2, median)

eta_true <- eta_pred <- array(, dim=c(nsite, nsite, nspec, nspec))
for (s1 in 1:nspec) {
  for (s2 in 1:nspec) {
    eta_true[,,s1,s2] <- exp(-1 * kappa    [s1,s2] * d)
    eta_pred[,,s1,s2] <- exp(-1 * kappa_med[s1,s2] * d)
  } # s2
} # s1

min_eta <- 0.2
max_lwd <- 3
titles <- c(expression(z["ant"]==0), expression(z["ant"]==1))

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(4,5,5,2))

for (s in 1:nspec) {

  plot(cov$lon, cov$lat, xlim=c(-10.2,10.2), ylim=c(-8.2,8.2), axes=F, xlab='', ylab='', 
       pch=21, cex=rowMeans(cov$size)*4.6, col='grey18', bg='white')
  axis(1, at=seq(-10,10,5), labels=rep('',5))
  if (s == 1) {
    axis(2, at=seq(-8,8,4), cex.axis=1.6, las=2)
  } else {
    axis(2, at=seq(-8,8,4), labels=rep('',5))
    axis(4, at=0, labels='True', cex.axis=2.5, tick=F, line=.5)
  }
  box()
  for (i in 1:(nsite-1)) {
    for (j in (i+1):nsite) {
      if (eta_true[i,j,s,2] >= min_eta) {
        lines(x=cov$lon[c(i,j)], y=cov$lat[c(i,j)], col='royalblue', lwd=eta_true[i,j,s,2]*max_lwd)
      }
    } # j
  } # i
  axis(3, at=0, labels=titles[s], cex.axis=2.5, tick=F, line=-.8)

} # s

for (s in 1:nspec) {

  plot(cov$lon, cov$lat, xlim=c(-10.2,10.2), ylim=c(-8.2,8.2), axes=F, xlab='', ylab='', 
       pch=21, cex=rowMeans(cov$size)*4.6, col='grey18', bg='white')
  axis(1, at=seq(-10,10,5), cex.axis=1.6, las=1)
  if (s == 1) {
    axis(2, at=seq(-8,8,4), cex.axis=1.6, las=2)
  } else {
    axis(2, at=seq(-8,8,4), labels=rep('',5))
    axis(4, at=0, labels='Estimated', cex.axis=2.5, tick=F, line=.5)
  }
  box()
  for (i in 1:(nsite-1)) {
    for (j in (i+1):nsite) {
      if (eta_pred[i,j,s,2] >= min_eta) {
        lines(x=cov$lon[c(i,j)], y=cov$lat[c(i,j)], col='lightcoral', lwd=eta_pred[i,j,s,2]*max_lwd)
      }
    } # j
  } # i

} # s

title(xlab='Easting' , line=2.6, cex.lab=2.8, outer=T)
title(ylab='Northing', line=2.4, cex.lab=2.8, outer=T)
title(main=expression(paste('Colonization Rate (', eta, ') of Butterfly', sep='')), line=3.2, cex.main=3, outer=T)

dev.off()


