#=================================================================================================
# Spatial dynamic N-mixture model with population growth for interacting species
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
nsite <- 360  # number of sites
nyear <- 20   # number of years
nreps <- 5    # number of within-season replicates
nspec <- 2    # number of species
ncovs <- 2    # number of environmental covariates
npcvs <- 2    # number of observational covariates

### For beta0 and beta_lambda, the first column is for the invasive red-eared slider, and the second column is for the native Blanding's turtle
beta0 <- matrix(c(3, 0.3, 0.3, 1.6, 0.6, -0.6), ncovs+1, nspec)                             # intercept and slopes for initial population size
beta_rho <- matrix(c(0, -0.4, 0, 0.1, 0.1, 0, -0.8, -0.2, 0.2, -0.2), nspec+ncovs+1, nspec) # intercept and slopes for local population growth rate
kappa <- c(3, 1)                                                                            # distance effect on movement
alpha <- c(0.4, -0.6, 0.3)                                                                  # intercept and slopes for detection probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log food patch size
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
lambda <- array(, dim=c(nsite, nyear, nspec)) # expectation of abundance
N <- array(, dim=c(nsite, nyear, nspec)) # abundance
for (s in 1:nspec) {
  lambda[,1,s] <- exp(cbind(1,x[,1,]) %*% beta0[,s]) # expectation of initial abundance
  N[,1,s] <- rpois(nsite, lambda[,1,s]) # initial abundance
} # s

lon <- cov$lon # longitude or easting
lat <- cov$lat # latitude or northing
d <- matrix(, nsite, nsite) # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1

eta   <- array(, dim=c(nsite, nsite, nspec)) # unstandardized colonization rate
theta <- array(, dim=c(nsite, nsite, nspec)) #   standardized colonization rate
for (s in 1:nspec) {
  eta[,,s] <- exp(-1 * kappa[s] * d)
  theta[,,s] <- eta[,,s] / rowSums(eta[,,s])
} # s

rho <- array(, dim=c(nsite, nyear-1, nspec)) # population growth rate
for (t in 2:nyear) {
  for (s in 1:nspec) {
    rho[,t-1,s] <- exp(cbind(1, (N[,t-1,]-lambda[,1,])/lambda[,1,], x[,t,]) %*% beta_rho[,s])
    lambda[,t,s] <- (N[,t-1,s] * rho[,t-1,s]) %*% theta[,,s]
    N[,t,s] <- rpois(nsite, lambda[,t,s])
  } # s
} # t

w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

p <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (t in 1:nyear) {
  for (j in 1:nreps) {
    p[,t,j] <- inv.logit(cbind(1,w[,t,j,]) %*% alpha)
  } # j
} # t

y <- array(, dim=c(nsite, nyear, nreps, nspec)) # count data
for (t in 1:nyear) {
  for (j in 1:nreps) {
    for (s in 1:nspec) {
      y[,t,j,s] <- rbinom(nsite, N[,t,s], p[,t,j])
    } # s
  } # j
} # t

#=======================
# Define MCMC algorithm
#=======================
spatial_growth_Nmix_ii_mcmc <- function(y, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  nspec <- dim(y)[4]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, c(1,2,4), max)

  beta0_save <- array(0, dim=c(ncovs+1, nspec, nmcmc))
  beta_rho_save <- array(0, dim=c(nspec+ncovs+1, nspec, nmcmc))
  kappa_save <- matrix(0, nspec, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc)
  N_save <- array(0, dim=c(nsite, nyear, nspec, 2000))

  # Priors
  beta0_mean <- matrix(0, ncovs+1, nspec)
  beta0_sd <- 2
  beta_rho_mean <- matrix(0, nspec+ncovs+1, nspec)
  beta_rho_sd <- 2
  log_kappa_mean <- rep(0, nspec)
  log_kappa_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta0 <- matrix(0, ncovs+1, nspec)
  beta_rho <- matrix(0, nspec+ncovs+1, nspec)
  kappa <- rep(1, nspec)
  alpha <- rep(0, npcvs+1)
  eta <- theta <- array(, dim=c(nsite, nsite, nspec))
  for (s in 1:nspec) {
    eta[,,s] <- exp(-1 * kappa[s] * d)
    theta[,,s] <- eta[,,s] / rowSums(eta[,,s])
  } # s
  rho <- array(1, dim=c(nsite, nyear-1, nspec))
  p <- array(0.5, dim=c(nsite, nyear, nreps))
  N <- round((ymax + 1) / 0.5)

  # Tuning factors
  beta0_tune <- matrix(0.035, ncovs+1, nspec)
  beta_rho_tune <- matrix(c(0.05, 0.07), nspec+ncovs+1, nspec, byrow=T)
  kappa_tune <- 0.15
  alpha_tune <- c(0.02, 0.02, 0.02)
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- array(rpois(nsite * nyear * nspec, N + N_tune), dim=c(nsite, nyear, nspec))
    lambda0 <- matrix(, nsite, nspec) # lambda for the first year
    lambdar <- array(, dim=c(nsite, nyear-1, nspec)) # lambda except for the first year
    mh1 <- mh2 <- array(, dim=c(nsite, nyear, nspec))
    for (s in 1:nspec) {
      lambda0[,s] <- exp(cbind(1, x[,1,]) %*% beta0[,s])
      for (t in 2:nyear) {
        lambdar[,t-1,s] <- (N[,t-1,s] * rho[,t-1,s]) %*% theta[,,s]
      } # t
      mh1[,,s] <- apply(dbinom(y[,,,s], N_star[,,s], p, log=T), 1:2, sum) + 
                  dpois(N_star[,,s], cbind(lambda0[,s], lambdar[,,s]), log=T) + 
                  dpois(N[,,s], N_star[,,s] + N_tune, log=T)
      mh2[,,s] <- apply(dbinom(y[,,,s], N     [,,s], p, log=T), 1:2, sum) + 
                  dpois(N     [,,s], cbind(lambda0[,s], lambdar[,,s]), log=T) + 
                  dpois(N_star[,,s], N[,,s] + N_tune, log=T)
    } # s
    mh <- exp(mh1 - mh2)
    Nkeep <- ((mh > runif(nsite*nyear*nspec)) & (N_star >= ymax))
    N[Nkeep] <- N_star[Nkeep]

    ### Sample beta0
    beta0_star <- matrix(rnorm((ncovs+1)*nspec, beta0, beta0_tune), ncovs+1, nspec)
    lambda0_star <- matrix(, nsite, nspec)
    for (s in 1:nspec) {
      lambda0_star[,s] <- exp(cbind(1, x[,1,]) %*% beta0_star[,s])
      mh1 <- sum(dpois(N[,1,s], exp(cbind(1, x[,1,]) %*% beta0_star[,s]), log=T)) + 
             sum(dnorm(beta0_star[,s], beta0_mean[,s], beta0_sd, log=T))
      mh2 <- sum(dpois(N[,1,s], exp(cbind(1, x[,1,]) %*% beta0     [,s]), log=T)) + 
             sum(dnorm(beta0     [,s], beta0_mean[,s], beta0_sd, log=T))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta0[,s] <- beta0_star[,s]
        lambda0[,s] <- lambda0_star[,s]
      }
    } # s

    ### Sample beta_rho
    beta_rho_star <- matrix(rnorm((nspec+ncovs+1)*nspec, beta_rho, beta_rho_tune), nspec+ncovs+1, nspec)
    rho_star <- lambdar <- lambdar_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (s in 1:nspec) {
      for (t in 2:nyear) {
        rho_star[,t-1,s] <- exp(cbind(1,(N[,t-1,]-lambda0)/lambda0,x[,t,]) %*% beta_rho_star[,s])
      } # t
      lambdar     [,,s] <- t(t(N[,-nyear,s] * rho     [,,s]) %*% theta[,,s])
      lambdar_star[,,s] <- t(t(N[,-nyear,s] * rho_star[,,s]) %*% theta[,,s])
      mh1 <- sum(dpois(N[,-1,s], lambdar_star[,,s], log=T)) + 
             sum(dnorm(beta_rho_star[,s], beta_rho_mean[,s], beta_rho_sd, log=T))
      mh2 <- sum(dpois(N[,-1,s], lambdar     [,,s], log=T)) + 
             sum(dnorm(beta_rho     [,s], beta_rho_mean[,s], beta_rho_sd, log=T))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        beta_rho[,s] <- beta_rho_star[,s]
        rho[,,s] <- rho_star[,,s]
      }
    } # s

    ### Sample kappa
    kappa_star <- exp(rnorm(nspec, log(kappa), kappa_tune))
    eta_star <- theta_star <- array(, dim=c(nsite, nsite, nspec))
    lambdar <- lambdar_star <- array(, dim=c(nsite, nyear-1, nspec))
    for (s in 1:nspec) {
      eta_star[,,s] <- exp(-1 * kappa_star[s] * d)
      theta_star[,,s] <- eta_star[,,s] / rowSums(eta_star[,,s])
      lambdar     [,,s] <- t(t(N[,-nyear,s] * rho[,,s]) %*% theta     [,,s])
      lambdar_star[,,s] <- t(t(N[,-nyear,s] * rho[,,s]) %*% theta_star[,,s])
      mh1 <- sum(dpois(N[,-1,s], lambdar_star[,,s], log=T)) + 
             sum(dnorm(log(kappa_star[s]), log_kappa_mean[s], log_kappa_sd, log=T))
      mh2 <- sum(dpois(N[,-1,s], lambdar     [,,s], log=T)) + 
             sum(dnorm(log(kappa     [s]), log_kappa_mean[s], log_kappa_sd, log=T))
      mh <- exp(mh1 - mh2)
      if (mh > runif(1)) {
        kappa[s] <- kappa_star[s]
        eta[,,s] <- eta_star[,,s]
        theta[,,s] <- theta_star[,,s]
      }
    } # s

    ### Sample alpha
    alpha_star <- rnorm(npcvs+1, alpha, alpha_tune)
    p_star <- array(0, dim=c(nsite, nyear, nreps))
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        p_star[,t,j] <- inv.logit(cbind(1,w[,t,j,]) %*% alpha_star)
      } # j
    } # t
    mh1t <- mh2t <- rep(0, nspec)
    for (s in 1:nspec) {
      mh1t[s] <- sum(dbinom(y[,,,s], N[,,s], p_star, log=T))
      mh2t[s] <- sum(dbinom(y[,,,s], N[,,s], p     , log=T))
    } # s
    mh1 <- sum(mh1t) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=T))
    mh2 <- sum(mh2t) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Save samples
    beta0_save[,,k] <- beta0
    beta_rho_save[,,k] <- beta_rho
    kappa_save[,k] <- kappa
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        N_save[,,,k-nmcmc+2000] <- N
      }
    }
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_rho_save=beta_rho_save, kappa_save=kappa_save, 
       alpha_save=alpha_save, 
       N_save=N_save)

} # spatial_growth_Nmix_ii_mcmc

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
  spatial_growth_Nmix_ii_mcmc(y=y, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.4 spatial population growth N-mixture model for interacting species_output.RData')

#==============
# Plot results
#==============
pdf(file='3.4.1 spatial population growth N-mixture model for interacting species_chains.pdf', width=10, height=10)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- matrix(
  c(expression(beta["0,red-eared"]^"[0]"), expression(beta["pond,red-eared"]^"[0]"), expression(beta["temperature,red-eared"]^"[0]"), 
    expression(beta["0,Blanding's"]^"[0]"), expression(beta["pond,Blanding's"]^"[0]"), expression(beta["temperature,Blanding's"]^"[0]")), 
  ncovs+1, nspec)
ylab_beta_rho <- matrix(
  c(expression(beta["0,red-eared"]^""["["*rho*"]"]), expression(beta["density,red-eared"]^""["["*rho*"]"]), 
    expression(beta["competition,red-eared"]^""["["*rho*"]"]), expression(beta["pond,red-eared"]^""["["*rho*"]"]), 
    expression(beta["temperature,red-eared"]^""["["*rho*"]"]), 
    expression(beta["0,Blanding's"]^""["["*rho*"]"]), expression(beta["competition,Blanding's"]^""["["*rho*"]"]), 
    expression(beta["density,Blanding's"]^""["["*rho*"]"]), expression(beta["pond,Blanding's"]^""["["*rho*"]"]), 
    expression(beta["temperature,Blanding's"]^""["["*rho*"]"])), 
  nspec+ncovs+1, nspec)
ylab_kappa <- c(expression(kappa["red-eared"]), expression(kappa["Blanding's"]))
ylab_alpha <- c(expression(alpha["0"]), expression(alpha["1"]), expression(alpha["2"]))

par(mfrow=c(7,3))
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
  ymin <- beta0[i,s] - yint * matrix(c(2,2,2,4,4,4),ncovs+1,nspec)[i,s]
  ymax <- beta0[i,s] + yint * matrix(c(2,2,2,4,4,4),ncovs+1,nspec)[i,s]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*matrix(c(1,1,1,2,2,2),ncovs+1,nspec)[i,s]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0[i,s], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i,s], cex.main=2, line=1.4)
  text(x=nmcmc*0.42, y=beta0[i,s]+yint*matrix(c(1.5,1.5,1.5,3,3,3),ncovs+1,nspec)[i,s], pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i
} # s

for (s in 1:nspec) {
for (i in 1:(ncovs+nspec+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_rho_save[i,s,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_rho[i,s] - yint * matrix(c(2,4,2,2,2,2,4,4,4,2),ncovs+nspec+1,nspec)[i,s]
  ymax <- beta_rho[i,s] + yint * matrix(c(2,4,2,2,2,2,4,4,4,2),ncovs+nspec+1,nspec)[i,s]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*matrix(c(1,2,1,1,1,1,2,2,2,1),ncovs+nspec+1,nspec)[i,s]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_rho[i,s], col='grey16', lwd=1.5)
  title(main=ylab_beta_rho[i,s], cex.main=2, line=1.1)
  text(x=nmcmc*0.42, y=beta_rho[i,s]+yint*matrix(c(1.5,3,1.5,1.5,1.5,1.5,3,3,3,1.5),ncovs+nspec+1,nspec)[i,s], pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i
} # s

for (s in 1:nspec) {
tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$kappa_save[s,]
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa[s] - yint * c(8,4)[s]
ymax <- kappa[s] + yint * c(8,4)[s]
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*c(4,2)[s]), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa[s], col='grey16', lwd=1.5)
title(main=ylab_kappa[s], cex.main=2, line=1.1)
text(x=nmcmc*0.42, y=kappa[s]+yint*c(6,3)[s], pos=4, cex=1.2, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # s

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
  if (s == 1) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1)
  text(x=nmcmc*0.42, y=alpha[i]+yint*1.5, pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=3, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=3, line=0.8, outer=T)

dev.off()

pdf(file='3.4.2 spatial population growth N-mixture model for interacting species_Nsum.pdf', width=10, height=8)

library(vioplot) # for making violin plots

N_post <- array(, dim=c(nsite, nyear, nspec, 2000*chain))
for (i in 1:chain) {
  N_post[,,,(i-1)*2000+1:2000] <- out[[i]]$N_save
} # i
Nsum <- apply(N, 2:3, sum)
Nsum_post <- apply(N_post, 2:4, sum)

par(mfrow=c(1,1))
par(mar=c(5,7.5,1,7))

plot(1, xlim=c(1,nyear), ylim=c(0, 16000), type='n', axes=F, xlab='', ylab='')
vioplot(t(Nsum_post[,1,]*1), col='firebrick' , rectCol=NA, lineCol=NA, border=NA, add=T)
vioplot(t(Nsum_post[,2,]*5), col='lightcoral', rectCol=NA, lineCol=NA, border=NA, add=T)
lines(Nsum[,1]*1, type='o', pch=16, cex=.8, lwd=1.2, col='navy'     )
lines(Nsum[,2]*5, type='o', pch=16, cex=.8, lwd=1.2, col='royalblue')
axis(1, at=seq(0,nyear,5), cex.axis=1.5)
axis(2, at=seq(0,16000,4000), cex.axis=1.5, las=2)
axis(4, at=seq(0,16000,4000), labels=seq(0,16000,4000)/5, cex.axis=1.5, las=2)
box()
title(xlab='Year', cex.lab=2.5, line=3.5)
axis(2, at=8000, labels="Total Population Size of Red-eared" , cex.axis=2.5, line=4.4, tick=F)
axis(4, at=8000, labels="Total Population Size of Blanding's", cex.axis=2.5, line=4.6, tick=F)

text(x= 5.5, y=16000, labels="Red-eared", pos=1, cex=2.5)
text(x=15.5, y=16000, labels="Blanding's", pos=1, cex=2.5)
points(x=2.5, y=14000, pch=16, cex=1, col='navy')
lines(x=c(2,3), y=c(14000,14000), lwd=1.2, col='navy')
text(x=3, y=14000, labels='True', pos=4, cex=1.5)
vioplot(rnorm(1e4,14000,100), at=5.5, col='firebrick', rectCol=NA, lineCol=NA, border=NA, add=T)
text(x=6, y=14000, labels='Estimated', pos=4, cex=1.5)
points(x=12.5, y=14000, pch=16, cex=1, col='royalblue')
lines(x=c(12,13), y=c(14000,14000), lwd=1.2, col='royalblue')
text(x=13, y=14000, labels='True', pos=4, cex=1.5)
vioplot(rnorm(1e4,14000,100), at=15.5, col='lightcoral', rectCol=NA, lineCol=NA, border=NA, add=T)
text(x=16, y=14000, labels='Estimated', pos=4, cex=1.5)

dev.off()

pdf(file='3.4.3 spatial population growth N-mixture model for interacting species_Nprop.pdf', width=9, height=8)

library(plotrix) # for making pie charts at selected positions

years_plot <- c(3, 18)
N_prop <- N[,years_plot,1] / apply(N[,years_plot,], 1:2, sum)
N_post_med <- apply(N_post, 1:3, median)
N_prop_post <- N_post_med[,years_plot,1] / apply(N_post_med[,years_plot,], 1:2, sum)
size <- rowMeans(cov$size)

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(4,5,5,2))

for (t in 1:2) {
  plot(1, xlim=c(-10, 10), ylim=c(-8, 8), type='n', axes=F, xlab='', ylab='')
  for (i in 1:nsite) {
    floating.pie(xpos=cov$lon[i], ypos=cov$lat[i], x=c(N_prop[i,t],1-N_prop[i,t]), 
                 radius=size[i]*0.9, col=c('royalblue','white'), border='grey18')
  } # i
  axis(1, at=seq(-10,10,5), labels=rep('',5))
  if (t == 1) {
    axis(2, at=seq(-8,8,4), las=2, cex.axis=1.6)
  } else {
    axis(2, at=seq(-8,8,4), labels=rep('',5))
  }
  axis(3, at=0, labels=paste('Year', years_plot[t], sep=' '), cex.axis=2.5, line=-.8, tick=F)
  if (t == 2) {
    axis(4, at=0, labels='True', cex.axis=2.5, line=0.4, tick=F)
  }
  box()
} # t

for (t in 1:2) {
  plot(1, xlim=c(-10, 10), ylim=c(-8, 8), type='n', axes=F, xlab='', ylab='')
  for (i in 1:nsite) {
    floating.pie(xpos=cov$lon[i], ypos=cov$lat[i], x=c(N_prop_post[i,t],1-N_prop_post[i,t]), 
                 radius=size[i]*0.9, col=c('lightcoral','white'), border='grey18')
  } # i
  axis(1, at=seq(-10,10,5), cex.axis=1.6)
  if (t == 1) {
    axis(2, at=seq(-8,8,4), las=2, cex.axis=1.6)
  } else {
    axis(2, at=seq(-8,8,4), labels=rep('',5))
  }
  if (t == 2) {
    axis(4, at=0, labels='Estimated', cex.axis=2.5, line=0.4, tick=F)
  }
  box()
} # t

title(xlab='Easting' , line=2.6, cex.lab=2.8, outer=T)
title(ylab='Northing', line=2.4, cex.lab=2.8, outer=T)
title(main=expression("Proportion of Red-eared Slider"), line=3.2, cex.main=3, outer=T)

dev.off()


