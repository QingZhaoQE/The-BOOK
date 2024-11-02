#=================================================================================================
# Spatial dynamic N-mixture model with population growth and distance-based movement
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
kappa <- 2                        # distance effect on movement
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

lon <- cov$lon # longitude or easting
lat <- cov$lat # latitude or northing
d <- matrix(, nsite, nsite) # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1
eta <- exp(-1 * kappa * d)  # unstandardized colonization rate
theta <- eta / rowSums(eta) #   standardized colonization rate

rho <- matrix(0, nsite, nyear-1) # population growth rate
for (t in 2:nyear) {
  rho[,t-1] <- exp(cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[,t,]) %*% beta_rho)
  lambda[,t] <- (N[,t-1] * rho[,t-1]) %*% theta
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
dist_growth_Nmix_mcmc <- function(y, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:2, max)

  beta0_save <- matrix(0, ncovs+1, nmcmc) 
  beta_rho_save <- matrix(0, ncovs+2, nmcmc)
  kappa_save <- rep(0, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc)
  ppd_save <- array(0, dim=c(nsite, nyear, nreps, 2000)) # save posterior predictive density for calculating WAIC

  # Priors
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_rho_mean <- rep(0, ncovs+2)
  beta_rho_sd <- 2
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_rho <- rep(0, ncovs+2)
  kappa <- 1
  alpha <- rep(0, npcvs+1)
  lambda0 <- exp(cbind(1,x[,1,]) %*% beta0) # expectation of initial abundance
  rho <- matrix(1, nsite, nyear-1)
  eta <- exp(-1 * kappa * d)
  theta <- eta / rowSums(eta)
  N <- round((ymax + 1) / 0.5)
  p <- array(0.5, dim=c(nsite, nyear, nreps))

  # Tuning factors
  beta0_tune <- 0.08
  beta_rho_tune <- 0.15
  kappa_tune <- 0.15
  alpha_tune <- 0.02
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- matrix(rpois(nsite*nyear, N + N_tune), nsite, nyear)
    mh1 <- apply(dbinom(y, N_star, p, log=T), 1:2, sum) + 
           dpois(N_star, cbind(lambda0, t(t(N[,-nyear] * rho) %*% theta)), log=T) + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , p, log=T), 1:2, sum) + 
           dpois(N     , cbind(lambda0, t(t(N[,-nyear] * rho) %*% theta)), log=T) + 
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
    mh1 <- sum(dpois(N[,-1], t(t(N[,-nyear] * rho_star) %*% theta), log=T)) + 
           sum(dnorm(beta_rho_star, beta_rho_mean, beta_rho_sd, log=T))
    mh2 <- sum(dpois(N[,-1], t(t(N[,-nyear] * rho     ) %*% theta), log=T)) + 
           sum(dnorm(beta_rho     , beta_rho_mean, beta_rho_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_rho <- beta_rho_star
      rho <- rho_star
    }

    ### Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    eta_star <- exp(-1 * kappa_star * d)
    theta_star <- eta_star / rowSums(eta_star)
    mh1 <- sum(dpois(N[,-1], t(t(N[,-nyear] * rho) %*% theta_star), log=T)) + 
           dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=T)
    mh2 <- sum(dpois(N[,-1], t(t(N[,-nyear] * rho) %*% theta     ), log=T)) + 
           dnorm(log(kappa     ), log_kappa_mean, log_kappa_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
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

    ### Calculate posterior predictive density for calculating WAIC
    ppd <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        ppd[i,t,] <- dbinom(y[i,t,], N[i,t], p[i,t,], log=FALSE)
      } # t
    } # i

    ### Save samples
    beta0_save[,k] <- beta0
    beta_rho_save[,k] <- beta_rho
    kappa_save[k] <- kappa
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        ppd_save[,,,k-nmcmc+2000] <- ppd # save posterior predictive density for calculating WAIC
      }
    }
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_rho_save=beta_rho_save, kappa_save=kappa_save, 
       alpha_save=alpha_save, 
       ppd_save=ppd_save)

} # dist_growth_Nmix_mcmc

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
  dist_growth_Nmix_mcmc(y=y, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.3.1 dist population growth N-mixture model_output.RData')

#==============
# Plot results
#==============
pdf(file='3.3.1.1 dist population growth N-mixture model_chains.pdf', width=10, height=8)

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
  text(x=nmcmc*0.6, y=ymax, pos=1, cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=0.4, outer=T)

dev.off()

pdf(file='3.3.1.2 dist population growth N-mixture model_Nsum.pdf', width=10, height=8)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    beta0_post <- out[[i]]$beta0_save[,(nmcmc-1999):nmcmc]
    beta_rho_post <- out[[i]]$beta_rho_save[,(nmcmc-1999):nmcmc]
    kappa_post <- out[[i]]$kappa_save[(nmcmc-1999):nmcmc]
  } else {
    beta0_post <- cbind(beta0_post, out[[i]]$beta0_save[,(nmcmc-1999):nmcmc])
    beta_rho_post <- cbind(beta_rho_post, out[[i]]$beta_rho_save[,(nmcmc-1999):nmcmc])
    kappa_post <- c(kappa_post, out[[i]]$kappa_save[(nmcmc-1999):nmcmc])
  }
} # i
npost <- length(kappa_post)

lambda_pred_spatial <- lambda_pred_nonspatial <- array(, dim=c(nsite, nyear, npost))
N_pred_spatial <- N_pred_nonspatial <- array(, dim=c(nsite, nyear, npost))
lambda_pred_spatial[,1,] <- lambda_pred_nonspatial[,1,] <- exp(cbind(1,x[,1,]) %*% beta0_post)
N_pred_spatial   [,1,] <- rpois(nsite*npost, lambda_pred_spatial   [,1,])
N_pred_nonspatial[,1,] <- rpois(nsite*npost, lambda_pred_nonspatial[,1,])

eta_pred <- theta_pred <- array(, dim=c(nsite, nsite, npost))
for (k in 1:npost) {
  eta_pred[,,k] <- exp(-1 * kappa_post[k] * d)
  theta_pred[,,k] <- eta_pred[,,k] / rowSums(eta_pred[,,k])
} # k

for (k in 1:npost) {
  for (t in 2:nyear) {
    rho_pred_spatial    <- exp(cbind(1, (N_pred_spatial   [,t-1,k]-lambda_pred_spatial   [,1,k])/lambda_pred_spatial   [,1,k], x[,t,]) %*% beta_rho_post[,k])
    rho_pred_nonspatial <- exp(cbind(1, (N_pred_nonspatial[,t-1,k]-lambda_pred_nonspatial[,1,k])/lambda_pred_nonspatial[,1,k], x[,t,]) %*% beta_rho_post[,k])
    lambda_pred_spatial   [,t,k] <- t(N_pred_spatial   [,t-1,k] * rho_pred_spatial   ) %*% theta_pred[,,k]
    lambda_pred_nonspatial[,t,k] <- t(N_pred_nonspatial[,t-1,k] * rho_pred_nonspatial)
    N_pred_spatial   [,t,k] <- rpois(nsite, lambda_pred_spatial   [,t,k])
    N_pred_nonspatial[,t,k] <- rpois(nsite, lambda_pred_nonspatial[,t,k])
  } # t
} # k

Nsum_pred_spatial    <- apply(N_pred_spatial   , 2:3, sum)
Nsum_pred_nonspatial <- apply(N_pred_nonspatial, 2:3, sum)

par(mfrow=c(1,1))
par(mar=c(5,7,1,1))

plot(1, xlim=c(1,nyear), ylim=c(0, 6000), type='n', axes=F, xlab='', ylab='')
vioplot(t(Nsum_pred_spatial   ), col='lightcoral', rectCol='grey36', lineCol='grey36', colMed='yellow', border=NA, side='left',  add=T)
vioplot(t(Nsum_pred_nonspatial), col='firebrick',  rectCol='grey16', lineCol='grey16', colMed='gold',   border=NA, side='right', add=T)
lines(colSums(N), type='o', pch=16, cex=.8, lwd=1.2, col='royalblue')
axis(1, at=seq(0,nyear,5), cex.axis=1.5)
axis(2, at=seq(0,6000,1000), cex.axis=1.5, las=2)
box()
title(xlab='Year', cex.lab=2.5, line=3.5)
title(ylab='Total Population Size', cex.lab=2.5, line=4.8)

points(x=6, y=5900, pch=16, cex=1, col='royalblue')
lines(x=c(5.5,6.5), y=c(5900,5900), lwd=1.2, col='royalblue')
text(x=6.5, y=5900, labels='True', pos=4, cex=1.8)
vioplot(rnorm(10000,5450,50), at=6, col='lightcoral', rectCol='grey36', lineCol='grey36', colMed='yellow', border=NA, side='left',  add=T)
text(x=6.5, y=5450, labels='Predicted with movement', pos=4, cex=1.8)
vioplot(rnorm(10000,5000,50), at=6, col='firebrick',  rectCol='grey16', lineCol='grey16', colMed='gold',   border=NA, side='right', add=T)
text(x=6.5, y=5000, labels='Predicted without movement', pos=4, cex=1.8)

dev.off()

#================
# Calculate WAIC
#================
ppd <- matrix(, nsite * nyear * nreps, 6000)
ppd[,   1:2000] <- out[[1]]$ppd_save
ppd[,2001:4000] <- out[[2]]$ppd_save
ppd[,4001:6000] <- out[[3]]$ppd_save
lppd <- sum(log(apply(ppd, 1, mean)))
pwaic <- sum(apply(log(ppd), 1, var))
waic <- -2 * lppd + 2 * pwaic
waic <- list(lppd=lppd, pwaic=pwaic, waic=waic)
waic

# Alternative approach
library(LaplacesDemon)
WAIC(log(ppd))


