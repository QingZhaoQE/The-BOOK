#=================================================================================================
# Spatial dynamic occupancy model with distance-based movement
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
npcvs <- 1   # number of observational covariates

beta0 <- c(0.4, 1, -0.4)       # intercept and slopes for initial occupancy
beta_phi <- c(-0.2, 0.6, -0.3) # intercept and slopes for persistence
kappa <- 0.5                   # decay parameter for colonization
alpha <- c(0.3, -0.6)          # intercept and slopes for detection probability

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
eta <- exp(-1 * kappa * d)  # unstandardized colonization rate
theta <- eta / rowSums(eta) #   standardized colonization rate

phi   <- matrix(, nsite, nyear-1) # probability of persistence
gamma <- matrix(, nsite, nyear-1) # probability of colonization
for (t in 2:nyear) {
  for (i in 1:nsite) {
    phi  [i,t-1] <- inv.logit(c(1,x[i,t,]) %*% beta_phi)
    gamma[i,t-1] <- sum(theta[i,] * z[,t-1])
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

# Generate missing values
pmiss <- 0.2  # proportion of observations with missing values
nmiss <- nsite * nyear * nreps * pmiss # number of observations with missing values
miss <- array(sample(c(rep(0,nmiss),rep(1,nsite*nyear*nreps-nmiss)), nsite*nyear*nreps, replace=F), dim=c(nsite, nyear, nreps)) # identifier for missing values (0 represents missing value)
y[which(miss == 0)] <- NA # create missing values in data

#=======================
# Define MCMC algorithm
#=======================
dist_dynam_occu_mcmc <- function(y, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3] 
  npcvs <- dim(w)[4] 
  ysum <- apply(y, 1:2, sum, na.rm=T) # need na.rm=T to deal with missing values in y

  beta0_save <- matrix(0, ncovs + 1, nmcmc)
  beta_phi_save <- matrix(0, ncovs + 1, nmcmc)
  kappa_save <- rep(0, nmcmc)
  alpha_save <- matrix(0, npcvs + 1, nmcmc)

  # Priors
  beta0_mean <- rep(0, ncovs + 1)
  beta0_sd <- 2
  beta_phi_mean <- rep(0, ncovs + 1)
  beta_phi_sd <- 2
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  alpha_mean <- rep(0, npcvs + 1)
  alpha_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs + 1)
  beta_phi <- rep(0, ncovs + 1)
  kappa <- 1
  alpha <- rep(0, npcvs + 1)
  phi <- matrix(0.5, nsite, nyear-1)
  eta <- exp(-1 * kappa * d)
  theta <- eta / rowSums(eta)
  z <- ifelse(ysum==0, 0, 1)
  gamma <- matrix(, nsite, nyear-1)
  for (t in 2:nyear) {
    for (i in 1:nsite) {
      gamma[i,t-1] <- sum(theta[i,] * z[,t-1])
    } # i
  } # t
  p <- array(0.5, dim=c(nsite, nyear, nreps))
  p[which(is.na(y))] <- 0 # force detection to be 0 for missing values; this is important for sampling z

  # Tuning factors
  beta0_tune <- 0.4
  beta_phi_tune <- 0.05
  kappa_tune <- 0.1
  alpha_tune <- 0.03

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

    ### Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    eta_star <- exp(-1 * kappa_star * d)
    theta_star <- eta_star / rowSums(eta_star)
    gamma <- gamma_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      for (i in 1:nsite) {
        gamma     [i,t-1] <- sum(theta     [i,] * z[,t-1])
        gamma_star[i,t-1] <- sum(theta_star[i,] * z[,t-1])
      } # i
    } # t
    mh1 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi + (1 - z[,-nyear]) * gamma_star, log=TRUE)) + 
           dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=TRUE)
    mh2 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi + (1 - z[,-nyear]) * gamma     , log=TRUE)) + 
           dnorm(log(kappa     ), log_kappa_mean, log_kappa_sd, log=TRUE)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
      eta <- eta_star
      theta <- theta_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs + 1, alpha, alpha_tune)
    p_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        p_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_star)
      } # t
    } # i
    p_star[which(is.na(y))] <- 0 # again, force detection to be 0 for missing values
    mh1t <- mh2t <- numeric(nyear)
    for (t in 1:nyear) {
      mh1t[t] <- sum(dbinom(y[which(z[,t]==1),t,], 1, p_star[which(z[,t]==1),t,], log=TRUE), na.rm=T) # missing values don't contribute to the likelihood
      mh2t[t] <- sum(dbinom(y[which(z[,t]==1),t,], 1, p     [which(z[,t]==1),t,], log=TRUE), na.rm=T) # missing values don't contribute to the likelihood
    } # t
    mh1 <- sum(mh1t) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(mh2t) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Sample z
    s_temp <- matrix(, nsite, nyear)
    s_temp[,1] <- inv.logit(cbind(1,x[,1,]) %*% beta0)
    for (t in 2:nyear) {
      s_temp[,t] <- z[,t-1] * phi[,t-1] + (1 - z[,t-1]) * gamma[,t-1]
    } # t
    num_temp <- s_temp * apply((1 - p), 1:2, prod)
    psi_temp <- num_temp / (num_temp + (1 - s_temp))
    z[which(ysum==0)] <- rbinom(length(which(ysum==0)), 1, psi_temp[which(ysum==0)])

    ### Save samples
    beta0_save[,k] <- beta0
    beta_phi_save[,k] <- beta_phi
    kappa_save[k] <- kappa
    alpha_save[,k] <- alpha
  } # k

  # Write output
  list(beta0_save=beta0_save, 
       beta_phi_save=beta_phi_save, kappa_save=kappa_save, 
       alpha_save=alpha_save)
} # dist_dynam_occu_mcmc

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
  dist_dynam_occu_mcmc(y=y, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='2.3.2 dist dynamic occupancy model_output.RData')

#==============
# Plot results
#==============
pdf(file='2.3.2.1 dist dynamic occupancy model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["meadow"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["meadow"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]))

par(mfrow=c(3,3))
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
  title(main=ylab_beta0[i], cex.main=2, line=1.4)
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
  tt[,j] <- out[[j]]$kappa_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa - yint * 8
ymax <- kappa + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa, col='grey16', lwd=1.5)
title(main=expression(kappa), cex.main=2, line=1.2)
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
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=1, outer=T)

dev.off()

pdf(file='2.3.2.2 dist dynamic occupancy model_kappa.pdf', width=8, height=6)

kappa_post <- matrix(, nmcmc*0.8, chain)
for (j in 1:chain) {
  kappa_post[,j] <- out[[j]]$kappa_save[(nmcmc*0.2+1):nmcmc]
} # j
kappa_post <- as.vector(kappa_post)
npred <- 100
d_pred <- seq(0, 15, length.out=npred)
eta_pred <- matrix(, npred, nmcmc*0.8*chain)
for (i in 1:npred) {
  eta_pred[i,] <- exp(-1 * kappa_post * d_pred[i])
} # i
nqt <- 200
eta_qt <- apply(eta_pred, 1, quantile, probs=seq(0.001, 0.999, length.out=nqt+1))
cols <- colorRampPalette(c('white','firebrick'))(nqt/2)
cols <- c(cols, rev(cols))

par(mfrow=c(1,1))
par(mar=c(5,6,1,1))

plot(1, xlim=c(0,15), ylim=c(0,1), type='n', xlab='', ylab='', axes=F)
for (i in 1:nqt) {
  polygon(x=c(d_pred,rev(d_pred)), y=c(eta_qt[i,],rev(eta_qt[i+1,])), 
          col=cols[i], border=NA)
} # i
curve(exp(-1 * kappa * x), col='royalblue', lwd=2, add=T)
axis(1, at=seq(0,15,5.0), cex.axis=1.2)
axis(2, at=seq(0, 1,0.2), cex.axis=1.2, las=2)
box()
lines(x=c(4,5), y=c(0.95,0.95), col='royalblue', lwd=2)
for (i in 1:nqt) {
  rect(8, seq(0.9,1,length.out=nqt+1)[i], 9, seq(0.9,1,length.out=nqt+1)[i+1], 
       col=cols[i], border=NA)
} # i
text(x=5, y=0.95, labels='True', pos=4, cex=1.6)
text(x=9, y=0.95, labels='Estimated', pos=4, cex=1.6)

title(xlab='Distance', line=3.5, cex.lab=2)
title(ylab=expression(paste('Colonization Rate (', eta, ')', sep='')), line=3.8, cex.lab=2)

dev.off()


