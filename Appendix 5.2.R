#=================================================================================================
# Dynamic (i.e., multi-season) spatial capture-recapture (SCR) model with survival
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 5_Spatial capture-recapture models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(8)

# Basic values
nindv <- 100 # number of individuals (with data augmentation)
nsite <- 360 # number of sites
nyear <- 6   # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 3   # number of observational covariates

psi <- 0.8                        # probability of first-year existence
kappa <- 2.5                      # distance effect on habitat use
beta_upsilon <- c(2, 0.6, -0.3)   # intercept and slopes for habitat use
beta_phi <- c(1, -0.4, 0.6, -0.2) # intercept and slopes for apparent survival
alpha <- c(-0.6, 0.4, -0.2)       # slopes for detection rate
lon_min <- -10                    # lower bound of longitude or easting
lon_max <- +10                    # upper bound of longitude or easting
lat_min <- - 8                    # lower bound of latitude or northing
lat_max <- + 8                    # upper bound of latitude or northing

# Define functions
crossdist <- function(ss, ll) { # function for calculating distance
  c1 <- complex(real=ss[,1], imaginary=ss[,2])
  c2 <- complex(real=ll[,1], imaginary=ll[,2])
  dist <- outer(c1, c2, function(z1, z2) Mod(z1 - z2))
} # crossdist

# Read in environmental covariates
load('data/covariates.RData')
ll <- cbind(cov$lon, cov$lat) # coordinates of sites

x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- array(, dim=c(nsite, nyear, ncovs)) # environmental covariates
x[,,1] <- x1[,1:nyear]
x[,,2] <- x2[,1:nyear]

# Simulate data
### Process
ss <- matrix(, nindv, 2) # cooridnates of animal territory centroids
ss[,1] <- runif(nindv, lon_min, lon_max)
ss[,2] <- runif(nindv, lat_min, lat_max)

d <- crossdist(ss, ll) # distance between territory centroids and sites

upsilon <- array(0, dim=c(nindv, nsite, nyear)) # habitat use (i.e., time spent at a site)
for (v in 1:nindv) {
  for (t in 1:nyear) {
    upsilon[v,,t] <- exp(-1 * kappa * d[v,] + cbind(1, x[,t,]) %*% beta_upsilon)
  } # t
} # v

z <- matrix(0, nindv, nyear) # indicator of individual existence
phi <- matrix(0, nindv, nyear-1) # apparent survival probability
z[,1] <- rbinom(nindv, 1, psi) # individual existence in the first year
for (v in 1:nindv) {
  for (t in 2:nyear) {
    phi[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1])) %*% beta_phi)
    z[v,t] <- rbinom(1, 1, z[v,t-1] * phi[v,t-1])
  } # t
} # v

### Spatial capture-recapture data
w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

p <- array(0, dim=c(nsite, nyear, nreps)) # detection rate
for (i in 1:nsite) {
  for (t in 1:nyear) {
    p[i,t,] <- exp(w[i,t,,] %*% alpha)
  } # t
} # i

o <- ifelse(ll[,1] <= 8 & ll[,1] >= -8 & ll[,2] <= 6 & ll[,2] >= -6, 1, 0) # indicator of sites with camera traps

y <- array(0, dim=c(nindv, nsite, nyear, nreps)) # frequency of capture data
for (v in 1:nindv) {
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      y[v,i,t,] <- rpois(nreps, upsilon[v,i,t] * z[v,t] * p[i,t,] * o[i])
    } # t
  } # i
} # v

#=======================
# Define MCMC algorithm
#=======================
scr_survival_mcmc <- function(y, x, w, o, ll, lon_min, lon_max, lat_min, lat_max, nmcmc) {

  # Setup variables
  nindv <- dim(y)[1]
  nsite <- dim(y)[2]
  nyear <- dim(y)[3]
  nreps <- dim(y)[4]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:3, max)
  obs01 <- ifelse(apply(y, c(1,3), sum) > 0, 1, 0) # indicate if an individual is observed at all or not

  psi_save <- rep(0, nmcmc)
  kappa_save <- rep(0, nmcmc)
  beta_upsilon_save <- matrix(0, ncovs+1, nmcmc)
  beta_phi_save <- matrix(0, ncovs+2, nmcmc)
  alpha_save <- matrix(0, npcvs, nmcmc)
  z_save <- array(, dim=c(nindv, nyear, 2000))

  # Priors
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  beta_upsilon_mean <- 0
  beta_upsilon_sd <- 2
  beta_phi_mean <- 0
  beta_phi_sd <- 2
  alpha_mean <- 0
  alpha_sd <- 2

  # Starting values
  psi <- 0.5
  kappa <- 1
  beta_upsilon <- rep(0, ncovs+1)
  beta_phi <- rep(0, ncovs+2)
  alpha <- rep(0, npcvs)
  ss <- matrix(0, nindv, 2)
  ss[,1] <- runif(nindv, lon_min, lon_max)
  ss[,2] <- runif(nindv, lat_min, lat_max)
  for (v in 1:nindv) {
    if (sum(obs01[v,]) > 0) {
      ss[v,1] <- mean(ll[which(rowSums(ymax[v,,]) > 0), 1])
      ss[v,2] <- mean(ll[which(rowSums(ymax[v,,]) > 0), 2])
    }
  } # v
  d <- crossdist(ss, ll)
  upsilon <- array(0, dim=c(nindv, nsite, nyear))
  for (v in 1:nindv) {
    for (t in 1:nyear) {
      upsilon[v,,t] <- exp(-1 * kappa * d[v,] + cbind(1, x[,t,]) %*% beta_upsilon)
    } # t
  } # v
  phi <- matrix(0.5, nindv, nyear-1)
  p <- array(0.5, dim=c(nsite, nyear, nreps))
  zi <- obs01
  for (v in 1:nindv) {
    if (sum(zi[v,]) > 0) {
      zi[v,1:max(which(zi[v,] == 1))] <- 1
    }
  } # v
  z <- zi

  # Tuning factors
  kappa_tune <- 0.15
  beta_upsilon_tune <- 0.05
  beta_phi_tune <- 0.25
  alpha_tune <- 0.05
  ss_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample psi
    psi <- rbeta(1, sum(z[,1]) + 0.001, sum(1-z[,1]) + 1)

    ### Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    upsilon_star <- array(0, dim=c(nindv, nsite, nyear))
    for (v in 1:nindv) {
      for (t in 1:nyear) {
        upsilon_star[v,,t] <- exp(-1 * kappa_star * d[v,] + cbind(1, x[,t,]) %*% beta_upsilon)
      } # t
    } # v
    prob1 <- prob2 <- matrix(0, nindv, nyear)
    for (v in 1:nindv) {
      for (t in 1:nyear) {
        prob1[v,t] <- sum(dpois(y[v,,t,], upsilon_star[v,,t] * z[v,t] * p[,t,] * o, log=TRUE))
        prob2[v,t] <- sum(dpois(y[v,,t,], upsilon     [v,,t] * z[v,t] * p[,t,] * o, log=TRUE))
      } # t
    } # v
    mh1 <- sum(prob1) + 
           dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=TRUE)
    mh2 <- sum(prob2) + 
           dnorm(log(kappa     ), log_kappa_mean, log_kappa_sd, log=TRUE)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
      upsilon <- upsilon_star
    }

    ### Sample beta_upsilon
    beta_upsilon_star <- rnorm(ncovs+1, beta_upsilon, beta_upsilon_tune)
    upsilon_star <- array(0, dim=c(nindv, nsite, nyear))
    for (v in 1:nindv) {
      for (t in 1:nyear) {
        upsilon_star[v,,t] <- exp(-1 * kappa * d[v,] + cbind(1, x[,t,]) %*% beta_upsilon_star)
      } # t
    } # v
    prob1 <- prob2 <- matrix(0, nindv, nyear)
    for (v in 1:nindv) {
      for (t in 1:nyear) {
        prob1[v,t] <- sum(dpois(y[v,,t,], upsilon_star[v,,t] * z[v,t] * p[,t,] * o, log=TRUE))
        prob2[v,t] <- sum(dpois(y[v,,t,], upsilon     [v,,t] * z[v,t] * p[,t,] * o, log=TRUE))
      } # t
    } # v
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_upsilon_star, beta_upsilon_mean, beta_upsilon_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta_upsilon     , beta_upsilon_mean, beta_upsilon_sd, log=TRUE))
    mhB <- exp(mh1 - mh2)
    if (mhB > runif(1)) {
      beta_upsilon <- beta_upsilon_star
      upsilon <- upsilon_star
    }

    ### Sample beta_phi
    beta_phi_star <- rnorm(ncovs+2, beta_phi, beta_phi_tune)
    phi_star <- matrix(0, nindv, nyear-1)
    for (v in 1:nindv) {
      for (t in 2:nyear) {
        phi_star[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1])) %*% beta_phi_star)
      } # t
    } # v
    mh1 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi_star, log=TRUE)) + 
           sum(dnorm(beta_phi_star, beta_phi_mean, beta_phi_sd, log=TRUE))
    mh2 <- sum(dbinom(z[,-1], 1, z[,-nyear] * phi     , log=TRUE)) + 
           sum(dnorm(beta_phi     , beta_phi_mean, beta_phi_sd, log=TRUE))
    mhP <- exp(mh1 - mh2)
    if (mhP > runif(1)) {
      beta_phi <- beta_phi_star
      phi <- phi_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs, alpha, alpha_tune)
    p_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        p_star[i,t,] <- exp(w[i,t,,] %*% alpha_star)
      } # t
    } # i
    prob1 <- prob2 <- matrix(0, nindv, nyear)
    for (v in 1:nindv) {
      for (t in 1:nyear) {
        prob1[v,t] <- sum(dpois(y[v,,t,], upsilon[v,,t] * z[v,t] * p_star[,t,] * o, log=TRUE))
        prob2[v,t] <- sum(dpois(y[v,,t,], upsilon[v,,t] * z[v,t] * p     [,t,] * o, log=TRUE))
      } # t
    } # v
    mh1 <- sum(prob1) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=TRUE))
    mhA <- exp(mh1 - mh2)
    if (mhA > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Sample ss
    ss_star <- matrix(rnorm(nindv*2, ss, ss_tune), nindv, 2)
    d_star <- crossdist(ss_star, ll)
    upsilon_star <- array(0, dim=c(nindv, nsite, nyear))
    for (v in 1:nindv) {
      for (t in 1:nyear) {
        upsilon_star[v,,t] <- exp(-1 * kappa * d_star[v,] + cbind(1, x[,t,]) %*% beta_upsilon)
      } # t
    } # v
    for (v in 1:nindv) {
      if (ss_star[v,1] > lon_min & ss_star[v,1] < lon_max & 
          ss_star[v,2] > lat_min & ss_star[v,2] < lat_max) {
        prob1 <- prob2 <- rep(0, nyear)
        for (t in 1:nyear) {
          prob1[t] <- sum(dpois(y[v,,t,], upsilon_star[v,,t] * z[v,t] * p[,t,] * o, log=TRUE))
          prob2[t] <- sum(dpois(y[v,,t,], upsilon     [v,,t] * z[v,t] * p[,t,] * o, log=TRUE))
        } # t
        mh1 <- sum(prob1)
        mh2 <- sum(prob2)
        mhS <- exp(mh1 - mh2)
        if (mhS > runif(1)) {
          ss[v,] <- ss_star[v,]
          d[v,] <- d_star[v,]
          upsilon[v,,] <- upsilon_star[v,,]
        }
      }
    } # v

    ### Sample z
    for (v in 1:nindv) {
      if (zi[v,1] == 0) {
        prod_temp <- psi * prod(exp(-1 * upsilon[v,,1] * p[,1,] * o))
        psi_temp <- prod_temp / (prod_temp + 1 - psi)
        z[v,1] <- rbinom(1, 1, psi_temp)
      }
      for (t in 2:nyear) {
        if (z[v,t-1] == 0 & zi[v,t] == 0) {
          z[v,t] <- 0
        } else if (z[v,t-1] == 1 & zi[v,t] == 0) {
          prod_temp <- phi[v,t-1] * prod(exp(-1 * upsilon[v,,t] * p[,t,] * o))
          psi_temp <- prod_temp / (prod_temp + 1 - phi[v,t-1])
          z[v,t] <- rbinom(1, 1, psi_temp)
        }
      } # t
    } # v

    ### Save samples
    psi_save[k] <- psi
    kappa_save[k] <- kappa
    beta_upsilon_save[,k] <- beta_upsilon
    beta_phi_save[,k] <- beta_phi
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        z_save[,,k-nmcmc+2000] <- z
      }
    }
  } # k

  ### Write output
  list(psi_save=psi_save, 
       kappa_save=kappa_save, beta_upsilon_save=beta_upsilon_save, 
       beta_phi_save=beta_phi_save, 
       alpha_save=alpha_save, 
       z_save=z_save)

} # scr_survival_mcmc

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
  scr_survival_mcmc(y=y, x=x, w=w, o=o, ll=ll, 
                    lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, 
                    nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='5.2.1 dynamic scr with survival_output.RData')

#==============
# Plot results
#==============
pdf(file='5.2.1.1 dynamic scr with survival_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta_upsilon <- c(expression(beta[0]^""["["*upsilon*"]"]), expression(beta["bamboo"]^""["["*upsilon*"]"]), 
                       expression(beta["temperature"]^""["["*upsilon*"]"]))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["density"]^""["["*phi*"]"]), 
                   expression(beta["bamboo"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]))
ylab_alpha <- c(expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))

par(mfrow=c(4,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,0))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$psi_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- psi - yint * 8
ymax <- psi + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=psi, col='grey16', lwd=1.5)
title(main=expression(psi), cex.main=2, line=1.2)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$kappa_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa - yint * 2
ymax <- kappa + yint * 2
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa, col='grey16', lwd=1.5)
title(main=expression(kappa), cex.main=2, line=1.2)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_upsilon_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_upsilon[i] - yint * c(2,2,2)[i]
  ymax <- beta_upsilon[i] + yint * c(2,2,2)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(1,1,1)[i]), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_upsilon[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_upsilon[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi[i] - yint * c(8,8,8,8)[i]
  ymax <- beta_phi[i] + yint * c(8,8,8,8)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(4,4,4,4)[i]), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:npcvs) {
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

pdf(file='5.2.1.2 dynamic scr with survival_z.pdf', width=8, height=8)

z_post <- array(, dim=c(nindv, nyear, 6000))
z_post[,,   1:2000] <- out[[1]]$z
z_post[,,2001:4000] <- out[[2]]$z
z_post[,,4001:6000] <- out[[3]]$z
z_mean <- apply(z_post, 1:2, mean)
nscale <- 100
z_scale <- round(z_mean * (nscale - 1)) + 1
col_scale <- colorRampPalette(c('white','lightcoral'))(nscale)
z_order <- z[rev(order(rowSums(z))),]
z_scale_order <- z_scale[rev(order(rowSums(z))),]

par(mfrow=c(2,1))
par(mar=c(1,0,2.5,1))
par(oma=c(4,5.5,0,0))

plot(1, xlim=c(.5,nyear+.5), ylim=c(.5,nindv+.5), type='n', axes=F)
for (i in 1:nindv) {
  for (t in 1:nyear) {
    rect(t-.5,i-.5,t+.5,i+.5, col=ifelse(z_order[i,t]==1, 'royalblue', 'white'), border=NA)
  } # t
} # i
axis(1, at=seq(1,nyear,1), labels=rep('',nyear))
axis(2, at=seq(10,nindv,20), cex.axis=1.5, las=2)
box()
title(main='True Survival History', line=.5, cex.main=2)

plot(1, xlim=c(.5,nyear+.5), ylim=c(.5,nindv+.5), type='n', axes=F)
for (i in 1:nindv) {
  for (t in 1:nyear) {
    rect(t-.5,i-.5,t+.5,i+.5, col=col_scale[z_scale_order[i,t]], border=NA)
  } # t
} # i
axis(1, at=seq(1,nyear,1), cex.axis=1.5)
axis(2, at=seq(10,nindv,20), cex.axis=1.5, las=2)
box()
title(main='Estimated Survival History', line=.5, cex.main=2)

title(ylab='Individual', cex.lab=3, line=3.2, outer=T)
title(xlab='Year', cex.lab=3, line=2.8, outer=T)

dev.off()


