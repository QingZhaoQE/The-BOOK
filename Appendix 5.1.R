#=================================================================================================
# Static (i.e., single-season) spatial capture-recapture (SCR) model
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
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 4   # number of observational covariates

psi <- 0.8                       # probability of existence
kappa <- 2.5                     # distance effect on habitat use
beta <- c(2, 0.6, -0.3)          # intercept and slopes for habitat use
alpha <- c(0.8, -0.4, 0.2, -0.1) # slopes for detection rate
lon_min <- -10                   # lower bound of longitude or easting
lon_max <- +10                   # upper bound of longitude or easting
lat_min <- - 8                   # lower bound of latitude or northing
lat_max <- + 8                   # upper bound of latitude or northing

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
x <- cbind(x1[,1], x2[,1]) # environmental covariates

# Simulate data
### Process
ss <- matrix(, nindv, 2) # cooridnates of animal territory centroids
ss[,1] <- runif(nindv, lon_min, lon_max)
ss[,2] <- runif(nindv, lat_min, lat_max)

d <- crossdist(ss, ll) # distance between territory centroids and sites

upsilon <- matrix(0, nindv, nsite) # habitat use (i.e., time spent at a site)
for (v in 1:nindv) {
  upsilon[v,] <- exp(-1 * kappa * d[v,] + cbind(1, x) %*% beta)
} # v

z <- rbinom(nindv, 1, psi) # indicator of individual existence
N <- sum(z) # population size, i.e., number of individuals that actually exist

### Spatial capture-recapture data
w <- array(rnorm(nsite * nreps * npcvs, 0, 1), dim=c(nsite, nreps, npcvs)) # observational covariates

p <- matrix(0, nsite, nreps) # detection rate
for (i in 1:nsite) {
  p[i,] <- exp(w[i,,] %*% alpha)
} # i

o <- ifelse(ll[,1] <= 8 & ll[,1] >= -8 & ll[,2] <= 6 & ll[,2] >= -6, 1, 0) # indicator of sites with camera traps

y <- array(0, dim=c(nindv, nsite, nreps)) # frequency of capture data
for (v in 1:nindv) {
  for (i in 1:nsite) {
    y[v,i,] <- rpois(nreps, upsilon[v,i] * z[v] * p[i,] * o[i])
  } # i
} # v

#=======================
# Define MCMC algorithm
#=======================
scr_mcmc <- function(y, x, w, o, ll, lon_min, lon_max, lat_min, lat_max, nmcmc) {

  # Setup variables
  nindv <- dim(y)[1]
  nsite <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[2]
  npcvs <- dim(w)[3]
  ymax <- apply(y, 1:2, max)
  obs01 <- ifelse(apply(y, 1, sum) > 0, 1, 0) # indicate if an individual is observed at all or not

  psi_save <- rep(0, nmcmc)
  kappa_save <- rep(0, nmcmc)
  beta_save <- matrix(0, ncovs+1, nmcmc)
  alpha_save <- matrix(0, npcvs, nmcmc)
  ss_save <- array(, dim=c(nindv, 2, 2000))
  N_save <- rep(0, 2000)

  # Priors
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  beta_mean <- 0
  beta_sd <- 2
  alpha_mean <- 0
  alpha_sd <- 2

  # Starting values
  psi <- 0.5
  kappa <- 1
  beta <- rep(0, ncovs+1)
  alpha <- rep(0, npcvs)
  ss <- matrix(0, nindv, 2)
  ss[,1] <- runif(nindv, lon_min, lon_max)
  ss[,2] <- runif(nindv, lat_min, lat_max)
  for (v in 1:nindv) {
    if (obs01[v] == 1) {
      ss[v,1] <- mean(ll[which(ymax[v,] > 0), 1])
      ss[v,2] <- mean(ll[which(ymax[v,] > 0), 2])
    }
  } # v
  d <- crossdist(ss, ll)
  upsilon <- matrix(, nindv, nsite)
  for (v in 1:nindv) {
    upsilon[v,] <- exp(-1 * kappa * d[v,] + cbind(1, x) %*% beta)
  } # v
  p <- matrix(1, nsite, nreps)
  z <- obs01

  # Tuning factors
  kappa_tune <- 0.2
  beta_tune <- 0.1
  alpha_tune <- 0.05
  ss_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample psi
    psi <- rbeta(1, sum(z) + 0.001, sum(1-z) + 1)

    ### Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    upsilon_star <- matrix(, nindv, nsite)
    for (v in 1:nindv) {
      upsilon_star[v,] <- exp(-1 * kappa_star * d[v,] + cbind(1, x) %*% beta)
    } # v
    prob1 <- prob2 <- rep(0, nindv)
    for (v in 1:nindv) {
      prob1[v] <- sum(dpois(y[v,,], upsilon_star[v,] * z[v] * p * o, log=TRUE))
      prob2[v] <- sum(dpois(y[v,,], upsilon     [v,] * z[v] * p * o, log=TRUE))
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

    ### Sample beta
    beta_star <- rnorm(ncovs+1, beta, beta_tune)
    upsilon_star <- matrix(, nindv, nsite)
    for (v in 1:nindv) {
      upsilon_star[v,] <- exp(-1 * kappa * d[v,] + cbind(1, x) %*% beta_star)
    } # v
    prob1 <- prob2 <- rep(0, nindv)
    for (v in 1:nindv) {
      prob1[v] <- sum(dpois(y[v,,], upsilon_star[v,] * z[v] * p * o, log=TRUE))
      prob2[v] <- sum(dpois(y[v,,], upsilon     [v,] * z[v] * p * o, log=TRUE))
    } # v
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_star, beta_mean, beta_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta     , beta_mean, beta_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta <- beta_star
      upsilon <- upsilon_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs, alpha, alpha_tune)
    p_star <- matrix(, nsite, nreps)
    for (i in 1:nsite) {
      p_star[i,] <- exp(w[i,,] %*% alpha_star)
    } # i
    prob1 <- prob2 <- matrix(0, nindv, nsite)
    prob1 <- prob2 <- rep(0, nindv)
    for (v in 1:nindv) {
      prob1[v] <- sum(dpois(y[v,,], upsilon[v,] * z[v] * p_star * o, log=TRUE))
      prob2[v] <- sum(dpois(y[v,,], upsilon[v,] * z[v] * p      * o, log=TRUE))
    } # v
    mh1 <- sum(prob1) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Sample ss
    ss_star <- matrix(rnorm(nindv*2, ss, ss_tune), nindv, 2)
    d_star <- crossdist(ss_star, ll)
    upsilon_star <- matrix(, nindv, nsite)
    for (v in 1:nindv) {
      upsilon_star[v,] <- exp(-1 * kappa * d_star[v,] + cbind(1, x) %*% beta)
    } # v
    for (v in 1:nindv) {
      if (ss_star[v,1] > lon_min & ss_star[v,1] < lon_max & 
          ss_star[v,2] > lat_min & ss_star[v,2] < lat_max) {
        mh1 <- sum(dpois(y[v,,], upsilon_star[v,] * z[v] * p * o, log=TRUE))
        mh2 <- sum(dpois(y[v,,], upsilon     [v,] * z[v] * p * o, log=TRUE))
        mh <- exp(mh1 - mh2)
        if (mh > runif(1)) {
          ss[v,] <- ss_star[v,]
          d[v,] <- d_star[v,]
          upsilon[v,] <- upsilon_star[v,]
        }
      }
    } # v

    ### Sample z
    for (v in 1:nindv) {
      if (obs01[v] == 0) {
        prod_temp <- psi * prod(exp(-1 * upsilon[v,] * p * o))
        psi_temp <- prod_temp / (prod_temp + 1 - psi)
        z[v] <- rbinom(1, 1, psi_temp)
      }
    } # v

    ### Save samples
    psi_save[k] <- psi
    kappa_save[k] <- kappa
    beta_save[,k] <- beta
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        ss_save[,,k-nmcmc+2000] <- ss
        N_save[k-nmcmc+2000] <- sum(z)
      }
    }
  } # k

  ### Write output
  list(psi_save=psi_save, 
       kappa_save=kappa_save, beta_save=beta_save, 
       alpha_save=alpha_save, 
       ss_save=ss_save, N_save=N_save)

} # scr_mcmc

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
  scr_mcmc(y=y, x=x, w=w, o=o, ll=ll, 
           lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, 
           nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='5.1 static scr_output.RData')

#==============
# Plot results
#==============
pdf(file='5.1.1 static scr_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta <- c(expression(beta[0]), expression(beta["bamboo"]), expression(beta["temperature"]))
ylab_alpha <- c(expression(alpha[1]), expression(alpha[2]), expression(alpha[3]), expression(alpha[4]))

par(mfrow=c(3,3))
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
ymin <- kappa - yint * 4
ymax <- kappa + yint * 4
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
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
    tt[,j] <- out[[j]]$beta_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta[i] - yint * 4
  ymax <- beta[i] + yint * 4
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
  if (i == 1) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
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

pdf(file='5.1.2 static scr_ss.pdf', width=8.5, height=8)

library(vioplot) # for making violin plots

ss_post <- array(, dim=c(nindv, 2, 6000))
ss_post[,,   1:2000] <- out[[1]]$ss
ss_post[,,2001:4000] <- out[[2]]$ss
ss_post[,,4001:6000] <- out[[3]]$ss
ss_sel <- ss_post[,,sample(1:6000, 100, replace=F)]
N_post <- c(out[[1]]$N, out[[2]]$N, out[[3]]$N)

layout(mat=matrix(c(rep(c(1,1,2,2),times=3),rep(3,8)), nrow=4, ncol=5, byrow=FALSE))
par(mar=c(1,2,3,1))
par(oma=c(4,4,0,6))

plot(1, xlim=c(-10,10), ylim=c(-8,8), axes=F, xlab='', ylab='', type='n')
for (v in 1:N) {
  col <- col2rgb(v)
  for (i in 1:nsite) {
    if (upsilon[v,i] > 0.5) {
      points(x=ll[i,1], y=ll[i,2], pch=16, cex=sqrt(upsilon[v,i])*1.4, 
             col=rgb(col[1],col[2],col[3],alpha=120,maxColorValue=255))
      lines(x=c(ss[v,1],ll[i,1]), y=c(ss[v,2],ll[i,2]), lwd=1, col=v)
    }
  } # i
  points(ss[v,1], ss[v,2], cex=ifelse(sum(y[v,,])==0,1.6,1.2), pch=16, col=v)
} # v
points(ll[which(o==1),], lwd=0.5, pch=3, col='grey50', cex=0.8)
axis(1, at=seq(-10,10,5), labels=rep('',5))
axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
title(main='True Centroid & Habitat Use', cex.main=2)
box()

plot(1, xlim=c(-10,10), ylim=c(-8,8), axes=F, xlab='', ylab='', type='n')
for (v in 1:N) {
  points(t(ss_sel[v,,]), pch=16, cex=3, col=rgb(240,128,128,alpha=6,maxColorValue=255))
} # v
for (v in 1:N) {
  points(ss[v,1], ss[v,2], cex=ifelse(sum(y[v,,])==0,1.6,1.2), pch=ifelse(sum(y[v,,])==0,21,16), col='royalblue', bg='white')
} # v
points(ll[which(o==1),], lwd=0.5, pch=3, col='grey50', cex=0.8)
axis(1, at=seq(-10,10,5), cex.axis=1.6)
axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
title(main='True & Estimated Centroid', cex.main=2)
box()

axis(1, at=0, labels='Easting', tick=F, line=2.5, cex.axis=3)
title(ylab='Northing', cex.lab=3, line=1.5, outer=T)

par(mar=c(1,1,3,1))
plot(1, ylim=c(75, 100), type='n', axes=F, xlab='', ylab='')
vioplot(N_post, col='lightcoral', rectCol='grey36', lineCol='grey36', border=NA, cex=3, add=T)
abline(h=N, col='royalblue', lwd=1.8)
axis(4, at=seq(75, 100, 5), cex.axis=1.6, las=2)
axis(4, at=87.5, label='Number of Individuals', line=4, cex.axis=3, tick=F)

dev.off()


