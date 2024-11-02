#=================================================================================================
# Multistate model
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 4_Integrated population models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions
library(LaplacesDemon) # for generating random numbers from categorical distributions (rcat)

set.seed(6)

# Basic values
nindv <- 10800 # number of captured individuals
nsite <- 60    # number of sites
nyear <- 10    # number of years
ncovs <- 2     # number of covariates
npcvs <- 1     # number of capture covariates

beta <-  c( 0.4,  0.6, -0.4) # intercept and slopes for survival probability
kappa <- 0.5                 # distance effect on movement
alpha <- c(-1.4, -0.8)       # intercept and slopes for capture probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- array(, dim=c(nsite, nyear, ncovs)) # environmental covariates
x[,,1] <- x1[1:nsite, 1:nyear]
x[,,2] <- x2[1:nsite, 1:nyear]

# Simulate data
### Process
phi <- matrix(, nsite, nyear-1) # survival probability
for (t in 1:(nyear-1)) {
  phi[,t] <- inv.logit(cbind(1,x[,t,]) %*% beta)
} # t

lon <- cov$lon[1:nsite] # longitude or easting
lat <- cov$lat[1:nsite] # latitude or northing
d <- matrix(, nsite, nsite)  # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1
eta <- exp(-1 * kappa * d)  # unstandardized colonization rate
theta <- eta / rowSums(eta) #   standardized colonization rate

year_f <- rep(1:(nyear-1), each=nindv/(nyear-1))                       # year of first capture
site_f <- rep(rep(1:nsite, each=nindv/nsite/(nyear-1)), times=nyear-1) # site of first capture

z <- matrix(0, nindv, nyear) # individual survival-location history
for (v in 1:nindv) {
  z[v,year_f[v]] <- site_f[v]
  for (t in (year_f[v]+1):nyear) {
    if (z[v,t-1] %in% c(1:nsite)) {
      z[v,t] <- rcat(1, p=c(phi[z[v,t-1],t-1] * theta[z[v,t-1],], 1 - phi[z[v,t-1],t-1]))
    } else {
      z[v,t] <- 0
    }
    if (!(z[v,t] %in% 1:nsite)) {
      z[v,t] <- 0
    }
  } # t
} # v

### Capture
w <- array(rnorm(nsite * (nyear-1) * npcvs, 0, 1), dim=c(nsite, nyear-1, npcvs)) # capture covariates

p <- matrix(, nsite, nyear-1) # capture probability
for (t in 1:(nyear-1)) {
  p[,t] <- inv.logit(cbind(1,w[,t,]) %*% alpha)
} # t

c <- matrix(0, nindv, nyear) # individual capture history
for (v in 1:nindv) {
  c[v,year_f[v]] <- z[v,year_f[v]]
  for (t in (year_f[v]+1):nyear) {
    if (z[v,t] %in% c(1:nsite)) {
      c[v,t] <- ifelse(rbinom(1, 1, p[z[v,t],t-1]) == 1, z[v,t], 0)
    } else {
      c[v,t] <- 0
    }
  } # t
} # v

# Extend individual capture history
for (v in 1:nindv) {
  t1 <- c[v,]
  ycap <- which(t1 %in% c(1:nsite))
  ncap <- length(ycap)
  if (ncap == 1) {
    t2 <- t1
  } else {
    t2 <- matrix(, ncap, nyear)
    for (h in 1:(ncap-1)) {
      t2[h, ycap[h]:ycap[h+1]] <- t1[ycap[h]:ycap[h+1]]
    } # h
    t2[ncap, ycap[ncap]] <- t1[ycap[ncap]]
  }
  if (v == 1) {
    c_ext <- t2
  } else {
    c_ext <- rbind(c_ext, t2)
  }
} # v
c_ext[which(is.na(c_ext))] <- 0

year_f_ext <- site_f_ext <- numeric(dim(c_ext)[1])
for (v in 1:dim(c_ext)[1]) {
  year_f_ext[v] <- min(which(c_ext[v,] %in% c(1:nsite)))
  site_f_ext[v] <- c_ext[v, year_f_ext[v]]
} # v

# Converting individual capture history to m-array
for (t in 2:nyear) {
  tt <- paste('site', '_', 1:nsite, ', year', t, sep='')
  if (t == 2) {
    nm <- tt
  } else { 
    nm <- c(nm, tt)
  }
} # t
nm <- c(nm, 'not seen')

m <- array(0, # m-array of capture-recapture
  dim=c(nyear-1, nsite*(nyear-1)+1, nsite), 
  dimnames=list(paste('year', 1:(nyear-1), sep='_'), nm, paste('site', 1:nsite, sep='_')))
for (t1 in 1:(nyear-1)) {
  temp1 <- c_ext[which(year_f_ext == t1),]
  for (i1 in 1:nsite) {
    temp2 <- temp1[which(temp1[,t1] == i1),]
    nband <- length(which(temp1[,t1] == i1))
    for (t2 in (t1+1):nyear) {
      for (i2 in 1:nsite) {
        m[t1,(t2-2)*nsite+i2,i1] <- length(which(temp2[,t2] == i2))
      } # i2
    } # t2
    m[t1,nsite*(nyear-1)+1,i1] <- nband - sum(m[t1,1:(nsite*(nyear-1)),i1])
  } # i1
} # t1

#=======================
# Define MCMC algorithm
#=======================
# Define the function to calculate the cell probability of m-array
cell_prob_fun <- function(nsite, nyear, phi, theta, p) {
  q <- 1 - p
  cell_prob_p <- cell_prob_q <- array(0, dim=c(nyear-1, nsite*(nyear-1)+1, nsite))
  for (i in 1:nsite) {
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,nsite*(t-1)+1:nsite,i] <- phi[i,t] * theta[i,] * p[,t]
      cell_prob_q[t,nsite*(t-1)+1:nsite,i] <- phi[i,t] * theta[i,] * q[,t]
    } # t
    for (t1 in 1:(nyear-2)) {
      for (t2 in (t1+1):(nyear-1)) {
        cell_prob_p[t1,nsite*(t2-1)+1:nsite,i] <- ((cell_prob_q[t1,nsite*(t2-2)+1:nsite,i] * phi[,t2]) %*% theta) * p[,t2]
        cell_prob_q[t1,nsite*(t2-1)+1:nsite,i] <- ((cell_prob_q[t1,nsite*(t2-2)+1:nsite,i] * phi[,t2]) %*% theta) * q[,t2]
      } # t2
    } # t1
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,nsite*(nyear-1)+1,i] <- 1 - sum(cell_prob_p[t,1:(nsite*(nyear-1)),i])
    } # t
  } # i
  cell_prob <- cell_prob_p

  return(cell_prob)
} # cell_prob_fun

multistate_mcmc <- function(m, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(x)[1]
  nyear <- dim(x)[2]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[3]

  beta_save <- matrix(0, ncovs+1, nmcmc)
  kappa_save <- rep(0, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc)

  # Prior
  beta_mean <- rep(0, ncovs+1)
  beta_sd <- 2
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta <- rep(0, ncovs+1)
  kappa <- 1
  alpha <- rep(0, npcvs+1)
  phi <- matrix(0.5, nsite, nyear-1)
  eta <- exp(-1 * kappa * d)
  theta <- eta / rowSums(eta)
  p <- matrix(0.5, nsite, nyear-1) # capture probability
  cell_prob <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi, theta=theta, p=p)

  # Tuning factor
  beta_tune <- 0.035
  kappa_tune <- 0.035
  alpha_tune <- 0.035

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    # Sample beta
    beta_star <- rnorm(ncovs+1, beta, beta_tune)
    phi_star <- matrix(, nsite, nyear-1)
    for (t in 1:(nyear-1)) {
      phi_star[,t] <- inv.logit(cbind(1,x[,t,]) %*% beta_star)
    } # t
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi_star, theta=theta, p=p)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_star, beta_mean, beta_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta     , beta_mean, beta_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta <- beta_star
      phi <- phi_star
      cell_prob <- cell_prob_star
    }

    # Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    eta_star <- exp(-1 * kappa_star * d)
    theta_star <- eta_star / rowSums(eta_star)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi, theta=theta_star, p=p)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=T)
    mh2 <- sum(prob2) + 
           dnorm(log(kappa     ), log_kappa_mean, log_kappa_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
      eta <- eta_star
      theta <- theta_star
      cell_prob <- cell_prob_star
    }

    # Sample alpha
    alpha_star <- rnorm(npcvs+1, alpha, alpha_tune)
    p_star <- matrix(, nsite, nyear-1)
    for (t in 1:(nyear-1)) {
      p_star[,t] <- inv.logit(cbind(1,w[,t,]) %*% alpha_star)
    } # t
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi, theta=theta, p=p_star)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
      cell_prob <- cell_prob_star
    }

    # Save samples
    beta_save[,k] <- beta
    kappa_save[k] <- kappa
    alpha_save[,k] <- alpha
  } # k

  list(beta_save=beta_save, kappa_save=kappa_save, 
       alpha_save=alpha_save)
} # multistate_mcmc

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
  multistate_mcmc(m=m, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='4.2 multistate model_output.RData')

#==============
# Plot results
#==============
pdf(file='4.2.1 multistate model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta <- c(expression(beta[0]), expression(beta["reef"]), expression(beta["temperature"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]))

par(mfrow=c(2,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,0))

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
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta[i], col='grey16', lwd=1.5)
  title(main=ylab_beta[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

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
axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa, col='grey16', lwd=1.5)
title(main=expression(kappa), cex.main=2, line=1)
text(x=nmcmc*0.6, y=ymax, pos=1,  cex=1.8, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i] - yint * 4
  ymax <- alpha[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.6, y=ymax, pos=1, cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=1, outer=T)

dev.off()

pdf(file='4.2.2 multistate model_eta.pdf', width=8, height=12)

kappa_post <- c(out[[1]]$kappa[(nmcmc*0.2+1):nmcmc], out[[2]]$kappa[(nmcmc*0.2+1):nmcmc], out[[3]]$kappa[(nmcmc*0.2+1):nmcmc])
kappa_med <- median(kappa_post)

eta_true <- exp(-1 * kappa     * d)
eta_pred <- exp(-1 * kappa_med * d)

min_eta <- 0.2
max_lwd <- 5

par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
par(oma=c(4,5,2,2))

plot(lon, lat, xlim=c(-10.2,10.2), ylim=c(-8.2,8.2), axes=F, xlab='', ylab='', 
     pch=21, cex=rowMeans(cov$size)[1:nsite]*4.6, col='grey18', bg='white')
axis(1, at=seq(-10,10,5), labels=rep('',5))
axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
box()
for (i1 in 1:(nsite-1)) {
  for (i2 in (i1+1):nsite) {
    if (eta_true[i1,i2] >= min_eta) {
      lines(x=lon[c(i1,i2)], y=lat[c(i1,i2)], col='royalblue', lwd=eta_true[i1,i2]*max_lwd)
    }
  } # i2
} # i1
axis(3, at=0, labels=expression(paste('Movement Rate (', eta, ')', sep='')), 
     cex.axis=2.5, tick=F, line=-0.6)
axis(4, at=0, labels='True', cex.axis=2.5, tick=F, line=0.5)

plot(lon, lat, xlim=c(-10.2,10.2), ylim=c(-8.2,8.2), axes=F, xlab='', ylab='', 
     pch=21, cex=rowMeans(cov$size)[1:nsite]*4.6, col='grey18', bg='white')
axis(1, at=seq(-10,10,5), cex.axis=1.6)
axis(2, at=seq(-08,08,4), cex.axis=1.6, las=2)
box()
for (i1 in 1:(nsite-1)) {
  for (i2 in (i1+1):nsite) {
    if (eta_pred[i1,i2] >= min_eta) {
      lines(x=lon[c(i1,i2)], y=lat[c(i1,i2)], col='lightcoral', lwd=eta_pred[i1,i2]*max_lwd)
    }
  } # i2
} # i1
axis(4, at=0, labels='Estimated', cex.axis=2.5, tick=F, line=0.5)

title(xlab='Easting' , line=2.6, cex.lab=2.8, outer=T)
title(ylab='Northing', line=2.4, cex.lab=2.8, outer=T)

dev.off()


