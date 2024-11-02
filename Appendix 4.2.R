#=================================================================================================
# Cormack-Jolly-Seber (CJS) model with data in m-array format
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2024 in Colorado
#=================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 4_Integrated population models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(6)

# Basic values
nindv <- 1800 # number of captured individuals
nyear <- 10   # number of years
ncovs <- 2    # number of environmental covariates
npcvs <- 2    # number of capture covariates

beta  <- c( 0.4,  0.1, -0.8) # intercept and slopes for apparent survival probability
alpha <- c(-1.4, -0.8,  0.4) # intercept and slopes for capture probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- colMeans(cov$size)[1:(nyear-1)] # mean habitat patch size
x2 <- colMeans(cov$temp)[1:(nyear-1)] # mean temperature
x1 <- (log(x1) - mean(log(x1))) / sd(log(x1))
x2 <- (x2 - mean(x2)) / sd(x2)
x <- cbind(x1, x2) # environmental covariates

# Simulate data
### Process
phi <- inv.logit(cbind(1,x) %*% beta) # apparent survival probability

f <- rep(1:(nyear-1), each=nindv/(nyear-1)) # year of first capture

z <- matrix(0, nindv, nyear) # individual survival history
for (v in 1:nindv) {
  z[v,f[v]] <- 1
  for (t in (f[v]+1):nyear) {
    z[v,t] <- rbinom(1, 1, z[v,t-1] * phi[t-1])
  } # t
} # v

### Capture
w <- matrix(rnorm((nyear-1)*npcvs, 0, 1), nyear-1, npcvs) # capture covariates

p <- inv.logit(cbind(1,w) %*% alpha) # capture probability

c <- matrix(0, nindv, nyear) # individual capture history
for (v in 1:nindv) {
  c[v,f[v]] <- z[v,f[v]]
  for (t in (f[v]+1):nyear) {
    c[v,t] <- rbinom(1, 1, z[v,t] * p[t-1])
  } # t
} # v

### Extend individual capture history
for (v in 1:nindv) {
  t1 <- c[v,]
  ycap <- which(t1 == 1)
  ncap <- length(ycap)
  if (ncap == 1) {
    t2 <- t1
  } else {
    t2 <- matrix(, ncap, nyear)
    for (h in 1:(ncap-1)) {
      t2[h, ycap[h]:ycap[h+1]] <- t1[ycap[h]:ycap[h+1]]
    } # h
    t2[ncap, ycap[ncap]] <- 1
  }
  if (v == 1) {
    c_ext <- t2
  } else {
    c_ext <- rbind(c_ext, t2)
  }
} # v
c_ext[which(is.na(c_ext))] <- 0

f_ext <- numeric(dim(c_ext)[1])
for (v in 1:dim(c_ext)[1]) {
  f_ext[v] <- min(which(c_ext[v,] == 1))
} # v

# Convert individual capture history to m-array
m <- matrix(0, nyear-1, nyear) # m-array of capture-recapture
for (t in 1:(nyear-1)) {
  tt <- c_ext[which(f_ext == t),]
  nband <- dim(tt)[1]
  m[t, t:(nyear-1)] <- colSums(tt, na.rm=T)[(t+1):nyear]
  m[t, nyear] <- nband - sum(m[t, 1:(nyear-1)])
} # t

#=======================
# Define MCMC algorithm
#=======================
# Define the function to calculate the cell probability of m-array
cell_prob_fun <- function(nyear, phi, p) {
  pmiss <- 1 - p
  cell_prob <- matrix(0, nyear-1, nyear)
  for (t in 1:(nyear-1)) {
    cell_prob[t,t] <- phi[t] * p[t]
  } # t
  for (t1 in 1:(nyear-2)) {
    for (t2 in (t1+1):(nyear-1)) {
      cell_prob[t1,t2] <- prod(phi[t1:t2]) * prod(pmiss[t1:(t2-1)]) * p[t2]
    } # t2
  } # t1
  for (t in 1:(nyear-1)) {
    cell_prob[t, nyear] <- 1 - sum(cell_prob[t, 1:(nyear-1)])
  } # t

  return(cell_prob)
} # cell_prob_fun

cjs_marr_mcmc <- function(m, x, w, nmcmc) {

  # Setup variables
  nyear <- dim(m)[2]
  ncovs <- dim(x)[2]
  npcvs <- dim(w)[2]

  beta_save <- matrix(0, ncovs+1, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc)
  phi_save <- matrix(0, nyear-1, nmcmc)

  # Prior
  beta_mean <- rep(0, ncovs+1)
  beta_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta <- rep(0, ncovs+1)
  alpha <- rep(0, npcvs+1)
  phi <- rep(0.5, nyear-1)
  p <- rep(0.5, nyear-1)
  cell_prob <- cell_prob_fun(nyear=nyear, phi=phi, p=p)

  # Tuning factor
  beta_tune <- 0.05
  alpha_tune <- 0.05

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    # Sample beta
    beta_star <- rnorm(ncovs+1, beta, beta_tune)
    phi_star <- inv.logit(cbind(1,x) %*% beta_star)
    cell_prob_star <- cell_prob_fun(nyear=nyear, phi=phi_star, p=p)
    prob1 <- prob2 <- numeric(nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[t] <- dmultinom(m[t,], size=sum(m[t,]), prob=cell_prob_star[t,], log=T)
      prob2[t] <- dmultinom(m[t,], size=sum(m[t,]), prob=cell_prob     [t,], log=T)
    } # t
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

    # Sample alpha
    alpha_star <- rnorm(npcvs+1, alpha, alpha_tune)
    p_star <- inv.logit(cbind(1,w) %*% alpha_star)
    cell_prob_star <- cell_prob_fun(nyear=nyear, phi=phi, p=p_star)
    prob1 <- prob2 <- numeric(nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[t] <- dmultinom(m[t,], size=sum(m[t,]), prob=cell_prob_star[t,], log=T)
      prob2[t] <- dmultinom(m[t,], size=sum(m[t,]), prob=cell_prob     [t,], log=T)
    } # t
    mh1 <- sum(prob1) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd))
    mh2 <- sum(prob2) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha <- alpha_star
      p <- p_star
      cell_prob <- cell_prob_star
    }

    # Save samples
    beta_save[,k] <- beta
    alpha_save[,k] <- alpha
    phi_save[,k] <- phi
  } # k

  list(beta_save=beta_save, alpha_save=alpha_save, phi_save=phi_save)
} # cjs_marr_mcmc

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
  cjs_marr_mcmc(m=m, x=x, w=w, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='4.1.2 cjs marr model_output.RData')

#==============
# Plot results
#==============
pdf(file='4.1.2.1 cjs marr model_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta <- c(expression(beta[0]), expression(beta["reef"]), expression(beta["temperature"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]))

par(mfrow=c(2,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,0))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta[i] - yint * 8
  ymax <- beta[i] + yint * 8
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

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i] - yint * 8
  ymax <- alpha[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
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

pdf(file='4.1.2.2 cjs marr model_phi.pdf', width=10, height=8)

library(vioplot) # for making violin plots

out_marr <- out
load('results/4.1.1 cjs indv model_output.RData')
out_indv <- output$out

phi_marr_post <- cbind(out_marr[[1]]$phi[,(nmcmc*0.2+1):nmcmc], out_marr[[2]]$phi[,(nmcmc*0.2+1):nmcmc], out_marr[[3]]$phi[,(nmcmc*0.2+1):nmcmc])
phi_indv_post <- cbind(out_indv[[1]]$phi[,(nmcmc*0.2+1):nmcmc], out_indv[[2]]$phi[,(nmcmc*0.2+1):nmcmc], out_indv[[3]]$phi[,(nmcmc*0.2+1):nmcmc])

par(mfrow=c(1,1))
par(mar=c(5,6,1,1))

plot(1, xlim=c(0.5,nyear-0.5), ylim=c(0,1), type='n', xlab='', ylab='', axes=F)
vioplot(t(phi_marr_post), col='lightcoral', rectCol='grey36', lineCol='grey36', colMed='yellow', border=NA, side='left',  add=T)
vioplot(t(phi_indv_post), col='firebrick' , rectCol='grey16', lineCol='grey16', colMed='gold'  , border=NA, side='right', add=T)
lines(phi ~ c(1:(nyear-1)), type='o', pch=16, cex=0.8, lwd=1.2, col='royalblue')
axis(1, at=seq(1,nyear,1), cex.axis=1.5)
axis(2, at=seq(0,1,0.2), cex.axis=1.5, las=2)
title(xlab='Year', cex.lab=2.5, line=3.5)
title(ylab='Apparent Survival', cex.lab=2.5, line=3.8)

points(x=1, y=0.1, pch=16, cex=1, col='royalblue')
lines(x=c(0.6,1.4), y=c(0.1,0.1), lwd=1.2, col='royalblue')
text(x=1.4, y=0.1, labels='True', pos=4, cex=1.6)
vioplot(phi_marr_post[1,]-median(phi_marr_post[1,])+0.1, at=3.4, col='lightcoral', rectCol='grey36', lineCol='grey36', colMed='yellow', border=NA, side='left',  add=T)
text(x=3.5, y=0.1, labels='Estimated from\nindividual encounter', pos=4, cex=1.6)
vioplot(phi_indv_post[1,]-median(phi_indv_post[1,])+0.1, at=6.9, col='firebrick' , rectCol='grey16', lineCol='grey16', colMed='gold'  , border=NA, side='right', add=T)
text(x=7.4, y=0.1, labels='Estimated from\nm-array', pos=4, cex=1.6)

dev.off()


