#=================================================================
# Gibbs sampling, in comparison to Metropolis-Hasting sampling
# for mean of a normal distribution with known standard deviation
# written by Qing Zhao, 2023 in Colorado
#=================================================================

#setwd('')

#===============
# Simulate data
#===============
set.seed(1)

# True values
n <- 100 # sample size
mu <- 1 # mean of y, we are going to estimate this
sigma <- 1.5 # SD of y, assumed to be known

# Simulate data
y <- rnorm(n, mu, sigma) # data

#================
# Gibbs sampling
#================
# Prior
mu0 <- 0 # mean of prior of mu
sigma0 <- 20 # SD of prior of mu

# Posterior
mu_var <- 1 / (n / sigma ^ 2 + 1 / sigma0 ^ 2) # variance of posterior of mu
mu_mean <- mu_var * (sum(y) / sigma ^ 2 + mu0 / sigma0 ^ 2) # mean of posterior of mu

# Sampling
library(foreach)    # for parallel computing
library(doParallel) # for parallel computing

numCores <- round(detectCores() / 2) # use half of the cores for parallel computing
registerDoParallel(numCores) # setup parallel computing

nmcmc <- 50000 # number of iterations
chain <- 3     # number of chains

start_time <- Sys.time() # start time of computing
mu_gibbs <- foreach (i=1:chain) %dopar% {
  rnorm(nmcmc, mu_mean, sqrt(mu_var)) # posterior samples of mu based on gibbs sampling
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

#=============================
# Metropolis-Hasting sampling
#=============================
# Define algorithm
mu_mcmc <- function(y, nmcmc) {

  mu_save <- rep(0, nmcmc) # vector to save the posterior samples of mu based on MH sampling
  mu_keep <- 0 # starting value of mu
  mu_tune <- 0.05 # tuning factor of mu

  for (k in 2:nmcmc) {
    mu_star <- rnorm(1, mu_keep, mu_tune) # proposed new value of mu
    mh1 <- sum(dnorm(y, mu_star, sigma, log=T)) + dnorm(mu_star, mu0, sigma0, log=T)
    mh2 <- sum(dnorm(y, mu_keep, sigma, log=T)) + dnorm(mu_keep, mu0, sigma0, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      mu_keep <- mu_star
    }

    mu_save[k] <- mu_keep # posterior samples of mu based on MH sampling
  } # k

  # Write output
  list(mu_save=mu_save)
} # logistic_mcmc

# Sampling
numCores <- round(detectCores() / 2) # use half of the cores for parallel computing
registerDoParallel(numCores) # setup parallel computing

nmcmc <- 50000 # number of iterations
chain <- 3     # number of chains

start_time <- Sys.time() # start time of computing
mu_mh <- foreach (i=1:chain) %dopar% {
  mu_mcmc(y=y, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

#==============
# Plot results
#==============
pdf(file='1.3 gibbs vs mh.pdf', width=10, height=8)

cols <- c('gold', 'tomato', 'royalblue')

par(mfrow=c(1,2))
par(mar=c(1,2,3,1))
par(oma=c(4,3,0,0))

plot(1, xlim=c(0,nmcmc), ylim=c(0, 2), type='n', xlab='', ylab='', axes=F)
axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.2)
axis(2, at=round(seq(0,2,0.5), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(mu_gibbs[[j]], col=cols[j])
} # j
abline(h=mu, col='grey16', lwd=1.5)
title(main='Gibbs', cex.main=2, line=0.6)

plot(1, xlim=c(0,nmcmc), ylim=c(0, 2), type='n', xlab='', ylab='', axes=F)
axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.2)
axis(2, at=round(seq(0,2,0.5), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(mu_mh[[j]][[1]], col=cols[j])
} # j
abline(h=mu, col='grey16', lwd=1.5)
title(main='MH', cex.main=2, line=0.6)

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=1.1, outer=T)

dev.off()


