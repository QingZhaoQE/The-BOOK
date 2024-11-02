#=================================================================================================
# Age-structured Dail-Madsen model
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
nsite <- 360 # number of sites
nyear <- 20  # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 1   # number of observational covariates

beta0_adu <- c(1.8, 0.5, -0.4)          # intercept and slopes for adult initial abundance
beta0_juv <- c(1.4,   1, -0.8)          # intercept and slopes for juvenile initial abundance
beta_phi_adu <- c(0.8, -0.3, 0.4, -0.2) # intercept and slopes for adult apparent survival
beta_phi_juv <- c(0.4, -0.5, 0.6, -0.3) # intercept and slopes for juvenile apparent survival
beta_omega   <- c(0.5, -0.1, 0.5, -0.3) # intercept and slopes for transition probability from juvenile to adult
beta_gamma   <- c(0.4, -0.6, 0.6, -0.3) # intercept and slopes for gain (recruitment)
alpha <- c(0.3, -0.6)                   # intercept and slopes for detection

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
lambda_adu <- matrix(, nsite, nyear) # expectation of adult abundance
lambda_juv <- matrix(, nsite, nyear) # expectation of juvenile abundance
N_adu <- matrix(, nsite, nyear) # adult abundance
N_juv <- matrix(, nsite, nyear) # juvenile abundance
N <- matrix(, nsite, nyear) # total abundance

lambda_adu[,1] <- exp(cbind(1,x[,1,]) %*% beta0_adu) # expectation of adult initial abundance
lambda_juv[,1] <- exp(cbind(1,x[,1,]) %*% beta0_juv) # expectation of juvenile initial abundance
lambda0 <- lambda_adu[,1] + lambda_juv[,1] # expectation of total initial abundance
N_adu[,1] <- rpois(nsite, lambda_adu[,1]) # adult initial abundance
N_juv[,1] <- rpois(nsite, lambda_juv[,1]) # juvenile initial abundance
N[,1] <- N_adu[,1] + N_juv[,1] # total initial abundance

phi_adu <- matrix(, nsite, nyear-1) # adult apparent survival
phi_juv <- matrix(, nsite, nyear-1) # juvenile apparent survival
omega   <- matrix(, nsite, nyear-1) # transition from juvenile to adult
gamma   <- matrix(, nsite, nyear-1) # recrutiment
for (t in 2:nyear) {
  phi_adu[,t-1] <- inv.logit(cbind(1, (N[,t-1]-lambda0)/lambda0, x[,t,]) %*% beta_phi_adu)
  phi_juv[,t-1] <- inv.logit(cbind(1, (N[,t-1]-lambda0)/lambda0, x[,t,]) %*% beta_phi_juv)
  omega[,t-1] <- inv.logit(cbind(1, (N[,t-1]-lambda0)/lambda0, x[,t,]) %*% beta_omega)
  gamma[,t-1] <- exp      (cbind(1, (N[,t-1]-lambda0)/lambda0, x[,t,]) %*% beta_gamma)
  lambda_adu[,t] <- N_adu[,t-1] * phi_adu[,t-1] + N_juv[,t-1] * phi_juv[,t-1] * omega[,t-1]
  lambda_juv[,t] <- N_juv[,t-1] * phi_juv[,t-1] * (1 - omega[,t-1]) + gamma[,t-1]
  N_adu[,t] <- rpois(nsite, lambda_adu[,t])
  N_juv[,t] <- rpois(nsite, lambda_juv[,t])
  N[,t] <- N_adu[,t] + N_juv[,t]
} # t

w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

p <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    p[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha)
  } # t
} # i

y_adu <- array(, dim=c(nsite, nyear, nreps)) # count data of adult
y_juv <- array(, dim=c(nsite, nyear, nreps)) # count data of juvenile
for (i in 1:nsite) {
  for (t in 1:nyear) {
    y_adu[i,t,] <- rbinom(nreps, N_adu[i,t], p[i,t,])
    y_juv[i,t,] <- rbinom(nreps, N_juv[i,t], p[i,t,])
  } # t
} # i

#=======================
# Define MCMC algorithm
#=======================
age_Dail_Madsen_mcmc <- function(y_adu, y_juv, x, w, nmcmc) {

  # Setup variables
  nsite <- dim(y_adu)[1]
  nyear <- dim(y_adu)[2]
  nreps <- dim(y_adu)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  y_adu_max <- apply(y_adu, 1:2, max)
  y_juv_max <- apply(y_juv, 1:2, max)

  beta0_adu_save <- matrix(0, ncovs+1, nmcmc) 
  beta0_juv_save <- matrix(0, ncovs+1, nmcmc) 
  beta_phi_adu_save <- matrix(0, ncovs+2, nmcmc)
  beta_phi_juv_save <- matrix(0, ncovs+2, nmcmc)
  beta_omega_save   <- matrix(0, ncovs+2, nmcmc)
  beta_gamma_save   <- matrix(0, ncovs+2, nmcmc)
  alpha_save <- matrix(0, npcvs+1, nmcmc) 
  N_adu_save <- matrix(0, nyear, nmcmc)
  N_juv_save <- matrix(0, nyear, nmcmc)
  phi_adu_save <- matrix(0, nyear-1, nmcmc)
  phi_juv_save <- matrix(0, nyear-1, nmcmc)
  omega_save   <- matrix(0, nyear-1, nmcmc)
  gamma_save   <- matrix(0, nyear-1, nmcmc)

  # Priors
  beta0_adu_mean <- rep(0, ncovs+1)
  beta0_adu_sd <- 2
  beta0_juv_mean <- rep(0, ncovs+1)
  beta0_juv_sd <- 2
  beta_phi_adu_mean <- rep(0, ncovs+2)
  beta_phi_adu_sd <- 2
  beta_phi_juv_mean <- rep(0, ncovs+2)
  beta_phi_juv_sd <- 2
  beta_omega_mean <- rep(0, ncovs+2)
  beta_omega_sd <- 2
  beta_gamma_mean <- rep(0, ncovs+2)
  beta_gamma_sd <- 2
  alpha_mean <- rep(0, npcvs+1)
  alpha_sd <- 2

  # Starting values
  beta0_adu <- rep(0, ncovs+1)
  beta0_juv <- rep(0, ncovs+1)
  beta_phi_adu <- rep(0, ncovs+2)
  beta_phi_juv <- rep(0, ncovs+2)
  beta_omega <- rep(0, ncovs+2)
  beta_gamma <- rep(0, ncovs+2)
  alpha <- rep(0, npcvs+1)
  lambda0_adu <- exp(cbind(1,x[,1,]) %*% beta0_adu)
  lambda0_juv <- exp(cbind(1,x[,1,]) %*% beta0_juv)
  lambda0 <- lambda0_adu + lambda0_juv
  phi_adu <- matrix(0.5, nsite, nyear-1)
  phi_juv <- matrix(0.5, nsite, nyear-1)
  omega <- matrix(0.5, nsite, nyear-1)
  gamma <- matrix(1, nsite, nyear-1)
  N_adu <- round((y_adu_max + 1) / 0.5)
  N_juv <- round((y_juv_max + 1) / 0.5)
  N <- N_adu + N_juv
  p <- array(0.5, dim=c(nsite, nyear, nreps))

  # Tuning factors
  beta0_adu_tune <- 0.05
  beta0_juv_tune <- 0.05
  beta_phi_adu_tune <- c(0.12, 0.100, 0.100, 0.100)
  beta_phi_juv_tune <- c(0.15, 0.100, 0.100, 0.100)
  beta_omega_tune <- c(0.08, 0.050, 0.050, 0.050)
  beta_gamma_tune <- c(0.05, 0.035, 0.035, 0.035)
  alpha_tune <- 0.01
  N_adu_tune <- 1
  N_juv_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N_adu
    N_adu_star <- matrix(rpois(nsite*nyear, N_adu + N_adu_tune), nsite, nyear)
    mh1 <- apply(dbinom(y_adu, N_adu_star, p, log=T), 1:2, sum) + 
           dpois(N_adu_star, cbind(lambda0_adu, N_adu[,-nyear] * phi_adu + N_juv[,-nyear] * phi_juv * omega  + 1e-6), log=T) + 
           dpois(N_adu, N_adu_star + N_adu_tune, log=T)
    mh2 <- apply(dbinom(y_adu, N_adu     , p, log=T), 1:2, sum) + 
           dpois(N_adu     , cbind(lambda0_adu, N_adu[,-nyear] * phi_adu + N_juv[,-nyear] * phi_juv * omega + 1e-6), log=T) + 
           dpois(N_adu_star, N_adu + N_adu_tune, log=T)
    mhNa <- exp(mh1 - mh2)
    Nakeep <- ((mhNa > runif(nsite*nyear)) & (N_adu_star >= y_adu_max))
    N_adu[Nakeep] <- N_adu_star[Nakeep]
    N <- N_adu + N_juv

    ### Sample N_juv
    N_juv_star <- matrix(rpois(nsite*nyear, N_juv + N_juv_tune), nsite, nyear)
    mh1 <- apply(dbinom(y_juv, N_juv_star, p, log=T), 1:2, sum) + 
           dpois(N_juv_star, cbind(lambda0_juv, N_juv[,-nyear] * phi_juv * (1 - omega) + gamma), log=T) + 
           dpois(N_juv, N_juv_star + N_juv_tune, log=T)
    mh2 <- apply(dbinom(y_juv, N_juv     , p, log=T), 1:2, sum) + 
           dpois(N_juv     , cbind(lambda0_juv, N_juv[,-nyear] * phi_juv * (1 - omega) + gamma), log=T) + 
           dpois(N_juv_star, N_juv + N_juv_tune, log=T)
    mhNj <- exp(mh1 - mh2)
    Njkeep <- ((mhNj > runif(nsite*nyear)) & (N_juv_star >= y_juv_max))
    N_juv[Njkeep] <- N_juv_star[Njkeep]
    N <- N_adu + N_juv

    ### Sample beta0_adu
    beta0_adu_star <- rnorm(ncovs+1, beta0_adu, beta0_adu_tune)
    lambda0_adu_star <- exp(cbind(1,x[,1,]) %*% beta0_adu_star)
    mh1 <- sum(dpois(N_adu[,1], lambda0_adu_star, log=T)) + 
           sum(dnorm(beta0_adu_star, beta0_adu_mean, beta0_adu_sd, log=T))
    mh2 <- sum(dpois(N_adu[,1], lambda0_adu     , log=T)) + 
           sum(dnorm(beta0_adu     , beta0_adu_mean, beta0_adu_sd, log=T))
    mh0a <- exp(mh1 - mh2)
    if (mh0a > runif(1)) {
      beta0_adu <- beta0_adu_star
      lambda0_adu<- lambda0_adu_star
    }

    ### Sample beta0_juv
    beta0_juv_star <- rnorm(ncovs+1, beta0_juv, beta0_juv_tune)
    lambda0_juv_star <- exp(cbind(1,x[,1,]) %*% beta0_juv_star)
    mh1 <- sum(dpois(N_juv[,1], lambda0_juv_star, log=T)) + 
           sum(dnorm(beta0_juv_star, beta0_juv_mean, beta0_juv_sd, log=T))
    mh2 <- sum(dpois(N_juv[,1], lambda0_juv     , log=T)) + 
           sum(dnorm(beta0_juv     , beta0_juv_mean, beta0_juv_sd, log=T))
    mh0j <- exp(mh1 - mh2)
    if (mh0j > runif(1)) {
      beta0_juv <- beta0_juv_star
      lambda0_juv <- lambda0_juv_star
    }
    lambda0 <- lambda0_adu + lambda0_juv

    ### Sample beta_phi_adu
    beta_phi_adu_star <- rnorm(ncovs+2, beta_phi_adu, beta_phi_adu_tune)
    phi_adu <- phi_adu_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi_adu     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_adu     )
      phi_adu_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_adu_star)
    } # t
    mh1 <- sum(dpois(N_adu[,-1], N_adu[,-nyear] * phi_adu_star + N_juv[,-nyear] * phi_juv * omega + 1e-6, log=T)) + 
           sum(dnorm(beta_phi_adu_star, beta_phi_adu_mean, beta_phi_adu_sd, log=T))
    mh2 <- sum(dpois(N_adu[,-1], N_adu[,-nyear] * phi_adu      + N_juv[,-nyear] * phi_juv * omega + 1e-6, log=T)) + 
           sum(dnorm(beta_phi_adu     , beta_phi_adu_mean, beta_phi_adu_sd, log=T))
    mhSa <- exp(mh1 - mh2)
    if (mhSa > runif(1)) {
      beta_phi_adu <- beta_phi_adu_star
    }

    ### Sample beta_phi_juv
    beta_phi_juv_star <- rnorm(ncovs+2, beta_phi_juv, beta_phi_juv_tune)
    phi_juv <- phi_juv_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi_juv     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_juv     )
      phi_juv_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_juv_star)
    } # t
    mh1 <- sum(dpois(N_adu[,-1], N_adu[,-nyear] * phi_adu + N_juv[,-nyear] * phi_juv_star * omega + 1e-6, log=T)) + 
           sum(dpois(N_juv[,-1], N_juv[,-nyear] * phi_juv_star * (1 - omega) + gamma, log=T)) + 
           sum(dnorm(beta_phi_juv_star, beta_phi_juv_mean, beta_phi_juv_sd, log=T))
    mh2 <- sum(dpois(N_adu[,-1], N_adu[,-nyear] * phi_adu + N_juv[,-nyear] * phi_juv      * omega + 1e-6, log=T)) + 
           sum(dpois(N_juv[,-1], N_juv[,-nyear] * phi_juv      * (1 - omega) + gamma, log=T)) + 
           sum(dnorm(beta_phi_juv     , beta_phi_juv_mean, beta_phi_juv_sd, log=T))
    mhSj <- exp(mh1 - mh2)
    if (mhSj > runif(1)) {
      beta_phi_juv <- beta_phi_juv_star
    }

    ### Sample beta_omega
    beta_omega_star <- rnorm(ncovs+2, beta_omega, beta_omega_tune)
    omega <- omega_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      omega     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_omega     )
      omega_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_omega_star)
    } # t
    mh1 <- sum(dpois(N_adu[,-1], N_adu[,-nyear] * phi_adu + N_juv[,-nyear] * phi_juv * omega_star + 1e-6, log=T)) + 
           sum(dpois(N_juv[,-1], N_juv[,-nyear] * phi_juv * (1 - omega_star) + gamma, log=T)) + 
           sum(dnorm(beta_omega_star, beta_omega_mean, beta_omega_sd, log=T))
    mh2 <- sum(dpois(N_adu[,-1], N_adu[,-nyear] * phi_adu + N_juv[,-nyear] * phi_juv * omega      + 1e-6, log=T)) + 
           sum(dpois(N_juv[,-1], N_juv[,-nyear] * phi_juv * (1 - omega     ) + gamma, log=T)) + 
           sum(dnorm(beta_omega     , beta_omega_mean, beta_omega_sd, log=T))
    mhD <- exp(mh1 - mh2)
    if (mhD > runif(1)) {
      beta_omega <- beta_omega_star
    }

    ### Sample beta_gamma
    beta_gamma_star <- rnorm(ncovs+2, beta_gamma, beta_gamma_tune)
    gamma <- gamma_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      gamma     [,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma     )
      gamma_star[,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma_star)
    } # t
    mh1 <- sum(dpois(N_juv[,-1], N_juv[,-nyear] * phi_juv * (1 - omega) + gamma_star, log=T)) + 
           sum(dnorm(beta_gamma_star, beta_gamma_mean, beta_gamma_sd, log=T))
    mh2 <- sum(dpois(N_juv[,-1], N_juv[,-nyear] * phi_juv * (1 - omega) + gamma     , log=T)) + 
           sum(dnorm(beta_gamma     , beta_gamma_mean, beta_gamma_sd, log=T))
    mhG <- exp(mh1 - mh2)
    if (mhG > runif(1)) {
      beta_gamma <- beta_gamma_star
    }

    ### Sample alpha
    alpha_star <- rnorm(npcvs+1, alpha, alpha_tune)
    p_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        p_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_star)
      } # t
    } # i
    mh1 <- sum(dbinom(y_adu, N_adu, p_star, log=TRUE), na.rm=T) + 
           sum(dbinom(y_juv, N_juv, p_star, log=TRUE), na.rm=T) + 
           sum(dnorm(alpha_star, alpha_mean, alpha_sd, log=TRUE))
    mh2 <- sum(dbinom(y_adu, N_adu, p     , log=TRUE), na.rm=T) + 
           sum(dbinom(y_juv, N_juv, p     , log=TRUE), na.rm=T) + 
           sum(dnorm(alpha     , alpha_mean, alpha_sd, log=TRUE))
    mhA <- exp(mh1 - mh2)
    if (mhA > runif(1)) {
      alpha <- alpha_star
      p <- p_star
    }

    ### Save samples
    beta0_adu_save[,k] <- beta0_adu
    beta0_juv_save[,k] <- beta0_juv
    beta_phi_adu_save[,k] <- beta_phi_adu
    beta_phi_juv_save[,k] <- beta_phi_juv
    beta_omega_save[,k] <- beta_omega
    beta_gamma_save[,k] <- beta_gamma
    alpha_save[,k] <- alpha
    N_adu_save[,k] <- colMeans(N_adu)
    N_juv_save[,k] <- colMeans(N_juv)
    phi_adu_save[,k] <- colMeans(phi_adu)
    phi_juv_save[,k] <- colMeans(phi_juv)
    omega_save[,k] <- colMeans(omega)
    gamma_save[,k] <- colMeans(gamma)
  } # k

  # Write output
  list(beta0_adu_save=beta0_adu_save, beta0_juv_save=beta0_juv_save, 
       beta_phi_adu_save=beta_phi_adu_save, beta_phi_juv_save=beta_phi_juv_save, 
       beta_omega_save=beta_omega_save, beta_gamma_save=beta_gamma_save, 
       alpha_save=alpha_save, 
       N_adu_save=N_adu_save, N_juv_save=N_juv_save, 
       phi_adu_save=phi_adu_save, phi_juv_save=phi_juv_save, 
       omega_save=omega_save, gamma_save=gamma_save)

} # age_Dail_Madsen_mcmc

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
  age_Dail_Madsen_mcmc(y_adu=y_adu, y_juv=y_juv, x=x, w=w, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='3.5.2 age Dail-Madsen model_output.RData')

#==============
# Plot results
#==============
pdf(file='3.5.2.1 age Dail-Madsen model_chains.pdf', width=10, height=10)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0_adu <- c(expression(beta["0,adu"]^"[0]"), expression(beta["pond,adu"]^"[0]"), expression(beta["temperature,adu"]^"[0]"))
ylab_beta0_juv <- c(expression(beta["0,juv"]^"[0]"), expression(beta["pond,juv"]^"[0]"), expression(beta["temperature,juv"]^"[0]"))
ylab_beta_phi_adu <- c(expression(beta["0,adu"]^""["["*phi*"]"]), expression(beta["density,adu"]^""["["*phi*"]"]), 
                     expression(beta["pond,adu"]^""["["*phi*"]"]), expression(beta["temperature,adu"]^""["["*phi*"]"]))
ylab_beta_phi_juv <- c(expression(beta["0,juv"]^""["["*phi*"]"]), expression(beta["density,juv"]^""["["*phi*"]"]), 
                     expression(beta["pond,juv"]^""["["*phi*"]"]), expression(beta["temperature,juv"]^""["["*phi*"]"]))
ylab_beta_omega <- c(expression(beta[0]^""["["*omega*"]"]), expression(beta["density"]^""["["*omega*"]"]), 
                     expression(beta["pond"]^""["["*omega*"]"]), expression(beta["temperature"]^""["["*omega*"]"]))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta["density"]^""["["*gamma*"]"]), 
                     expression(beta["pond"]^""["["*gamma*"]"]), expression(beta["temperature"]^""["["*gamma*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]), 
                expression(alpha[3]), expression(alpha[4]))

par(mfrow=c(8,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,1))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta0_adu_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0_adu[i] - yint * 4
  ymax <- beta0_adu[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0_adu[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0_adu[i], cex.main=2, line=1.4)
  text(x=nmcmc*0.42, y=beta0_adu[i]+yint*3, pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta0_juv_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0_juv[i] - yint * c(4,4,8)[i]
  ymax <- beta0_juv[i] + yint * c(4,4,8)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(2,2,4)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0_juv[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0_juv[i], cex.main=2, line=1.4)
  text(x=nmcmc*0.42, y=beta0_juv[i]+yint*c(3,3,6)[i], pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_adu_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi_adu[i] - yint * 4
  ymax <- beta_phi_adu[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi_adu[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi_adu[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_phi_adu[i]+yint*3, pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_juv_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi_juv[i] - yint * 8
  ymax <- beta_phi_juv[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi_juv[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi_juv[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_phi_juv[i]+yint*6, pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_omega_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_omega[i] - yint * 8
  ymax <- beta_omega[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_omega[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_omega[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_omega[i]+yint*6, pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_gamma_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma[i] - yint * c(4,4,2,2)[i]
  ymax <- beta_gamma[i] + yint * c(4,4,2,2)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  if (i < (ncovs+2)) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*c(2,2,1,1)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_gamma[i]+yint*c(3,3,1.5,1.5)[i], pos=4, cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

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

pdf(file='3.5.2.2 age Dail-Madsen model_N phi omega gamma.pdf', width=8, height=8)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    N_adu_mean <- out[[1]]$N_adu_save
    N_juv_mean <- out[[1]]$N_juv_save
    phi_adu_mean <- out[[1]]$phi_adu_save
    phi_juv_mean <- out[[1]]$phi_juv_save
    omega_mean <- out[[1]]$omega_save
    gamma_mean <- out[[1]]$gamma_save
  } else {
    N_adu_mean <- cbind(N_adu_mean, out[[1]]$N_adu_save)
    N_juv_mean <- cbind(N_juv_mean, out[[1]]$N_juv_save)
    phi_adu_mean <- cbind(phi_adu_mean, out[[1]]$phi_adu_save)
    phi_juv_mean <- cbind(phi_juv_mean, out[[1]]$phi_juv_save)
    omega_mean <- cbind(omega_mean, out[[1]]$omega_save)
    gamma_mean <- cbind(gamma_mean, out[[1]]$gamma_save)
  }
} # i

par(mfrow=c(3,2))
par(mar=c(1.5,6.5,0,0))
par(oma=c(4,0,1,1))

plot(1, xlim=c(1,nyear), ylim=c(0,10), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_adu_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(N_adu), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), labels=rep('',5))
axis(2, at=seq(0,10,2), cex.axis=1.6, las=2)
axis(2, at=5, labels='Adult Abundance', cex.axis=2.2, line=2.6, tick=F)

points(x=4, y=9, pch=16, cex=1, col='royalblue')
lines(x=c(3.5,4.5), y=c(9,9), lwd=1.2, col='royalblue')
vioplot(rnorm(10000,9,0.25), at=10.5, add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
text(x=4.5, y=9, labels='True', pos=4, cex=2.2)
text(x=11, y=9, labels='Estimated', pos=4, cex=2.2)

plot(1, xlim=c(1,nyear), ylim=c(0,10), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_juv_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(N_juv), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), labels=rep('',5))
axis(2, at=seq(0,10,2), cex.axis=1.6, las=2)
axis(2, at=5, labels='Juvenile Abundance', cex.axis=2.2, line=2.6, tick=F)

plot(1, xlim=c(1,nyear), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_adu_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(phi_adu), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), labels=rep('',5))
axis(2, at=seq(0,1,0.2), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Adult Survival', cex.axis=2.2, line=2.6, tick=F)

plot(1, xlim=c(1,nyear), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_juv_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(phi_juv), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), labels=rep('',5))
axis(2, at=seq(0,1,0.2), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Juvenile Survival', cex.axis=2.2, line=2.6, tick=F)

plot(1, xlim=c(1,nyear), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(omega_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(omega), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), cex.axis=1.6)
axis(2, at=seq(0,1,0.2), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Transition', cex.axis=2.2, line=2.6, tick=F)

plot(1, xlim=c(1,nyear), ylim=c(0,5), type='l', axes=F, xlab='', ylab='')
vioplot(t(gamma_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(gamma), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,5), cex.axis=1.6)
axis(2, at=seq(0,5,1), cex.axis=1.6, las=2)
axis(2, at=2.5, labels='Recruitment', cex.axis=2.2, line=2.6, tick=F)

title(xlab='Year', cex.lab=3, line=2.4, outer=T)

dev.off()


