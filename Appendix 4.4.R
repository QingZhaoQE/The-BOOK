#=================================================================================================
# Spatial integrated population model (SIPM) without age-specific survival and movement
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
nsite <- 60 # number of sites
nyear <- 10 # number of years
nreps <- 5  # number of within-season replicates
ncovs <- 2  # number of environmental covariates
npcvs <- 2  # number of observational covariates

beta0 <- c(3, 1, -0.8)                 # intercept and slopes for initial abundance (i.e., number of adults)
beta_gamma <- c(-0.8, -0.3, 0.4, -0.2) # intercept and slopes for reproduction rate
beta_phi   <- c( 0.5, -0.5, 0.6, -0.3) # intercept and slopes for survival probability
kappa <- 0.5                           # distance effect on movement rate
alpha_det <- c(0.9, -0.4, 0.2)         # intercept and slopes for detection probability
alpha_cap <- c(-3, -0.8, 0.4)          # intercept and slopes for capture probability

# Read in environmental covariates
load('data/covariates.RData')
x1 <- (log(cov$size) - mean(log(cov$size))) / sd(log(cov$size)) # standardized log habitat patch size
x2 <- (cov$temp - mean(cov$temp)) / sd(cov$temp) # standardized temperature
x <- array(, dim=c(nsite, nyear, ncovs)) # environmental covariates
x[,,1] <- x1[1:nsite, 1:nyear]
x[,,2] <- x2[1:nsite, 1:nyear]
### Standardize covariates for each site
for (i in 1:nsite) {
  for (h in 1:ncovs) {
    x[i,2:nyear,h] <- (x[i,2:nyear,h] - mean(x[i,,h])) / sd(x[i,,h])
  } # h
} # i

# Simulate data
N <- matrix(, nsite, nyear) # abundance
lambda0 <- exp(cbind(1,x[,1,]) %*% beta0) # expectation of initial abundance
N[,1] <- rpois(nsite, lambda0) # initial abundance

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

nindv <- 15000 # maximum number of individuals for simulation
z <- matrix(0, nindv, nyear) # individual survival-location history
for (i in 1:nsite) {
  if (i == 1) { tt <- rep(i, N[i,1]) } else { tt <- c(tt, rep(i, N[i,1])) }
}
z[1:length(tt),1] <- tt

gamma <- matrix(, nsite, nyear-1) # reproduction rate
phi   <- matrix(, nsite, nyear-1) # survival probability
R     <- matrix(, nsite, nyear-1) # number of offspring

for (t in 2:nyear) {
  gamma[,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma)
  R[,t-1] <- rpois(nsite, N[,t-1] * gamma[,t-1])

  for (i in 1:nsite) {
    if (i == 1) { tt <- rep(i, R[i,t-1]) } else { tt <- c(tt, rep(i, R[i,t-1])) }
  } # i
  z[max(which(rowMeans(z,na.rm=T)>0))+c(1:length(tt)), t-1] <- tt

  phi[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi)
  transition_matrix <- matrix(, nsite, nsite+1) # matrix of transition (survival & movement) probability
  for (i in 1:nsite) {
    transition_matrix[i,1:nsite] <- phi[i,t-1] * theta[i,]
    transition_matrix[i,nsite+1] <- 1 - phi[i,t-1]
  } # i

  for (v in 1:nindv) {
    if (z[v,t-1] %in% 1:nsite) {
      z[v,t] <- rcat(1, p=transition_matrix[z[v,t-1],])
    }
  } # v
  z[which(!(z %in% 1:nsite))] <- 0

  for (i in 1:nsite) {
    N[i,t] <- length(which(z[,t] == i))
  }
} # t
z[which(!(z %in% 1:nsite))] <- 0
z <- z[which(rowMeans(z, na.rm=T) != 0),] # remove individuals that never exist

# Simulate count data
w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

pdet <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    pdet[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_det)
  } # t
} # i

y <- array(, dim=c(nsite, nyear, nreps)) # count data of breeding population
for (i in 1:nsite) {
  for (t in 1:nyear) {
    y[i,t,] <- rbinom(nreps, N[i,t], pdet[i,t,])
  } # t
} # i

r <- array(, dim=c(nsite, nyear-1, nreps)) # count data of offspring
for (i in 1:nsite) {
  for (t in 1:(nyear-1)) {
    r[i,t,] <- rbinom(nreps, R[i,t], pdet[i,t,])
  } # t
} # i

# simulate capture-recapture data
pcap_per_visit <- array(, dim=c(nsite, nyear, nreps)) # per visit capture probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    pcap_per_visit[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_cap)
  } # t
} # i

c_all <- matrix(, dim(z)[1], nyear) # capture history of all individuals, 
                                    # including individuals that have never been captured
for (v in 1:dim(z)[1]) {
  for (t in 1:nyear) {
    if (z[v,t] == 0) {
      c_all[v,t] <- 0
    } else {
      cap_per_visit <- rbinom(nreps, 1, pcap_per_visit[z[v,t],t,])
      c01 <- ifelse(sum(cap_per_visit)==0, 0, 1)
      c_all[v,t] <- c01 * z[v,t]
    }
  } # t
} # v
if_cap <- ifelse(rowSums(c_all) > 0, 1, 0)

c <- c_all[which(if_cap==1),] # capture history of individuals that have ever been captured
f <- numeric(dim(c)[1]) # year of first capture
for (v in 1:length(f)) {
  f[v] <- min(which(c[v,] %in% c(1:nsite)))
} # v
c <- c[which(f %in% c(1:(nyear-1))),]

# Extend individual capture history
for (v in 1:dim(c)[1]) {
  t1 <- c[v,]
  ycap <- which(t1 %in% c(1:nsite))
  ncap <- length(ycap)
  if (ncap == 1) {
    t2 <- t1
  } else {
    t2 <- matrix(, ncap, nyear)
    for (k in 1:(ncap-1)) {
      t2[k, ycap[k]:ycap[k+1]] <- t1[ycap[k]:ycap[k+1]]
    } # k
    t2[ncap, ycap[ncap]] <- t1[ycap[ncap]]
  }
  if (v == 1) {
    c_ext <- t2
  } else {
    c_ext <- rbind(c_ext, t2)
  }
} # v
c_ext[which(is.na(c_ext))] <- 0

f_ext <- numeric(dim(c_ext)[1]) # year of first capture for the extended individual capture history
for (v in 1:length(f_ext)) {
  f_ext[v] <- min(which(c_ext[v,] %in% c(1:nsite)))
} # i

# Converting individual capture history to m-array format
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
  temp1 <- c_ext[which(f_ext == t1),]
  for (i1 in 1:nsite) {
    temp2 <- temp1[which(temp1[,t1] == i1),]
    nband <- length(which(temp1[,t1] == i1))
    if (nband == 1) {
      temp2 <- matrix(temp2, nrow=1, ncol=length(temp2))
    }
    for (t2 in (t1+1):nyear) {
      for (i2 in 1:nsite) {
        m[t1,(t2-2)*nsite+i2,i1] <- length(which(temp2[,t2] == i2))
      } # j
    } # s
    m[t1,nsite*(nyear-1)+1,i1] <- nband - sum(m[t1,1:(nsite*(nyear-1)),i1])
  } # i1
} # t1

#=======================
# Define MCMC algorithm
#=======================
# Function to calculate the cell probability of m-array
cell_prob_fun <- function(nsite, nyear, phi, theta, pcap) {
  qcap <- 1 - pcap
  cell_prob_p <- cell_prob_q <- array(0, dim=c(nyear-1, nsite*(nyear-1)+1, nsite))
  for (i in 1:nsite) {
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,nsite*(t-1)+1:nsite,i] <- phi[i,t] * theta[i,] * pcap[,t]
      cell_prob_q[t,nsite*(t-1)+1:nsite,i] <- phi[i,t] * theta[i,] * qcap[,t]
    } # t
    for (t1 in 1:(nyear-2)) {
      for (t2 in (t1+1):(nyear-1)) {
        cell_prob_p[t1,nsite*(t2-1)+1:nsite,i] <- ((cell_prob_q[t1,nsite*(t2-2)+1:nsite,i] * phi[,t2]) %*% theta) * pcap[,t2]
        cell_prob_q[t1,nsite*(t2-1)+1:nsite,i] <- ((cell_prob_q[t1,nsite*(t2-2)+1:nsite,i] * phi[,t2]) %*% theta) * qcap[,t2]
      } # t2
    } # t1
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,nsite*(nyear-1)+1,i] <- 1 - sum(cell_prob_p[t,1:(nsite*(nyear-1)),i])
    } # t
  } # i
  cell_prob <- cell_prob_p

  return(cell_prob)
} # cell_prob_fun

sipm_mcmc <- function(y, r, m, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:2, max)
  site_large <- max(which(rowMeans(x[,,1]) == max(rowMeans(x[,,1]))))
  site_small <- min(which(rowMeans(x[,,1]) == min(rowMeans(x[,,1]))))

  beta0_save <- matrix(0, ncovs+1, nmcmc)
  beta_gamma_save <- matrix(0, ncovs+2, nmcmc)
  beta_phi_save <- matrix(0, ncovs+2, nmcmc)
  kappa_save <- rep(0, nmcmc)
  alpha_det_save <- matrix(0, npcvs+1, nmcmc)
  alpha_cap_save <- matrix(0, npcvs+1, nmcmc)
  N_mean_save <- matrix(0, nyear, nmcmc)
  N_large_save <- matrix(0, nyear, nmcmc)
  N_small_save <- matrix(0, nyear, nmcmc)
  phi_mean_save <- matrix(0, nyear-1, nmcmc)
  phi_large_save <- matrix(0, nyear-1, nmcmc)
  phi_small_save <- matrix(0, nyear-1, nmcmc)
  gamma_mean_save <- matrix(0, nyear-1, nmcmc)
  gamma_large_save <- matrix(0, nyear-1, nmcmc)
  gamma_small_save <- matrix(0, nyear-1, nmcmc)

  # Prior
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_gamma_mean <- rep(0, ncovs+2)
  beta_gamma_sd <- 2
  beta_phi_mean <- rep(0, ncovs+2)
  beta_phi_sd <- 2
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  alpha_det_mean <- rep(0, npcvs+1)
  alpha_det_sd <- 2
  alpha_cap_mean <- rep(0, npcvs+1)
  alpha_cap_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_gamma <- rep(0, ncovs+2)
  beta_phi   <- rep(0, ncovs+2)
  kappa <- 1
  alpha_det <- rep(0, npcvs+1)
  alpha_cap <- rep(0, npcvs+1)
  lambda0 <- exp(cbind(1,x[,1,]) %*% beta0)
  gamma <- matrix(0.5, nsite, nyear-1)
  phi   <- matrix(0.5, nsite, nyear-1)
  eta <- exp(-1 * kappa * d)
  theta <- eta / rowSums(eta)
  N <- round((ymax + 1) / .5)
  pdet <- array(.5, dim=c(nsite, nyear, nreps))
  pcap_per_visit <- array(.5, dim=c(nsite, nyear-1, nreps))
  pcap <- matrix(, nsite, nyear-1)
  for (i in 1:nsite) {
    for (t in 1:(nyear-1)) {
      pcap[i,t] <- 1 - prod(1 - pcap_per_visit[i,t,])
    } # t
  } # i
  cell_prob <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi, theta=theta, pcap=pcap)

  # Tuning factor
  beta0_tune <- 0.08
  beta_gamma_tune <- 0.02
  beta_phi_tune <- 0.08
  kappa_tune <- 0.035
  alpha_det_tune <- 0.025
  alpha_cap_tune <- 0.05
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- matrix(rpois(nsite*nyear, N + N_tune), nsite, nyear)
    mh1 <- apply(dbinom(y, N_star, pdet, log=T), 1:2, sum) + 
           dpois(N_star, cbind(lambda0, t(t(N[,-nyear] * (1 + gamma) * phi) %*% theta)), log=T) + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , pdet, log=T), 1:2, sum) + 
           dpois(N     , cbind(lambda0, t(t(N[,-nyear] * (1 + gamma) * phi) %*% theta)), log=T) + 
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

    # Sample beta_gamma
    beta_gamma_star <- rnorm(ncovs+2, beta_gamma, beta_gamma_tune)
    gamma <- gamma_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      gamma     [,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma     )
      gamma_star[,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma_star)
    } # t
    mh1 <- sum(dbinom(r, r + y[,-nyear,], gamma_star / (gamma_star + 1), log=T)) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * (1 + gamma_star) * phi) %*% theta), log=T)) + 
           sum(dnorm(beta_gamma_star, beta_gamma_mean, beta_gamma_sd, log=T))
    mh2 <- sum(dbinom(r, r + y[,-nyear,], gamma      / (gamma      + 1), log=T)) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * (1 + gamma     ) * phi) %*% theta), log=T)) + 
           sum(dnorm(beta_gamma     , beta_gamma_mean, beta_gamma_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_gamma <- beta_gamma_star
    }

    # Sample beta_phi
    beta_phi_star <- rnorm(ncovs+2, beta_phi, beta_phi_tune)
    phi <- phi_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi)
      phi_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_star)
    } # t
    cell_prob      <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi     , theta=theta, pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi_star, theta=theta, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * (1 + gamma) * phi_star) %*% theta), log=T)) + 
           sum(dnorm(beta_phi_star, beta_phi_mean, beta_phi_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * (1 + gamma) * phi     ) %*% theta), log=T)) + 
           sum(dnorm(beta_phi     , beta_phi_mean, beta_phi_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_phi <- beta_phi_star
      phi <- phi_star
      cell_prob <- cell_prob_star
    }

    # Sample kappa
    kappa_star <- exp(rnorm(1, log(kappa), kappa_tune))
    eta_star <- exp(-1 * kappa_star * d)
    theta_star <- eta_star / rowSums(eta_star)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi, theta=theta_star, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * (1 + gamma) * phi) %*% theta_star), log=T)) + 
           dnorm(log(kappa_star), log_kappa_mean, log_kappa_sd, log=T)
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * (1 + gamma) * phi) %*% theta     ), log=T)) + 
           dnorm(log(kappa     ), log_kappa_mean, log_kappa_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa <- kappa_star
      eta <- eta_star
      theta <- theta_star
      cell_prob <- cell_prob_star
    }

    ### Sample alpha_det
    alpha_det_star <- rnorm(npcvs+1, alpha_det, alpha_det_tune)
    pdet_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        pdet_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_det_star)
      } # t
    } # i
    mh1 <- sum(dbinom(y, N, pdet_star, log=TRUE), na.rm=T) + 
           sum(dnorm(alpha_det_star, alpha_det_mean, alpha_det_sd, log=TRUE))
    mh2 <- sum(dbinom(y, N, pdet     , log=TRUE), na.rm=T) + 
           sum(dnorm(alpha_det     , alpha_det_mean, alpha_det_sd, log=TRUE))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha_det <- alpha_det_star
      pdet <- pdet_star
    }

    # Sample alpha_pcap
    alpha_cap_star <- rnorm(npcvs+1, alpha_cap, alpha_cap_tune)
    pcap_per_visit_star <- array(, dim=c(nsite, nyear-1, nreps))
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        pcap_per_visit_star[i,t,] <- inv.logit(cbind(1,w[i,t+1,,]) %*% alpha_cap_star)
      } # t
    } # i
    pcap_star <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        pcap_star[i,t] <- 1 - prod(1 - pcap_per_visit_star[i,t,])
      } # t
    } # i
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi=phi, theta=theta, pcap=pcap_star)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dnorm(alpha_cap_star, alpha_cap_mean, alpha_cap_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dnorm(alpha_cap     , alpha_cap_mean, alpha_cap_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      alpha_cap <- alpha_cap_star
      pcap_per_visit <- pcap_per_visit_star
      pcap <- pcap_star
      cell_prob <- cell_prob_star
    }

    # Save samples
    beta0_save[,k] <- beta0
    beta_gamma_save[,k] <- beta_gamma
    beta_phi_save[,k] <- beta_phi
    kappa_save[k] <- kappa
    alpha_det_save[,k] <- alpha_det
    alpha_cap_save[,k] <- alpha_cap
    N_mean_save[,k] <- colMeans(N)
    N_large_save[,k] <- N[site_large,]
    N_small_save[,k] <- N[site_small,]
    phi_mean_save[,k] <- colMeans(phi)
    phi_large_save[,k] <- phi[site_large,]
    phi_small_save[,k] <- phi[site_small,]
    gamma_mean_save[,k] <- colMeans(gamma)
    gamma_large_save[,k] <- gamma[site_large,]
    gamma_small_save[,k] <- gamma[site_small,]
  } # k

  list(beta0_save=beta0_save, 
       beta_gamma_save=beta_gamma_save, beta_phi_save=beta_phi_save, kappa_save=kappa_save, 
       alpha_det_save=alpha_det_save, alpha_cap_save=alpha_cap_save, 
       N_mean_save=N_mean_save, N_large_save=N_large_save, N_small_save=N_small_save, 
       phi_mean_save=phi_mean_save, phi_large_save=phi_large_save, phi_small_save=phi_small_save, 
       gamma_mean_save=gamma_mean_save, gamma_large_save=gamma_large_save, gamma_small_save=gamma_small_save)
} # sipm_mcmc

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
  sipm_mcmc(y=y, r=r, m=m, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='4.3.1 sipm_output.RData')

#==============
# Plot results
#==============
pdf(file='4.3.1.1 sipm_chains.pdf', width=10, height=10)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["reef"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta["density"]^""["["*gamma*"]"]), 
                     expression(beta["reef"]^""["["*gamma*"]"]), expression(beta["temperature"]^""["["*gamma*"]"]))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["density"]^""["["*phi*"]"]), 
                   expression(beta["reef"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]))
ylab_alpha_det <- c(expression(alpha[0]^"[det]"), expression(alpha[1]^"[det]"), expression(alpha[2]^"[det]"))
ylab_alpha_cap <- c(expression(alpha[0]^"[cap]"), expression(alpha[1]^"[cap]"), expression(alpha[2]^"[cap]"))

par(mfrow=c(6,3))
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
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
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
    tt[,j] <- out[[j]]$beta_gamma_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma[i] - yint * 2
  ymax <- beta_gamma[i] + yint * 2
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma[i], cex.main=2, line=1.2)
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
  ymin <- beta_phi[i] - yint * 8
  ymax <- beta_phi[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
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
ymin <- kappa - yint * 2
ymax <- kappa + yint * 2
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.2, las=2)
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
    tt[,j] <- out[[j]]$alpha_det_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha_det[i] - yint * c(4,2,2)[i]
  ymax <- alpha_det[i] + yint * c(4,2,2)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(2,1,1)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha_det[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha_det[i], cex.main=2, line=1.5)
  text(x=nmcmc*0.6, y=ymax, pos=1, cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_cap_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha_cap[i] - yint * 8
  ymax <- alpha_cap[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha_cap[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha_cap[i], cex.main=2, line=1.5)
  text(x=nmcmc*0.6, y=ymax, pos=1, cex=1.8, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=3, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=3, line=0.8, outer=T)

dev.off()

pdf(file='4.3.1.2 sipm_N phi gamma.pdf', width=10, height=8)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    N_mean <- out[[1]]$N_mean_save
    N_large <- out[[1]]$N_large_save
    N_small <- out[[1]]$N_small_save
    phi_mean <- out[[1]]$phi_mean_save
    phi_large <- out[[1]]$phi_large_save
    phi_small <- out[[1]]$phi_small_save
    gamma_mean <- out[[1]]$gamma_mean_save
    gamma_large <- out[[1]]$gamma_large_save
    gamma_small <- out[[1]]$gamma_small_save
  } else {
    N_mean <- cbind(N_mean, out[[1]]$N_mean_save)
    N_large <- cbind(N_large, out[[1]]$N_large_save)
    N_small <- cbind(N_small, out[[1]]$N_small_save)
    phi_mean <- cbind(phi_mean, out[[1]]$phi_mean_save)
    phi_large <- cbind(phi_large, out[[1]]$phi_large_save)
    phi_small <- cbind(phi_small, out[[1]]$phi_small_save)
    gamma_mean <- cbind(gamma_mean, out[[1]]$gamma_mean_save)
    gamma_large <- cbind(gamma_large, out[[1]]$gamma_large_save)
    gamma_small <- cbind(gamma_small, out[[1]]$gamma_small_save)
  }
} # i

site_large <- max(which(rowMeans(x[,,1]) == max(rowMeans(x[,,1]))))
site_small <- min(which(rowMeans(x[,,1]) == min(rowMeans(x[,,1]))))

par(mfrow=c(3,3))
par(mar=c(1.5,4,0,0))
par(oma=c(4,3,3,1))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,30), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(N), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,30,length.out=6), cex.axis=1.6, las=2)
axis(2, at=15, labels='Population Size', cex.axis=2.5, line=3.2, tick=F)
axis(3, at=5.5, labels='Mean', cex.axis=2.5, line=-.5, tick=F)

points(x=2, y=28, pch=16, cex=1, col='royalblue')
lines(x=c(1.5,2.5), y=c(28,28), lwd=1.2, col='royalblue')
vioplot(rnorm(10000,28,.5), at=5.5, add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
text(x=2.5, y=28, labels='True', pos=4, cex=1.8)
text(x=6, y=28, labels='Estimated', pos=4, cex=1.8)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,75), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_large), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(N[site_large,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,75,length.out=6), cex.axis=1.6, las=2)
axis(3, at=5.5, labels='Largest Site', cex.axis=2.5, line=-.5, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,15), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_small), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(N[site_small,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,15,length.out=6), cex.axis=1.6, las=2)
axis(3, at=5.5, labels='Smallest Site', cex.axis=2.5, line=-.5, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(phi), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Survival', cex.axis=2.5, line=3.2, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_large), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(phi[site_large,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), labels=rep('',6))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_small), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(phi[site_small,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), labels=rep('',6))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1.5), type='l', axes=F, xlab='', ylab='')
vioplot(t(gamma_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(gamma), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), cex.axis=1.6)
axis(2, at=seq(0,1.5,length.out=6), cex.axis=1.6, las=2)
axis(2, at=0.75, labels='Reproduction', cex.axis=2.5, line=3.2, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1.5), type='l', axes=F, xlab='', ylab='')
vioplot(t(gamma_large), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(gamma[site_large,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), cex.axis=1.6)
axis(2, at=seq(0,1.5,length.out=6), labels=rep('',6))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1.5), type='l', axes=F, xlab='', ylab='')
vioplot(t(gamma_small), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(gamma[site_small,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), cex.axis=1.6)
axis(2, at=seq(0,1.5,length.out=6), labels=rep('',6))

title(xlab='Year', cex.lab=3, line=2.4, outer=T)

dev.off()

pdf(file='4.3.1.3 sipm_correlation.pdf', width=8, height=8)

for (i in 1:chain) {
  if (i == 1) {
    beta_phi0_post <- out[[i]]$beta_phi[1,(nmcmc*0.2+1):nmcmc]
    beta_gamma0_post <- out[[i]]$beta_gamma[1,(nmcmc*0.2+1):nmcmc]
  } else {
    beta_phi0_post <- c(beta_phi0_post, out[[i]]$beta_phi[1,(nmcmc*0.2+1):nmcmc])
    beta_gamma0_post <- c(beta_gamma0_post, out[[i]]$beta_gamma[1,(nmcmc*0.2+1):nmcmc])
  }
} # i

nplot <- min(c(10000, nmcmc*0.2*chain))
sel <- sort(sample(1:(nmcmc*0.2*chain), nplot, replace=F))
a <- beta_phi0_post[sel]
b <- beta_gamma0_post[sel]
a2 <- (a - min(a)) / (max(a) - min(a))
b2 <- (b - min(b)) / (max(b) - min(b))
ab <- matrix(, nplot, nplot)
for (i in 1:nplot) {
  for (j in 1:nplot) {
    ab[i,j] <- sqrt((a2[i] - a2[j]) ^ 2 + (b2[i] - b2[j]) ^ 2)
  } # j
} # i

ab01 <- ifelse(ab < 0.08, 1, 0)
ab_scale <- rowSums(ab01) - min(rowSums(ab01)) + 1
col2 <- c(colorRampPalette(c('navy','seagreen'))(floor(max(ab_scale)*0.6)), 
          colorRampPalette(c('seagreen','yellow'))(ceiling(max(ab_scale)*0.4)))

par(mfrow=c(1,1))
par(mar=c(5,6.5,1,1))

plot(a, b, col=col2[ab_scale], pch=16, cex=1.4, xlim=c(0.2,0.8), ylim=c(-0.9,-0.7), axes=F, xlab='', ylab='')
axis(1, at=seq( 0.2, 0.8,length.out=5), cex.axis=1.5)
axis(2, at=seq(-0.9,-0.7,length.out=5), cex.axis=1.5, las=2)
axis(1, at= 0.5, labels=ylab_beta_phi  [1], cex.axis=2, line=2.8, tick=F)
axis(2, at=-0.8, labels=ylab_beta_gamma[1], cex.axis=2, line=3.6, tick=F, las=2)
box()

text(x=0.5, y=-0.7, labels=paste('correlation = ', round(cor(a, b), digits=2), sep=''), pos=1, cex=1.8)

dev.off()


