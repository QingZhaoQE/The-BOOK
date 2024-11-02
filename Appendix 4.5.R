#=================================================================================================
# Spatial integrated population model (SIPM) with age-specific survival and movement
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 4_Integrated population models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions
library(LaplacesDemon) # for generating random numbers from categorical distributions (rcat)

set.seed(7)

# Basic values
nsite <- 60 # number of sites
nyear <- 10 # number of years
nreps <- 5  # number of within-season replicates
ncovs <- 2  # number of environmental covariates
npcvs <- 1  # number of observational covariates

beta0 <- c(3.4, 1, -0.8)                # intercept and slopes for initial abundance (i.e., number of adults)
beta_gamma <-  c(-0.7, -0.3, 0.4, -0.2) # intercept and slopes for reproduction
beta_phi_adu <- c(1  , -0.5, 0.6, -0.3) # intercept and slopes for adult survival
beta_phi_juv <- c(0.4, -0.7, 0.8, -0.4) # intercept and slopes for juvenile survival
kappa_adu <- 2                          # distance effect on adult movement
kappa_juv <- 0.5                        # distance effect on juvenile movement
alpha_det <- c( 0.9, -0.4)              # intercept and slopes for detection probability
alpha_cap <- c(-3  , -0.8)              # intercept and slopes for capture probability

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
eta_adu <- exp(-1 * kappa_adu * d)      # unstandardized adult colonization rate
theta_adu <- eta_adu / rowSums(eta_adu) #   standardized adult colonization rate
eta_juv <- exp(-1 * kappa_juv * d)      # unstandardized juvenile colonization rate
theta_juv <- eta_juv / rowSums(eta_juv) #   standardized juvenile colonization rate

nindv <- 15000 # maximum number of individuals that ever exist
z <- matrix(0, nindv, nyear) # individual survival-location history
for (i in 1:nsite) {
  if (i == 1) { tt <- rep(i+nsite, N[i,1]) } else { tt <- c(tt, rep(i+nsite, N[i,1])) }
}
z[1:length(tt),1] <- tt

gamma <- matrix(, nsite, nyear-1) # reproduction rate
phi_adu <- matrix(, nsite, nyear-1) # adult survival probability
phi_juv <- matrix(, nsite, nyear-1) # juvenile survival probability
R <- matrix(, nsite, nyear-1) # number of offspring

for (t in 2:nyear) {
  gamma[,t-1] <- exp(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_gamma)
  R[,t-1] <- rpois(nsite, N[,t-1] * gamma[,t-1])

  for (i in 1:nsite) {
    if (i == 1) { tt <- rep(i, R[i,t-1]) } else { tt <- c(tt, rep(i, R[i,t-1])) }
  } # i
  z[max(which(rowMeans(z,na.rm=T)>0))+c(1:length(tt)), t-1] <- tt

  phi_adu[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_adu)
  phi_juv[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_juv)
  transition_matrix <- matrix(, nsite*2, nsite*2+1) # matrix of transition (survival & movement) probability
  for (i in 1:nsite) {
    transition_matrix[i,1:nsite] <- 0
    transition_matrix[i,(nsite+1):(nsite*2)] <- phi_juv[i,t-1] * theta_juv[i,]
    transition_matrix[i,nsite*2+1] <- 1 - phi_juv[i,t-1]
    transition_matrix[nsite+i,1:nsite] <- 0
    transition_matrix[nsite+i,(nsite+1):(nsite*2)] <- phi_adu[i,t-1] * theta_adu[i,]
    transition_matrix[nsite+i,nsite*2+1] <- 1 - phi_adu[i,t-1]
  } # i

  for (v in 1:nindv) {
    if (z[v,t-1] %in% 1:(nsite*2)) {
      z[v,t] <- rcat(1, p=transition_matrix[z[v,t-1],])
    }
  } # v
  z[which(!(z %in% 1:(nsite*2)))] <- 0

  for (i in 1:nsite) {
    N[i,t] <- length(which(z[,t] %in% c(i,nsite+i)))
  }
} # t
z[which(!(z %in% 1:(nsite*2)))] <- 0
z <- z[which(rowMeans(z, na.rm=T) != 0),]

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
pcap_per_visit <- array(, dim=c(nsite*2, nyear, nreps)) # per visit capture probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    pcap_per_visit[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_cap)
  } # t
} # i
pcap_per_visit[(nsite+1):(nsite*2),,] <- pcap_per_visit[1:nsite,,]

c_all <- matrix(, dim(z)[1], nyear) # capture history of all individuals
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
if_cap <- ifelse(apply(c_all, 1, sum) > 0, 1, 0)

c <- c_all[which(if_cap==1),] # capture history of individuals that have ever been captured
f <- numeric(dim(c)[1]) # year of first capture
for (v in 1:length(f)) {
  f[v] <- min(which(c[v,] %in% c(1:(nsite*2))))
} # v
c <- c[which(f %in% c(1:(nyear-1))),]

# Extend individual capture history
for (v in 1:dim(c)[1]) {
  t1 <- c[v,]
  ycap <- which(t1 %in% c(1:(nsite*2)))
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
  f_ext[v] <- min(which(c_ext[v,] %in% c(1:(nsite*2))))
} # v

# Converting individual capture history to m-array format
for (t in 2:nyear) {
  t1 <- paste('juv_site', '_', 1:nsite, ', year', t, sep='')
  t2 <- paste('adu_site', '_', 1:nsite, ', year', t, sep='')
  tt <- c(t1, t2)
  if (t == 2) {
    nm1 <- tt
  } else { 
    nm1 <- c(nm1, tt)
  }
} # t
nm1 <- c(nm1, 'not seen')
nm2 <- c(paste('juv_site', 1:nsite, sep='_'), paste('adu_site', 1:nsite, sep='_'))

m <- array(0, # m-array of capture-recapture
  dim=c(nyear-1, nsite*2*(nyear-1)+1, nsite*2), 
  dimnames=list(paste('year', 1:(nyear-1), sep='_'), nm1, nm2))
for (t1 in 1:(nyear-1)) {
  temp1 <- c_ext[which(f_ext == t1),]
  for (i1 in 1:(nsite*2)) {
    temp2 <- temp1[which(temp1[,t1] == i1),]
    nband <- length(which(temp1[,t1] == i1))
    if (nband == 1) {
      temp2 <- matrix(temp2, nrow=1, ncol=length(temp2))
    }
    for (t2 in (t1+1):nyear) {
      for (i2 in 1:(nsite*2)) {
        m[t1,(t2-2)*(nsite*2)+i2,i1] <- length(which(temp2[,t2] == i2))
      } # i2
    } # t2
    m[t1,(nsite*2)*(nyear-1)+1,i1] <- nband - sum(m[t1,1:((nsite*2)*(nyear-1)),i1])
  } # i1
} # t1

#=======================
# Define MCMC algorithm
#=======================
# Function to calculate the cell probability of m-array
cell_prob_fun <- function(nsite, nyear, phi_adu, phi_juv, theta_adu, theta_juv, pcap) { 
  qcap <- 1 - pcap
  phi_first_year <- rbind(phi_juv, phi_adu)
  phi_from_second_year <- rbind(phi_adu, phi_adu)
  theta_first_year <- rbind(cbind(matrix(0,nsite,nsite),theta_juv), cbind(matrix(0,nsite,nsite),theta_adu))
  theta_from_second_year <- rbind(cbind(matrix(0,nsite,nsite),theta_adu), cbind(matrix(0,nsite,nsite),theta_adu))
  cell_prob_p <- cell_prob_q <- array(0, dim=c(nyear-1, (nsite*2)*(nyear-1)+1, (nsite*2)))
  for (i in 1:(nsite*2)) {
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,(nsite*2)*(t-1)+1:(nsite*2),i] <- phi_first_year[i,t] * theta_first_year[i,] * pcap[,t]
      cell_prob_q[t,(nsite*2)*(t-1)+1:(nsite*2),i] <- phi_first_year[i,t] * theta_first_year[i,] * qcap[,t]
    } # t
    for (t1 in 1:(nyear-2)) {
      for (t2 in (t1+1):(nyear-1)) {
        cell_prob_p[t1,(nsite*2)*(t2-1)+1:(nsite*2),i] <- 
          ((cell_prob_q[t1,(nsite*2)*(t2-2)+1:(nsite*2),i] * phi_from_second_year[,t2]) %*% theta_from_second_year) * pcap[,t2]
        cell_prob_q[t1,(nsite*2)*(t2-1)+1:(nsite*2),i] <- 
          ((cell_prob_q[t1,(nsite*2)*(t2-2)+1:(nsite*2),i] * phi_from_second_year[,t2]) %*% theta_from_second_year) * qcap[,t2]
      } # t2
    } # t1
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,(nsite*2)*(nyear-1)+1,i] <- 1 - sum(cell_prob_p[t,1:((nsite*2)*(nyear-1)),i])
    } # t
  } # i
  cell_prob <- cell_prob_p

  return(cell_prob)
} # cell_prob_fun

age_sipm_mcmc <- function(y, r, m, x, w, d, nmcmc) {

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
  beta_phi_adu_save <- matrix(0, ncovs+2, nmcmc)
  beta_phi_juv_save <- matrix(0, ncovs+2, nmcmc)
  kappa_adu_save <- rep(0, nmcmc)
  kappa_juv_save <- rep(0, nmcmc)
  alpha_det_save <- matrix(0, npcvs+1, nmcmc)
  alpha_cap_save <- matrix(0, npcvs+1, nmcmc)
  N_mean_save <- matrix(0, nyear, nmcmc)
  N_large_save <- matrix(0, nyear, nmcmc)
  N_small_save <- matrix(0, nyear, nmcmc)
  phi_adu_mean_save <- matrix(0, nyear-1, nmcmc)
  phi_adu_large_save <- matrix(0, nyear-1, nmcmc)
  phi_adu_small_save <- matrix(0, nyear-1, nmcmc)
  phi_juv_mean_save <- matrix(0, nyear-1, nmcmc)
  phi_juv_large_save <- matrix(0, nyear-1, nmcmc)
  phi_juv_small_save <- matrix(0, nyear-1, nmcmc)
  gamma_mean_save <- matrix(0, nyear-1, nmcmc)
  gamma_large_save <- matrix(0, nyear-1, nmcmc)
  gamma_small_save <- matrix(0, nyear-1, nmcmc)
  N_save <- array(0, dim=c(nsite, nyear, 200)) # for posterior predictive check of y
  pdet_save <- array(0, dim=c(nsite, nyear, nreps, 200)) # for posterior predictive check of y
  gamma_save <- array(0, dim=c(nsite, nyear-1, 200)) # for posterior predictive check of r
  cell_prob_save <- array(0, dim=c(nyear-1, nsite*2*(nyear-1)+1, nsite*2, 200)) # for posterior predictive check for M-array

  # Prior
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_gamma_mean <- rep(0, ncovs+2)
  beta_gamma_sd <- 2
  beta_phi_adu_mean <- rep(0, ncovs+2)
  beta_phi_adu_sd <- 2
  beta_phi_juv_mean <- rep(0, ncovs+2)
  beta_phi_juv_sd <- 2
  log_kappa_adu_mean <- 0
  log_kappa_adu_sd <- 2
  log_kappa_juv_mean <- 0
  log_kappa_juv_sd <- 2
  alpha_det_mean <- rep(0, npcvs+1)
  alpha_det_sd <- 2
  alpha_cap_mean <- rep(0, npcvs+1)
  alpha_cap_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_gamma <- rep(0, ncovs+2)
  beta_phi_adu <- rep(0, ncovs+2)
  beta_phi_juv <- rep(0, ncovs+2)
  kappa_adu <- 1
  kappa_juv <- 1
  alpha_det <- rep(0, npcvs+1)
  alpha_cap <- rep(0, npcvs+1)
  lambda0 <- exp(cbind(1,x[,1,]) %*% beta0)
  gamma <- matrix(0.5, nsite, nyear-1)
  phi_adu <- matrix(0.5, nsite, nyear-1)
  phi_juv <- matrix(0.5, nsite, nyear-1)
  eta_adu <- exp(-1 * kappa_adu * d)
  theta_adu <- eta_adu / rowSums(eta_adu)
  eta_juv <- exp(-1 * kappa_juv * d)
  theta_juv <- eta_juv / rowSums(eta_juv)
  N <- round((ymax + 1) / 0.5)
  pdet <- array(0.5, dim=c(nsite, nyear, nreps))
  pcap_per_visit <- array(0.5, dim=c(nsite*2, nyear-1, nreps))
  pcap <- matrix(, nsite*2, nyear-1)
  for (i in 1:(nsite*2)) {
    for (t in 1:(nyear-1)) {
      pcap[i,t] <- 1 - prod(1 - pcap_per_visit[i,t,])
    } # t
  } # i
  cell_prob <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_adu=phi_adu, phi_juv=phi_juv, 
                             theta_adu=theta_adu, theta_juv=theta_juv, pcap=pcap)

  # Tuning factor
  beta0_tune <- 0.08
  beta_gamma_tune <- 0.02
  beta_phi_adu_tune <- c(0.10, 0.10, 0.07, 0.07)
  beta_phi_juv_tune <- c(0.10, 0.10, 0.07, 0.07)
  kappa_adu_tune <- 0.035
  kappa_juv_tune <- 0.035
  alpha_det_tune <- c(0.050, 0.035)
  alpha_cap_tune <- 0.050
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- matrix(rpois(nsite*nyear, N + N_tune), nsite, nyear)
    mh1 <- apply(dbinom(y, N_star, pdet, log=T), 1:2, sum) + 
           dpois(N_star, cbind(lambda0, t(t(N[,-nyear] * gamma * phi_juv) %*% theta_juv + t(N[,-nyear] * phi_adu) %*% theta_adu)), log=T) + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(dbinom(y, N     , pdet, log=T), 1:2, sum) + 
           dpois(N     , cbind(lambda0, t(t(N[,-nyear] * gamma * phi_juv) %*% theta_juv + t(N[,-nyear] * phi_adu) %*% theta_adu)), log=T) + 
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
           sum(dpois(N[,-1], t(t(N[,-nyear] * gamma_star * phi_juv) %*% theta_juv + t(N[,-nyear] * phi_adu) %*% theta_adu), log=T)) + 
           sum(dnorm(beta_gamma_star, beta_gamma_mean, beta_gamma_sd, log=T))
    mh2 <- sum(dbinom(r, r + y[,-nyear,], gamma      / (gamma      + 1), log=T)) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * gamma      * phi_juv) %*% theta_juv + t(N[,-nyear] * phi_adu) %*% theta_adu), log=T)) + 
           sum(dnorm(beta_gamma     , beta_gamma_mean, beta_gamma_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_gamma <- beta_gamma_star
    }

    # Sample beta_phi_adu & beta_phi_juv
    beta_phi_adu_star <- rnorm(ncovs+2, beta_phi_adu, beta_phi_adu_tune)
    beta_phi_juv_star <- rnorm(ncovs+2, beta_phi_juv, beta_phi_juv_tune)
    phi_adu <- phi_adu_star <- phi_juv <- phi_juv_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi_adu     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_adu     )
      phi_adu_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_adu_star)
      phi_juv     [,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_juv     )
      phi_juv_star[,t-1] <- inv.logit(cbind(1,(N[,t-1]-lambda0)/lambda0,x[,t,]) %*% beta_phi_juv_star)
    } # t
    cell_prob      <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_adu=phi_adu     , phi_juv=phi_juv     , theta_adu=theta_adu, theta_juv=theta_juv, pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_adu=phi_adu_star, phi_juv=phi_juv_star, theta_adu=theta_adu, theta_juv=theta_juv, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * gamma * phi_juv_star) %*% theta_juv + t(N[,-nyear] * phi_adu_star) %*% theta_adu), log=T)) + 
           sum(dnorm(beta_phi_adu_star, beta_phi_adu_mean, beta_phi_adu_sd, log=T)) + 
           sum(dnorm(beta_phi_juv_star, beta_phi_juv_mean, beta_phi_juv_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * gamma * phi_juv     ) %*% theta_juv + t(N[,-nyear] * phi_adu     ) %*% theta_adu), log=T)) + 
           sum(dnorm(beta_phi_adu     , beta_phi_adu_mean, beta_phi_adu_sd, log=T)) + 
           sum(dnorm(beta_phi_juv     , beta_phi_juv_mean, beta_phi_juv_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_phi_adu <- beta_phi_adu_star
      beta_phi_juv <- beta_phi_juv_star
      phi_adu <- phi_adu_star
      phi_juv <- phi_juv_star
      cell_prob <- cell_prob_star
    }

    # Sample kappa_adu & kappa_juv
    kappa_adu_star <- exp(rnorm(1, log(kappa_adu), kappa_adu_tune))
    eta_adu_star <- exp(-1 * kappa_adu_star * d)
    theta_adu_star <- eta_adu_star / rowSums(eta_adu_star)
    kappa_juv_star <- exp(rnorm(1, log(kappa_juv), kappa_juv_tune))
    eta_juv_star <- exp(-1 * kappa_juv_star * d)
    theta_juv_star <- eta_juv_star / rowSums(eta_juv_star)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_adu=phi_adu, phi_juv=phi_juv, theta_adu=theta_adu_star, theta_juv=theta_juv_star, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * gamma * phi_juv) %*% theta_juv_star + t(N[,-nyear] * phi_adu) %*% theta_adu_star), log=T)) + 
           dnorm(log(kappa_adu_star), log_kappa_adu_mean, log_kappa_adu_sd, log=T) + 
           dnorm(log(kappa_juv_star), log_kappa_juv_mean, log_kappa_juv_sd, log=T)
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1], t(t(N[,-nyear] * gamma * phi_juv) %*% theta_juv      + t(N[,-nyear] * phi_adu) %*% theta_adu     ), log=T)) + 
           dnorm(log(kappa_adu     ), log_kappa_adu_mean, log_kappa_adu_sd, log=T) + 
           dnorm(log(kappa_juv     ), log_kappa_juv_mean, log_kappa_juv_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa_adu <- kappa_adu_star
      eta_adu <- eta_adu_star
      theta_adu <- theta_adu_star
      kappa_juv <- kappa_juv_star
      eta_juv <- eta_juv_star
      theta_juv <- theta_juv_star
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
    pcap_per_visit_star <- array(, dim=c(nsite*2, nyear-1, nreps))
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        pcap_per_visit_star[i,t,] <- inv.logit(cbind(1,w[i,t+1,,]) %*% alpha_cap_star)
      } # t
    } # i
    pcap_per_visit_star[(nsite+1):(nsite*2),,] <- pcap_per_visit_star[1:nsite,,]
    pcap_star <- matrix(, nsite*2, nyear-1)
    for (i in 1:(nsite*2)) {
      for (t in 1:(nyear-1)) {
        pcap_star[i,t] <- 1 - prod(1 - pcap_per_visit_star[i,t,])
      } # t
    } # i
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_adu=phi_adu, phi_juv=phi_juv, theta_adu=theta_adu, theta_juv=theta_juv, pcap=pcap_star)
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
    beta_phi_adu_save[,k] <- beta_phi_adu
    beta_phi_juv_save[,k] <- beta_phi_juv
    kappa_adu_save[k] <- kappa_adu
    kappa_juv_save[k] <- kappa_juv
    alpha_det_save[,k] <- alpha_det
    alpha_cap_save[,k] <- alpha_cap
    N_mean_save[,k] <- colMeans(N)
    N_large_save[,k] <- N[site_large,]
    N_small_save[,k] <- N[site_small,]
    phi_adu_mean_save[,k] <- colMeans(phi_adu)
    phi_adu_large_save[,k] <- phi_adu[site_large,]
    phi_adu_small_save[,k] <- phi_adu[site_small,]
    phi_juv_mean_save[,k] <- colMeans(phi_juv)
    phi_juv_large_save[,k] <- phi_juv[site_large,]
    phi_juv_small_save[,k] <- phi_juv[site_small,]
    gamma_mean_save[,k] <- colMeans(gamma)
    gamma_large_save[,k] <- gamma[site_large,]
    gamma_small_save[,k] <- gamma[site_small,]
    if (nmcmc > 200) {
      if (k > nmcmc-200) {
        N_save[,,k-nmcmc+200] <- N
        pdet_save[,,,k-nmcmc+200] <- pdet
        gamma_save[,,k-nmcmc+200] <- gamma
        cell_prob_save[,,,k-nmcmc+200] <- cell_prob
      }
    }
  } # k

  list(beta0_save=beta0_save, 
       beta_gamma_save=beta_gamma_save, 
       beta_phi_adu_save=beta_phi_adu_save, beta_phi_juv_save=beta_phi_juv_save, 
       kappa_adu_save=kappa_adu_save, kappa_juv_save=kappa_juv_save, 
       alpha_det_save=alpha_det_save, alpha_cap_save=alpha_cap_save, 
       N_mean_save=N_mean_save, N_large_save=N_large_save, N_small_save=N_small_save, 
       phi_adu_mean_save=phi_adu_mean_save, phi_adu_large_save=phi_adu_large_save, phi_adu_small_save=phi_adu_small_save, 
       phi_juv_mean_save=phi_juv_mean_save, phi_juv_large_save=phi_juv_large_save, phi_juv_small_save=phi_juv_small_save, 
       gamma_mean_save=gamma_mean_save, gamma_large_save=gamma_large_save, gamma_small_save=gamma_small_save, 
       N_save=N_save, pdet_save=pdet_save, gamma_save=gamma_save, cell_prob_save=cell_prob_save)
} # age_sipm_mcmc

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
  age_sipm_mcmc(y=y, r=r, m=m, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='4.3.2 age sipm_output.RData')

#==============
# Plot results
#==============
pdf(file='4.3.2.1 age sipm_chains.pdf', width=10, height=10)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["reef"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta["density"]^""["["*gamma*"]"]), 
                     expression(beta["reef"]^""["["*gamma*"]"]), expression(beta["temperature"]^""["["*gamma*"]"]))
ylab_beta_phi_adu <- c(expression(beta["0,adu"]^""["["*phi*"]"]), expression(beta["density,adu"]^""["["*phi*"]"]), 
                       expression(beta["reef,adu"]^""["["*phi*"]"]), expression(beta["temperature,adu"]^""["["*phi*"]"]))
ylab_beta_phi_juv <- c(expression(beta["0,juv"]^""["["*phi*"]"]), expression(beta["density,juv"]^""["["*phi*"]"]), 
                       expression(beta["reef,juv"]^""["["*phi*"]"]), expression(beta["temperature,juv"]^""["["*phi*"]"]))
ylab_alpha_det <- c(expression(alpha[0]^"[det]"), expression(alpha[1]^"[det]"))
ylab_alpha_cap <- c(expression(alpha[0]^"[cap]"), expression(alpha[1]^"[cap]"))

par(mfrow=c(7,3))
par(mar=c(1,3,3,1))
par(oma=c(4,3,0,0))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta0_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0[i] - yint * c(8,4,8)[i]
  ymax <- beta0[i] + yint * c(8,4,8)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(4,2,4)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta0[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i], cex.main=2, line=1.4)
  text(x=nmcmc*0.4, y=beta0[i]+yint*c(6,3,6)[i], pos=4,  cex=1.6, 
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
  text(x=nmcmc*0.4, y=beta_gamma[i]+yint*1.5, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_adu_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi_adu[i] - yint * 8
  ymax <- beta_phi_adu[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi_adu[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi_adu[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_phi_adu[i]+yint*6, pos=4,  cex=1.6, 
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
  text(x=nmcmc*0.4, y=beta_phi_juv[i]+yint*6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$kappa_adu_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa_adu - yint * 8
ymax <- kappa_adu + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa_adu, col='grey16', lwd=1.5)
title(main=expression(kappa[adu]), cex.main=2, line=1.2)
text(x=nmcmc*0.4, y=kappa_adu+yint*6, pos=4,  cex=1.6, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$kappa_juv_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa_juv - yint * 2
ymax <- kappa_juv + yint * 2
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa_juv, col='grey16', lwd=1.5)
title(main=expression(kappa[juv]), cex.main=2, line=1.2)
text(x=nmcmc*0.4, y=kappa_juv+yint*1.5, pos=4,  cex=1.6, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_det_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha_det[i] - yint * c(4,2)[i]
  ymax <- alpha_det[i] + yint * c(4,2)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  if (i == npcvs+1) {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  }
  axis(2, at=round(seq(ymin, ymax, yint*c(2,1)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha_det[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha_det[i], cex.main=2, line=1.5)
  text(x=nmcmc*0.4, y=alpha_det[i]+yint*c(3,1.5)[i], pos=4,  cex=1.6, 
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
  text(x=nmcmc*0.4, y=alpha_cap[i]+yint*6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=3, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=3, line=0.8, outer=T)

dev.off()

pdf(file='4.3.2.2 age sipm_N phi_adu phi_juv gamma.pdf', width=10, height=10)

library(vioplot) # for making violin plots

for (i in 1:chain) {
  if (i == 1) {
    N_mean <- out[[1]]$N_mean_save
    N_large <- out[[1]]$N_large_save
    N_small <- out[[1]]$N_small_save
    phi_adu_mean <- out[[1]]$phi_adu_mean_save
    phi_adu_large <- out[[1]]$phi_adu_large_save
    phi_adu_small <- out[[1]]$phi_adu_small_save
    phi_juv_mean <- out[[1]]$phi_juv_mean_save
    phi_juv_large <- out[[1]]$phi_juv_large_save
    phi_juv_small <- out[[1]]$phi_juv_small_save
    gamma_mean <- out[[1]]$gamma_mean_save
    gamma_large <- out[[1]]$gamma_large_save
    gamma_small <- out[[1]]$gamma_small_save
  } else {
    N_mean <- cbind(N_mean, out[[1]]$N_mean_save)
    N_large <- cbind(N_large, out[[1]]$N_large_save)
    N_small <- cbind(N_small, out[[1]]$N_small_save)
    phi_adu_mean <- cbind(phi_adu_mean, out[[1]]$phi_adu_mean_save)
    phi_adu_large <- cbind(phi_adu_large, out[[1]]$phi_adu_large_save)
    phi_adu_small <- cbind(phi_adu_small, out[[1]]$phi_adu_small_save)
    phi_juv_mean <- cbind(phi_juv_mean, out[[1]]$phi_juv_mean_save)
    phi_juv_large <- cbind(phi_juv_large, out[[1]]$phi_juv_large_save)
    phi_juv_small <- cbind(phi_juv_small, out[[1]]$phi_juv_small_save)
    gamma_mean <- cbind(gamma_mean, out[[1]]$gamma_mean_save)
    gamma_large <- cbind(gamma_large, out[[1]]$gamma_large_save)
    gamma_small <- cbind(gamma_small, out[[1]]$gamma_small_save)
  }
} # i

site_large <- max(which(rowMeans(x[,,1]) == max(rowMeans(x[,,1]))))
site_small <- min(which(rowMeans(x[,,1]) == min(rowMeans(x[,,1]))))

par(mfrow=c(4,3))
par(mar=c(1.5,4,0,0))
par(oma=c(4,3,3,1))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,60), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(N), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,60,length.out=6), cex.axis=1.6, las=2)
axis(2, at=30, labels='Population Size', cex.axis=2.5, line=3.2, tick=F)
axis(3, at=5.5, labels='Mean', cex.axis=2.5, line=-.5, tick=F)

points(x=2, y=56, pch=16, cex=1, col='royalblue')
lines(x=c(1.5,2.5), y=c(56,56), lwd=1.2, col='royalblue')
vioplot(rnorm(10000,56,.5), at=5.5, add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
text(x=2.5, y=56, labels='True', pos=4, cex=1.8)
text(x=6, y=56, labels='Estimated', pos=4, cex=1.8)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,100), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_large), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(N[site_large,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,100,length.out=6), cex.axis=1.6, las=2)
axis(3, at=5.5, labels='Largest Site', cex.axis=2.5, line=-.5, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,30), type='l', axes=F, xlab='', ylab='')
vioplot(t(N_small), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(N[site_small,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,30,length.out=6), cex.axis=1.6, las=2)
axis(3, at=5.5, labels='Smallest Site', cex.axis=2.5, line=-.5, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_adu_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(phi_adu), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Adult Survival', cex.axis=2.5, line=3.2, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_adu_large), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(phi_adu[site_large,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), labels=rep('',6))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_adu_small), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(phi_adu[site_small,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), labels=rep('',6))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_juv_mean), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(colMeans(phi_juv), type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), cex.axis=1.6, las=2)
axis(2, at=0.5, labels='Juvenile Survival', cex.axis=2.5, line=3.2, tick=F)

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_juv_large), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(phi_juv[site_large,], type='o', pch=16, col='royalblue')
axis(1, at=seq(0,nyear,2), labels=rep('',6))
axis(2, at=seq(0,1,length.out=6), labels=rep('',6))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0,1), type='l', axes=F, xlab='', ylab='')
vioplot(t(phi_juv_small), add=T, col='lightcoral', rectCol=NA, lineCol=NA, border=NA)
lines(phi_juv[site_small,], type='o', pch=16, col='royalblue')
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

#============================
# Posterior predictive check
#============================
# Count data
N_post <- array(, dim=c(nsite, nyear, 600))
N_post[,,  1:200] <- out[[1]]$N_save
N_post[,,201:400] <- out[[2]]$N_save
N_post[,,401:600] <- out[[3]]$N_save
pdet_post <- array(, dim=c(nsite, nyear, nreps, 600))
pdet_post[,,,  1:200] <- out[[1]]$pdet_save
pdet_post[,,,201:400] <- out[[2]]$pdet_save
pdet_post[,,,401:600] <- out[[3]]$pdet_save
y_exp <- y_sim <- y_dobs <- y_dsim <- array(, dim=c(nsite, nyear, nreps, 600))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (j in 1:nreps) {
      y_exp[i,t,j,] <- N_post[i,t,] * pdet_post[i,t,j,]
      y_sim[i,t,j,] <- rbinom(600, N_post[i,t,], pdet_post[i,t,j,])
      y_dobs[i,t,j,] <- (y    [i,t,j ] - y_exp[i,t,j,]) ^ 2 / var(y_exp[i,t,j,])
      y_dsim[i,t,j,] <- (y_sim[i,t,j,] - y_exp[i,t,j,]) ^ 2 / var(y_exp[i,t,j,])
    } # j
  } # t
} # i
y_ppp <- length(which(y_dsim > y_dobs)) / length(which(y_dsim != y_dobs))

# Age-ratio data
gamma_post <- array(, dim=c(nsite, nyear-1, 600))
gamma_post[,,  1:200] <- out[[1]]$gamma_save
gamma_post[,,201:400] <- out[[2]]$gamma_save
gamma_post[,,401:600] <- out[[3]]$gamma_save
r_exp <- r_sim <- r_dobs <- r_dsim <- array(, dim=c(nsite, nyear-1, nreps, 600))
for (i in 1:nsite) {
  for (t in 1:(nyear-1)) {
    for (j in 1:nreps) {
      r_exp[i,t,j,] <- (r[i,t,j] + y[i,t,j]) * gamma_post[i,t,] / (gamma_post[i,t,] + 1)
      r_sim[i,t,j,] <- rbinom(600, r[i,t,j] + y[i,t,j], gamma_post[i,t,] / (gamma_post[i,t,] + 1))
      r_dobs[i,t,j,] <- (r    [i,t,j ] - r_exp[i,t,j,]) ^ 2 / var(r_exp[i,t,j,])
      r_dsim[i,t,j,] <- (r_sim[i,t,j,] - r_exp[i,t,j,]) ^ 2 / var(r_exp[i,t,j,])
    } # j
  } # t
} # i
r_ppp <- length(which(r_dsim > r_dobs)) / length(which(r_dsim != r_dobs))

# Capture-recapture data
cell_prob_post <- array(, dim=c(nyear-1, nsite*2*(nyear-1)+1, nsite*2, 600))
cell_prob_post[,,,  1:200] <- out[[1]]$cell_prob_save
cell_prob_post[,,,201:400] <- out[[2]]$cell_prob_save
cell_prob_post[,,,401:600] <- out[[3]]$cell_prob_save
nband <- apply(m, c(1,3), sum)
m_exp <- m_sim <- array(, dim=c(nyear-1, nsite*2*(nyear-1)+1, nsite*2, 600))
for (i in 1:(nsite*2)) {
  for (t in 1:(nyear-1)) {
    for (k in 1:600) {
      m_exp[t,,i,k] <- nband[t,i] * cell_prob_post[t,,i,k]
      m_sim[t,,i,k] <- rmultinom(1, nband[t,i], cell_prob_post[t,,i,k])
    } # k
  } # t
} # i
m_juv_dobs <- m_juv_dsim <- m_adu_dobs <- m_adu_dsim <- array(, dim=c(nsite, nyear-1, 600))
for (i in 1:nsite) {
  for (t in 1:(nyear-1)) {
    for (k in 1:600) {
      m_juv_dobs[i,t,k] <- sum((sqrt(m    [t,,i  ]) - sqrt(m_exp[t,,i,k])) ^ 2)
      m_juv_dsim[i,t,k] <- sum((sqrt(m_sim[t,,i,k]) - sqrt(m_exp[t,,i,k])) ^ 2)
      m_adu_dobs[i,t,k] <- sum((sqrt(m    [t,,nsite+i  ]) - sqrt(m_exp[t,,nsite+i,k])) ^ 2)
      m_adu_dsim[i,t,k] <- sum((sqrt(m_sim[t,,nsite+i,k]) - sqrt(m_exp[t,,nsite+i,k])) ^ 2)
    } # k
  } # t
} # i
m_juv_ppp <- length(which(m_juv_dsim > m_juv_dobs)) / length(which(m_juv_dsim != m_juv_dobs))
m_adu_ppp <- length(which(m_adu_dsim > m_adu_dobs)) / length(which(m_adu_dsim != m_adu_dobs))

pdf(file='4.3.2.3 age sipm_ppp.pdf', width=8, height=8)

col_point <- 'lightcoral'
col_line  <- 'royalblue'

par(mfrow=c(2,2))
par(mar=c(2,4,3,1))
par(oma=c(3,3,0,0))

y_sample <- array(c(rep(1,10000),rep(0,nsite*nyear*nreps*600-10000)), dim=c(nsite, nyear, nreps, 600))
plot(y_dsim[which(y_sample == 1)] ~ y_dobs[which(y_sample == 1)], 
     xlim=c(0,600), ylim=c(0,600), axes=F, xlab='', ylab='', pch=16, cex=1.2, col=col_point)
abline(0, 1, col=col_line, lwd=1.8)
axis(1, at=seq(0,600,length.out=4), cex.axis=1.4)
axis(2, at=seq(0,600,length.out=4), cex.axis=1.4, las=2)
box()
title(main='Adult Count', cex.main=1.8, line=.5)
text(x=168, y=584, pos=4, label=paste('PPP = ', format(round(y_ppp, digits=2), nsmall=2), sep=''), cex=1.6)

r_sample <- array(c(rep(1,10000),rep(0,nsite*(nyear-1)*nreps*600-10000)), dim=c(nsite, nyear-1, nreps, 600))
plot(r_dsim[which(r_sample == 1)] ~ r_dobs[which(r_sample == 1)], 
     xlim=c(0,15000), ylim=c(0,15000), axes=F, xlab='', ylab='', pch=16, cex=1.2, col=col_point)
abline(0, 1, col=col_line, lwd=1.8)
axis(1, at=seq(0,15000,length.out=4), cex.axis=1.4)
axis(2, at=seq(0,15000,length.out=4), cex.axis=1.4, las=2)
box()
title(main='Juvenile Count', cex.main=1.8, line=.5)
text(x=4200, y=14600, pos=4, label=paste('PPP = ', format(round(r_ppp, digits=2), nsmall=2), sep=''), cex=1.6)

m_adu_sample <- array(c(rep(1,10000),rep(0,nsite*(nyear-1)*600-10000)), dim=c(nsite, nyear-1, 600))
plot(m_adu_dsim[which(m_adu_sample == 1)] ~ m_adu_dobs[which(m_adu_sample == 1)], 
     xlim=c(0,12), ylim=c(0,12), axes=F, xlab='', ylab='', pch=16, cex=1.2, col=col_point)
abline(0, 1, col=col_line, lwd=1.8)
axis(1, at=seq(0,12,length.out=4), cex.axis=1.4)
axis(2, at=seq(0,12,length.out=4), cex.axis=1.4, las=2)
box()
title(main='Adult Capture', cex.main=1.8, line=.5)
text(x=3.78, y=11.68, pos=4, label=paste('PPP = ', format(round(m_adu_ppp, digits=2), nsmall=2), sep=''), cex=1.6)

m_juv_sample <- array(c(rep(1,10000),rep(0,nsite*(nyear-1)*600-10000)), dim=c(nsite, nyear-1, 600))
plot(m_juv_dsim[which(m_juv_sample == 1)] ~ m_juv_dobs[which(m_juv_sample == 1)], 
     xlim=c(0,15), ylim=c(0,15), axes=F, xlab='', ylab='', pch=16, cex=1.2, col=col_point)
abline(0, 1, col=col_line, lwd=1.8)
axis(1, at=seq(0,15,length.out=4), cex.axis=1.4)
axis(2, at=seq(0,15,length.out=4), cex.axis=1.4, las=2)
box()
title(main='Juvenile Capture', cex.main=1.8, line=.5)
text(x=4.2, y=14.6, pos=4, label=paste('PPP = ', format(round(m_juv_ppp, digits=2), nsmall=2), sep=''), cex=1.6)

title(xlab=expression(D^"[obs]"), cex.lab=2, outer=T, line=1.6)
title(ylab=expression(D^"[sim]"), cex.lab=2, outer=T, line=0.8)

dev.off()


