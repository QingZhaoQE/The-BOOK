#========================================================================================================
# Spatial integrated population model (SIPM) with body-class-specific reproduction, survival and movement
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#========================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 4_Integrated population models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions
library(LaplacesDemon) # for generating random numbers from categorical distributions (rcat)

set.seed(7)

# Basic values
nsite <- 30 # number of sites
nyear <- 10 # number of years
nreps <- 5  # number of within-season replicates
ncovs <- 2  # number of environmental covariates
npcvs <- 1  # number of observational covariates

beta0 <- c(3.4, 1, -0.8)                  # intercept and slopes for initial abundance
beta_gamma_small <-  c(-1.0, -0.4, 0.5)   # intercept and slopes for reproduction of small individuals
beta_gamma_large <-  c(-0.6, -0.3, 0.4)   # intercept and slopes for reproduction of large individuals
beta_phi_small <- c(0.8, -0.7, 0.8, -0.4) # intercept and slopes for survival of small individuals
beta_phi_large <- c(1.0, -0.5, 0.6, -0.3) # intercept and slopes for survival of large individuals
omega <- 0.7                              # growth probability from small to large body class
kappa_small <- 0.5                        # distance effect on movement of small individuals
kappa_large <- 2                          # distance effect on movement of large individuals
alpha_det <- c( 0.9, -0.4)                # intercept and slopes for detection probability
alpha_cap <- c(-3.0, -0.8)                # intercept and slopes for capture probability

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
lambda0 <- exp(cbind(1,x[,1,]) %*% beta0) # expectation of initial abundance
N_small <- matrix(, nsite, nyear) # abundance of small individuals
N_large <- matrix(, nsite, nyear) # abundance of large individuals
N_small[,1] <- rpois(nsite, lambda0) # initial abundance of small individuals
N_large[,1] <- rpois(nsite, lambda0) # initial abundance of large individuals

lon <- cov$lon[1:nsite] # longitude or easting
lat <- cov$lat[1:nsite] # latitude or northing
d <- matrix(, nsite, nsite)  # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1
eta_small <- exp(-1 * kappa_small * d)        # unstandardized colonization rate of small individuals
theta_small <- eta_small / rowSums(eta_small) #   standardized colonization rate of small individuals
eta_large <- exp(-1 * kappa_large * d)        # unstandardized colonization rate of large individuals
theta_large <- eta_large / rowSums(eta_large) #   standardized colonization rate of large individuals

nindv <- 20000 # maximum number of individuals that ever exist
z <- matrix(0, nindv, nyear) # individual survival-location history
N_acc <- rep(0, nyear) # number of individuals that ever exist
for (i in 1:nsite) {
  if (i == 1) {
    tt <- c(rep(i, N_small[i,1]), rep(i+nsite, N_large[i,1]))
  } else {
    tt <- c(tt, rep(i, N_small[i,1]), rep(i+nsite, N_large[i,1]))
  }
}
z[1:length(tt),1] <- tt
N_acc[1] <- length(tt)

gamma_small <- matrix(, nsite, nyear-1) # reproduction rate or small individuals
gamma_large <- matrix(, nsite, nyear-1) # reproduction rate or large individuals
phi_small <- matrix(, nsite, nyear-1) # survival probability of small individuals
phi_large <- matrix(, nsite, nyear-1) # survival probability of large individuals
R_small <- matrix(, nsite, nyear-1) # number of individuals reproduced by small individuals
R_large <- matrix(, nsite, nyear-1) # number of individuals reproduced by large individuals

for (t in 2:nyear) {
  gamma_small[,t-1] <- exp(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,1]) %*% beta_gamma_small)
  gamma_large[,t-1] <- exp(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,1]) %*% beta_gamma_large)
  R_small[,t-1] <- rpois(nsite, N_small[,t-1] * gamma_small[,t-1])
  R_large[,t-1] <- rpois(nsite, N_large[,t-1] * gamma_large[,t-1])

  # All juveniles are small
  for (i in 1:nsite) {
    if (i == 1) { tt <- rep(i, R_small[i,t-1] + R_large[i,t-1]) } else { tt <- c(tt, rep(i, R_small[i,t-1] + R_large[i,t-1])) }
  } # i
  z[N_acc[t-1]+1:length(tt), t-1] <- tt

  phi_small[,t-1] <- inv.logit(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,]) %*% beta_phi_small)
  phi_large[,t-1] <- inv.logit(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,]) %*% beta_phi_large)
  transition_matrix <- matrix(, nsite*2, nsite*2+1) # matrix of transition (survival, growth & movement) probability
  for (i in 1:nsite) {
    transition_matrix[i,1:nsite] <- phi_small[i,t-1] * (1 - omega) * theta_small[i,]
    transition_matrix[i,nsite+1:nsite] <- phi_small[i,t-1] * omega * theta_large[i,]
    transition_matrix[nsite+i,1:nsite] <- 0
    transition_matrix[nsite+i,nsite+1:nsite] <- phi_large[i,t-1] * theta_large[i,]
  } # i
  transition_matrix[,nsite*2+1] <- 1 - rowSums(transition_matrix[,1:(nsite*2)])

  for (v in 1:dim(z)[1]) {
    if (z[v,t-1] == 0) {
      z[v,t] <- 0 
    } else {
      z[v,t] <- rcat(1, p=transition_matrix[z[v,t-1],])
    }
  } # v
  z[which(!(z %in% c(1:(nsite*2))))] <- 0

  for (i in 1:nsite) {
    N_small[i,t] <- length(which(z[,t] == i))
    N_large[i,t] <- length(which(z[,t] == nsite+i))
  }
  N_acc[t] <- max(which(rowMeans(z, na.rm=T) > 0))
} # t
z[which(!(z %in% 1:(nsite*2)))] <- 0
z <- z[which(rowMeans(z, na.rm=T) != 0),] # only keep individuals that ever exist

# Simulate count data
w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

pdet <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    pdet[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_det)
  } # t
} # i

y <- array(, dim=c(nsite, nyear, nreps, 2)) # count data of breeders
for (i in 1:nsite) {
  for (t in 1:nyear) {
    y[i,t,,1] <- rbinom(nreps, N_small[i,t], pdet[i,t,]) # count data of small breeders
    y[i,t,,2] <- rbinom(nreps, N_large[i,t], pdet[i,t,]) # count data of large breeders
  } # t
} # i

r <- array(, dim=c(nsite, nyear-1, nreps, 2)) # count data of offspring
for (i in 1:nsite) {
  for (t in 1:(nyear-1)) {
    r[i,t,,1] <- rbinom(nreps, R_small[i,t], pdet[i,t,]) # count data of offspring reproduced by small breeders
    r[i,t,,2] <- rbinom(nreps, R_large[i,t], pdet[i,t,]) # count data of offspring reproduced by large breeders
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
    t2 <- matrix(0, ncap, nyear)
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

# Converting individual capture history to m-array format
f_ext <- numeric(dim(c_ext)[1]) # extended year of first capture
for (v in 1:length(f_ext)) {
  f_ext[v] <- min(which(c_ext[v,] %in% c(1:(nsite*2))))
} # v

for (t in 2:nyear) {
  t1 <- paste('small_site', '_', 1:nsite, ', year', t, sep='')
  t2 <- paste('large_site', '_', 1:nsite, ', year', t, sep='')
  tt <- c(t1, t2)
  if (t == 2) {
    nm1 <- tt
  } else { 
    nm1 <- c(nm1, tt)
  }
} # t
nm1 <- c(nm1, 'not seen')
nm2 <- c(paste('small_site', 1:nsite, sep='_'), paste('large_site', 1:nsite, sep='_'))

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
cell_prob_fun <- function(nsite, nyear, phi_small, phi_large, omega, theta_small, theta_large, pcap) { 
  qcap <- 1 - pcap
  phi <- rbind(phi_small, phi_large)
  theta <- rbind(cbind((1-omega)*theta_small,omega*theta_large), cbind(matrix(0,nsite,nsite),theta_large))
  cell_prob_p <- cell_prob_q <- array(0, dim=c(nyear-1, (nsite*2)*(nyear-1)+1, (nsite*2)))
  for (i in 1:(nsite*2)) {
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,(nsite*2)*(t-1)+1:(nsite*2),i] <- phi[i,t] * theta[i,] * pcap[,t]
      cell_prob_q[t,(nsite*2)*(t-1)+1:(nsite*2),i] <- phi[i,t] * theta[i,] * qcap[,t]
    } # t
    for (t1 in 1:(nyear-2)) {
      for (t2 in (t1+1):(nyear-1)) {
        cell_prob_p[t1,(nsite*2)*(t2-1)+1:(nsite*2),i] <- 
          ((cell_prob_q[t1,(nsite*2)*(t2-2)+1:(nsite*2),i] * phi[,t2]) %*% theta) * pcap[,t2]
        cell_prob_q[t1,(nsite*2)*(t2-1)+1:(nsite*2),i] <- 
          ((cell_prob_q[t1,(nsite*2)*(t2-2)+1:(nsite*2),i] * phi[,t2]) %*% theta) * qcap[,t2]
      } # t2
    } # t1
    for (t in 1:(nyear-1)) {
      cell_prob_p[t,(nsite*2)*(nyear-1)+1,i] <- 1 - sum(cell_prob_p[t,1:((nsite*2)*(nyear-1)),i])
    } # t
  } # i
  cell_prob <- cell_prob_p

  return(cell_prob)
} # cell_prob_fun

body_sipm_mcmc <- function(y, r, m, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  y_small_max <- apply(y[,,,1], 1:2, max)
  y_large_max <- apply(y[,,,2], 1:2, max)

  beta0_save <- matrix(0, ncovs+1, nmcmc)
  beta_gamma_small_save <- matrix(0, ncovs+1, nmcmc)
  beta_gamma_large_save <- matrix(0, ncovs+1, nmcmc)
  beta_phi_small_save <- matrix(0, ncovs+2, nmcmc)
  beta_phi_large_save <- matrix(0, ncovs+2, nmcmc)
  omega_save <- rep(0, nmcmc)
  kappa_small_save <- rep(0, nmcmc)
  kappa_large_save <- rep(0, nmcmc)
  alpha_det_save <- matrix(0, npcvs+1, nmcmc)
  alpha_cap_save <- matrix(0, npcvs+1, nmcmc)
  N_small_save <- array(0, dim=c(nsite, nyear, 200))
  N_large_save <- array(0, dim=c(nsite, nyear, 200))

  # Prior
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_gamma_small_mean <- rep(0, ncovs+1)
  beta_gamma_small_sd <- 2
  beta_gamma_large_mean <- rep(0, ncovs+1)
  beta_gamma_large_sd <- 2
  beta_phi_small_mean <- rep(0, ncovs+2)
  beta_phi_small_sd <- 2
  beta_phi_large_mean <- rep(0, ncovs+2)
  beta_phi_large_sd <- 2
  logit_omega_mean <- 0
  logit_omega_sd <- 2
  log_kappa_small_mean <- 0
  log_kappa_small_sd <- 2
  log_kappa_large_mean <- 0
  log_kappa_large_sd <- 2
  alpha_det_mean <- rep(0, npcvs+1)
  alpha_det_sd <- 2
  alpha_cap_mean <- rep(0, npcvs+1)
  alpha_cap_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_gamma_small <- rep(0, ncovs+1)
  beta_gamma_large <- rep(0, ncovs+1)
  beta_phi_small <- rep(0, ncovs+2)
  beta_phi_large <- rep(0, ncovs+2)
  omega <- 0.5
  kappa_small <- 1
  kappa_large <- 1
  alpha_det <- rep(0, npcvs+1)
  alpha_cap <- rep(0, npcvs+1)
  lambda0 <- exp(cbind(1,x[,1,]) %*% beta0)
  gamma_small <- matrix(0.5, nsite, nyear-1)
  gamma_large <- matrix(0.5, nsite, nyear-1)
  phi_small <- matrix(0.5, nsite, nyear-1)
  phi_large <- matrix(0.5, nsite, nyear-1)
  eta_small <- exp(-1 * kappa_small * d)
  theta_small <- eta_small / rowSums(eta_small)
  eta_large <- exp(-1 * kappa_large * d)
  theta_large <- eta_large / rowSums(eta_large)
  N_small <- round((y_small_max + 1) / 0.5)
  N_large <- round((y_large_max + 1) / 0.5)
  pdet <- array(0.5, dim=c(nsite, nyear, nreps))
  pcap_per_visit <- array(0.5, dim=c(nsite*2, nyear-1, nreps))
  pcap <- matrix(, nsite*2, nyear-1)
  for (i in 1:(nsite*2)) {
    for (t in 1:(nyear-1)) {
      pcap[i,t] <- 1 - prod(1 - pcap_per_visit[i,t,])
    } # t
  } # i
  cell_prob <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, 
                             omega=omega, theta_small=theta_small, theta_large=theta_large, pcap=pcap)

  # Tuning factor
  beta0_tune <- 0.08
  beta_gamma_small_tune <- c(0.02, 0.02, 0.02)
  beta_gamma_large_tune <- c(0.02, 0.02, 0.02)
  beta_phi_small_tune <- c(0.10, 0.10, 0.07, 0.07)
  beta_phi_large_tune <- c(0.10, 0.10, 0.07, 0.07)
  omega_tune <- 0.05
  kappa_small_tune <- 0.035
  kappa_large_tune <- 0.035
  alpha_det_tune <- c(0.050, 0.035)
  alpha_cap_tune <- 0.050
  N_small_tune <- 1
  N_large_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N_small
    N_small_star <- matrix(rpois(nsite*nyear, N_small + N_small_tune), nsite, nyear)
    mh1 <- apply(dbinom(y[,,,1], N_small_star, pdet, log=T), 1:2, sum) + 
           dpois(N_small_star, cbind(lambda0, t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega)) %*% theta_small)), log=T) + 
           dpois(N_small, N_small_star+N_small_tune, log=T)
    mh2 <- apply(dbinom(y[,,,1], N_small     , pdet, log=T), 1:2, sum) + 
           dpois(N_small     , cbind(lambda0, t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega)) %*% theta_small)), log=T) + 
           dpois(N_small_star, N_small+N_small_tune, log=T)
    mh <- exp(mh1 - mh2)
    N_small_keep <- ((mh > runif(nsite*nyear)) & (N_small_star >= y_small_max))
    N_small[N_small_keep] <- N_small_star[N_small_keep]

    ### Sample N_large
    N_large_star <- matrix(rpois(nsite*nyear, N_large + N_large_tune), nsite, nyear)
    mh1 <- apply(dbinom(y[,,,2], N_large_star, pdet, log=T), 1:2, sum) + 
           dpois(N_large_star, cbind(lambda0, t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large)), log=T) + 
           dpois(N_large, N_large_star+N_large_tune, log=T)
    mh2 <- apply(dbinom(y[,,,2], N_large     , pdet, log=T), 1:2, sum) + 
           dpois(N_large     , cbind(lambda0, t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large)), log=T) + 
           dpois(N_large_star, N_large+N_large_tune, log=T)
    mh <- exp(mh1 - mh2)
    N_large_keep <- ((mh > runif(nsite*nyear)) & (N_large_star >= y_large_max))
    N_large[N_large_keep] <- N_large_star[N_large_keep]

    ### Sample beta0
    beta0_star <- rnorm(ncovs+1, beta0, beta0_tune)
    lambda0_star <- exp(cbind(1,x[,1,]) %*% beta0_star)
    mh1 <- sum(dpois(N_small[,1], lambda0_star, log=T)) + 
           sum(dpois(N_large[,1], lambda0_star, log=T)) + 
           sum(dnorm(beta0_star, beta0_mean, beta0_sd, log=T))
    mh2 <- sum(dpois(N_small[,1], lambda0     , log=T)) + 
           sum(dpois(N_large[,1], lambda0     , log=T)) + 
           sum(dnorm(beta0     , beta0_mean, beta0_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta0 <- beta0_star
      lambda0 <- lambda0_star
    }

    # Sample beta_gamma_small
    beta_gamma_small_star <- rnorm(ncovs+1, beta_gamma_small, beta_gamma_small_tune)
    gamma_small <- gamma_small_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      gamma_small     [,t-1] <- exp(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,1]) %*% beta_gamma_small     )
      gamma_small_star[,t-1] <- exp(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,1]) %*% beta_gamma_small_star)
    } # t
    mh1 <- sum(dbinom(r[,,,1], r[,,,1] + y[,-nyear,,1], gamma_small_star / (gamma_small_star + 1), log=T)) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small_star) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small_star) + N_large[,-nyear] * gamma_large) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large), log=T)) + 
           sum(dnorm(beta_gamma_small_star, beta_gamma_small_mean, beta_gamma_small_sd, log=T))
    mh2 <- sum(dbinom(r[,,,1], r[,,,1] + y[,-nyear,,1], gamma_small      / (gamma_small      + 1), log=T)) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small     ) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small     ) + N_large[,-nyear] * gamma_large) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large), log=T)) + 
           sum(dnorm(beta_gamma_small     , beta_gamma_small_mean, beta_gamma_small_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_gamma_small <- beta_gamma_small_star
    }

    # Sample beta_gamma_large
    beta_gamma_large_star <- rnorm(ncovs+1, beta_gamma_large, beta_gamma_large_tune)
    gamma_large <- gamma_large_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      gamma_large     [,t-1] <- exp(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,1]) %*% beta_gamma_large     )
      gamma_large_star[,t-1] <- exp(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,1]) %*% beta_gamma_large_star)
    } # t
    mh1 <- sum(dbinom(r[,,,2], r[,,,2] + y[,-nyear,,2], gamma_large_star / (gamma_large_star + 1), log=T)) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large_star) * phi_small * (1 - omega)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large_star) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large), log=T)) + 
           sum(dnorm(beta_gamma_large_star, beta_gamma_large_mean, beta_gamma_large_sd, log=T))
    mh2 <- sum(dbinom(r[,,,2], r[,,,2] + y[,-nyear,,2], gamma_large      / (gamma_large      + 1), log=T)) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large     ) * phi_small * (1 - omega)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large     ) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large), log=T)) + 
           sum(dnorm(beta_gamma_large     , beta_gamma_large_mean, beta_gamma_large_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_gamma_large <- beta_gamma_large_star
    }

    # Sample beta_phi_small & beta_phi_large
    beta_phi_small_star <- rnorm(ncovs+2, beta_phi_small, beta_phi_small_tune)
    beta_phi_large_star <- rnorm(ncovs+2, beta_phi_large, beta_phi_large_tune)
    phi_small <- phi_small_star <- phi_large <- phi_large_star <- matrix(, nsite, nyear-1)
    for (t in 2:nyear) {
      phi_small     [,t-1] <- inv.logit(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,]) %*% beta_phi_small     )
      phi_small_star[,t-1] <- inv.logit(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,]) %*% beta_phi_small_star)
      phi_large     [,t-1] <- inv.logit(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,]) %*% beta_phi_large     )
      phi_large_star[,t-1] <- inv.logit(cbind(1,(N_small[,t-1]+N_large[,t-1]-lambda0*2)/(lambda0*2),x[,t,]) %*% beta_phi_large_star)
    } # t
    cell_prob      <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small     , phi_large=phi_large     , omega=omega, theta_small=theta_small, theta_large=theta_large, pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small_star, phi_large=phi_large_star, omega=omega, theta_small=theta_small, theta_large=theta_large, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small_star * (1 - omega)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small_star * omega + N_large[,-nyear] * phi_large_star) %*% theta_large), log=T)) + 
           sum(dnorm(beta_phi_small_star, beta_phi_small_mean, beta_phi_small_sd, log=T)) + 
           sum(dnorm(beta_phi_large_star, beta_phi_large_mean, beta_phi_large_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small      * (1 - omega)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small      * omega + N_large[,-nyear] * phi_large     ) %*% theta_large), log=T)) + 
           sum(dnorm(beta_phi_small     , beta_phi_small_mean, beta_phi_small_sd, log=T)) + 
           sum(dnorm(beta_phi_large     , beta_phi_large_mean, beta_phi_large_sd, log=T))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      beta_phi_small <- beta_phi_small_star
      beta_phi_large <- beta_phi_large_star
      phi_small <- phi_small_star
      phi_large <- phi_large_star
    }

    ### Sample omega
    omega_star <- inv.logit(rnorm(1, logit(omega), omega_tune))
    cell_prob      <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, omega=omega     , theta_small=theta_small, theta_large=theta_large, pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, omega=omega_star, theta_small=theta_small, theta_large=theta_large, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega_star)) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * omega_star + N_large[,-nyear] * phi_large) %*% theta_large), log=T)) + 
           dnorm(logit(omega_star), logit_omega_mean, logit_omega_sd, log=T)
    mh2 <- sum(prob2) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega     )) %*% theta_small), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * omega      + N_large[,-nyear] * phi_large) %*% theta_large), log=T)) + 
           dnorm(logit(omega     ), logit_omega_mean, logit_omega_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      omega <- omega_star
    }

    # Sample kappa_small & kappa_large
    kappa_small_star <- exp(rnorm(1, log(kappa_small), kappa_small_tune))
    eta_small_star <- exp(-1 * kappa_small_star * d)
    theta_small_star <- eta_small_star / rowSums(eta_small_star)
    kappa_large_star <- exp(rnorm(1, log(kappa_large), kappa_large_tune))
    eta_large_star <- exp(-1 * kappa_large_star * d)
    theta_large_star <- eta_large_star / rowSums(eta_large_star)
    cell_prob      <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, omega=omega, theta_small=theta_small     , theta_large=theta_large     , pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, omega=omega, theta_small=theta_small_star, theta_large=theta_large_star, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        prob1[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob_star[t,,i], log=T)
        prob2[i,t] <- dmultinom(m[t,,i], size=sum(m[t,,i]), prob=cell_prob     [t,,i], log=T)
      } # t
    } # i
    mh1 <- sum(prob1) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega)) %*% theta_small_star), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large_star), log=T)) + 
           dnorm(log(kappa_small_star), log_kappa_small_mean, log_kappa_small_sd, log=T) + 
           dnorm(log(kappa_large_star), log_kappa_large_mean, log_kappa_large_sd, log=T)
    mh2 <- sum(prob2) + 
           sum(dpois(N_small[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * (1 - omega)) %*% theta_small     ), log=T)) + 
           sum(dpois(N_large[,-1], t(t((N_small[,-nyear] * (1 + gamma_small) + N_large[,-nyear] * gamma_large) * phi_small * omega + N_large[,-nyear] * phi_large) %*% theta_large     ), log=T)) + 
           dnorm(log(kappa_small     ), log_kappa_small_mean, log_kappa_small_sd, log=T) + 
           dnorm(log(kappa_large     ), log_kappa_large_mean, log_kappa_large_sd, log=T)
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      kappa_small <- kappa_small_star
      eta_small <- eta_small_star
      theta_small <- theta_small_star
      kappa_large <- kappa_large_star
      eta_large <- eta_large_star
      theta_large <- theta_large_star
    }

    ### Sample alpha_det
    alpha_det_star <- rnorm(npcvs+1, alpha_det, alpha_det_tune)
    pdet_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        pdet_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_det_star)
      } # t
    } # i
    mh1 <- sum(dbinom(y[,,,1], N_small, pdet_star, log=TRUE), na.rm=T) + 
           sum(dbinom(y[,,,2], N_large, pdet_star, log=TRUE), na.rm=T) + 
           sum(dnorm(alpha_det_star, alpha_det_mean, alpha_det_sd, log=TRUE))
    mh2 <- sum(dbinom(y[,,,1], N_small, pdet     , log=TRUE), na.rm=T) + 
           sum(dbinom(y[,,,2], N_large, pdet     , log=TRUE), na.rm=T) + 
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
    cell_prob      <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, omega=omega, theta_small=theta_small, theta_large=theta_large, pcap=pcap     )
    cell_prob_star <- cell_prob_fun(nsite=nsite, nyear=nyear, phi_small=phi_small, phi_large=phi_large, omega=omega, theta_small=theta_small, theta_large=theta_large, pcap=pcap_star)
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
    }

    # Save samples
    beta0_save[,k] <- beta0
    beta_gamma_small_save[,k] <- beta_gamma_small
    beta_gamma_large_save[,k] <- beta_gamma_large
    beta_phi_small_save[,k] <- beta_phi_small
    beta_phi_large_save[,k] <- beta_phi_large
    omega_save[k] <- omega
    kappa_small_save[k] <- kappa_small
    kappa_large_save[k] <- kappa_large
    alpha_det_save[,k] <- alpha_det
    alpha_cap_save[,k] <- alpha_cap
    if (nmcmc > 200) {
      if (k > nmcmc-200) {
        N_small_save[,,k-nmcmc+200] <- N_small
        N_large_save[,,k-nmcmc+200] <- N_large
      }
    }
  } # k

  list(beta0_save=beta0_save, 
       beta_gamma_small_save=beta_gamma_small_save, beta_gamma_large_save=beta_gamma_large_save, 
       beta_phi_small_save=beta_phi_small_save, beta_phi_large_save=beta_phi_large_save, 
       omega_save=omega_save, 
       kappa_small_save=kappa_small_save, kappa_large_save=kappa_large_save, 
       alpha_det_save=alpha_det_save, alpha_cap_save=alpha_cap_save, 
       N_small_save=N_small_save, N_large_save=N_large_save)
} # body_sipm_mcmc

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
  body_sipm_mcmc(y=y, r=r, m=m, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='4.3.3 body sipm_output.RData')

#==============
# Plot results
#==============
pdf(file='4.3.3.1 body sipm_chains.pdf', width=10, height=10)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["reef"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_gamma_small <- c(expression(beta["0,small"]^""["["*gamma*"]"]), expression(beta["density,small"]^""["["*gamma*"]"]), 
                           expression(beta["reef,small"]^""["["*gamma*"]"]))
ylab_beta_gamma_large <- c(expression(beta["0,large"]^""["["*gamma*"]"]), expression(beta["density,large"]^""["["*gamma*"]"]), 
                           expression(beta["reef,large"]^""["["*gamma*"]"]))
ylab_beta_phi_small <- c(expression(beta["0,small"]^""["["*phi*"]"]), expression(beta["density,small"]^""["["*phi*"]"]), 
                         expression(beta["reef,small"]^""["["*phi*"]"]), expression(beta["temperature,small"]^""["["*phi*"]"]))
ylab_beta_phi_large <- c(expression(beta["0,large"]^""["["*phi*"]"]), expression(beta["density,large"]^""["["*phi*"]"]), 
                         expression(beta["reef,large"]^""["["*phi*"]"]), expression(beta["temperature,large"]^""["["*phi*"]"]))
ylab_alpha_det <- c(expression(alpha[0]^"[det]"), expression(alpha[1]^"[det]"))
ylab_alpha_cap <- c(expression(alpha[0]^"[cap]"), expression(alpha[1]^"[cap]"))

par(mfrow=c(8,3))
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
  text(x=nmcmc*0.4, y=beta0[i]+yint*c(8,4,8)[i]*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_gamma_small_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma_small[i] - yint * c(4,8,4)[i]
  ymax <- beta_gamma_small[i] + yint * c(4,8,4)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(2,4,2)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma_small[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma_small[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_gamma_small[i]+yint*c(4,8,4)[i]*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_gamma_large_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma_large[i] - yint * 2
  ymax <- beta_gamma_large[i] + yint * 2
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma_large[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma_large[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_gamma_large[i]+yint*2*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_small_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi_small[i] - yint * 8
  ymax <- beta_phi_small[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi_small[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi_small[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_phi_small[i]+yint*8*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+2)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_phi_large_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi_large[i] - yint * 8
  ymax <- beta_phi_large[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi_large[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi_large[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_phi_large[i]+yint*8*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$omega_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- omega - yint * 2
ymax <- omega + yint * 2
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=omega, col='grey16', lwd=1.5)
title(main=expression(omega), cex.main=2, line=1.2)
text(x=nmcmc*0.4, y=omega+yint*2*0.6, pos=4,  cex=1.6, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$kappa_small_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa_small - yint * 4
ymax <- kappa_small + yint * 4
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa_small, col='grey16', lwd=1.5)
title(main=expression(kappa[small]), cex.main=2, line=1.2)
text(x=nmcmc*0.4, y=kappa_small+yint*4*0.6, pos=4,  cex=1.6, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$kappa_large_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- kappa_large - yint * 8
ymax <- kappa_large + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=kappa_large, col='grey16', lwd=1.5)
title(main=expression(kappa[large]), cex.main=2, line=1.2)
text(x=nmcmc*0.4, y=kappa_large+yint*8*0.6, pos=4,  cex=1.6, 
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
  text(x=nmcmc*0.4, y=alpha_det[i]+yint*c(4,2)[i]*0.6, pos=4,  cex=1.6, 
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
  text(x=nmcmc*0.4, y=alpha_cap[i]+yint*8*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=3, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=3, line=0.8, outer=T)

dev.off()

pdf(file='4.3.3.1 body sipm_N.pdf', width=9, height=8)

library(plotrix) # for making pie charts at selected positions

N_small_post <- N_large_post <- array(, dim=c(nsite, nyear, 200*chain))
N_small_post[,,  1:200] <- out[[1]]$N_small_save
N_small_post[,,201:400] <- out[[2]]$N_small_save
N_small_post[,,401:600] <- out[[3]]$N_small_save
N_large_post[,,  1:200] <- out[[1]]$N_large_save
N_large_post[,,201:400] <- out[[2]]$N_large_save
N_large_post[,,401:600] <- out[[3]]$N_large_save
N_small_med <- apply(N_small_post, 1:2, median)
N_large_med <- apply(N_large_post, 1:2, median)
years_plot <- c(1, 8)

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(4,5,5,2))

for (t in 1:2) {
  plot(1, xlim=c(-11, 11), ylim=c(-8.8, 8.8), type='n', axes=F, xlab='', ylab='')
  for (i in 1:nsite) {
    if (N_small[i,years_plot[t]] + N_large[i,years_plot[t]] > 0) {
      prop <- c(N_small[i,years_plot[t]], N_large[i,years_plot[t]]) / (N_small[i,years_plot[t]] + N_large[i,years_plot[t]])
      size <- (N_small[i,years_plot[t]] + N_large[i,years_plot[t]]) ^ 0.67 / 20
      floating.pie(xpos=cov$lon[i], ypos=cov$lat[i], x=prop, radius=size, col=c('white','royalblue'), border='grey18')
    }
  } # i
  axis(1, at=seq(-10,10,5), labels=rep('',5))
  if (t == 1) {
    axis(2, at=seq(-8,8,4), las=2, cex.axis=1.6)
  } else {
    axis(2, at=seq(-8,8,4), labels=rep('',5))
  }
  axis(3, at=0, labels=paste('Year', years_plot[t], sep=' '), cex.axis=2.5, line=-0.8, tick=F)
  if (t == 2) {
    axis(4, at=0, labels='True', cex.axis=2.5, line=0.5, tick=F)
  }
  box()
} # t

for (t in 1:2) {
  plot(1, xlim=c(-11, 11), ylim=c(-8.8, 8.8), type='n', axes=F, xlab='', ylab='')
  for (i in 1:nsite) {
    if (N_small_med[i,years_plot[t]] + N_large_med[i,years_plot[t]] > 0) {
      prop <- c(N_small_med[i,years_plot[t]], N_large_med[i,years_plot[t]]) / (N_small_med[i,years_plot[t]] + N_large_med[i,years_plot[t]])
      size <- (N_small_med[i,years_plot[t]] + N_large_med[i,years_plot[t]]) ^ 0.67 / 20
      floating.pie(xpos=cov$lon[i], ypos=cov$lat[i], x=prop, radius=size, col=c('white','lightcoral'), border='grey18')
    }
  } # i
  axis(1, at=seq(-10,10,5), cex.axis=1.6)
  if (t == 1) {
    axis(2, at=seq(-8,8,4), las=2, cex.axis=1.6)
  } else {
    axis(2, at=seq(-8,8,4), labels=rep('',5))
  }
  if (t == 2) {
    axis(4, at=0, labels='Estimated', cex.axis=2.5, line=0.5, tick=F)
  }
  box()
} # t

title(xlab='Easting' , cex.lab=2.8, line=2.6, outer=T)
title(ylab='Northing', cex.lab=2.8, line=2.4, outer=T)
title(main=expression('Proportion of Body Class'), line=3.2, cex.main=3, outer=T)

dev.off()


