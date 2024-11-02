#=================================================================================================
# Dynamic (i.e., multi-season) spatial capture-recapture (SCR) model with indiviudal genetic info
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#=================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 5_Spatial capture-recapture models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

set.seed(10)

# Basic values
nindv <- 100 # number of individuals (with data augmentation)
nsite <- 360 # number of sites
nyear <- 6   # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 4   # number of observational covariates

psi <- 0.4                                   # probability of first-year existence
kappa <- 2.5                                 # distance effect on habitat use
beta_upsilon <- c(2, 0.6, -0.3)              # intercept and slopes for habitat use
beta_phi   <- c( 1.0, -0.4, 0.6, -0.2,  0.8) # intercept and slopes for apparent survival
beta_gamma <- c(-1.0, -0.6, 0.8, -0.4, -0.2) # intercept and slopes for recruitment
mu_gene <- 2                                 # mean of genetic score
sigma_gene <- 0.5                            # standard deviation of genetic score
alpha <- c(-0.6, 0.4, -0.2, 0)               # intercept and slopes for detection probability
lon_min <- -10                               # lower bound of longitude or easting
lon_max <- +10                               # upper bound of longitude or easting
lat_min <- - 8                               # lower bound of latitude or northing
lat_max <- + 8                               # upper bound of latitude or northing

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

z <- matrix(0, nindv, nyear) # indicator of individual existence at a given year
e <- matrix(0, nindv, nyear) # indicator of individual existence so far
g <- rnorm(nindv, mu_gene, sigma_gene) # genetic score for each individual
phi   <- matrix(0, nindv, nyear-1) # apparent survival probability
gamma <- matrix(0, nindv, nyear-1) # recruitment probability
z[,1] <- rbinom(nindv, 1, psi) # individual existence in the first year
e[,1] <- z[,1]
for (v in 1:nindv) {
  for (t in 2:nyear) {
    phi  [v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g[v] - mu_gene) / sigma_gene) %*% beta_phi  )
    gamma[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g[v] - mu_gene) / sigma_gene) %*% beta_gamma)
    z[v,t] <- rbinom(1, 1, z[v,t-1] * phi[v,t-1] + (1 - e[v,t-1]) * gamma[v,t-1])
    e[v,t] <- ifelse(sum(z[v,1:t]) == 0, 0, 1)
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

q <- g # observed genetic score
q[which(apply(y,1,sum) == 0)] <- NA # genetic score is known only for individuals that are captured for at least once

#=======================
# Define MCMC algorithm
#=======================
gene_scr_mcmc <- function(y, q, x, w, o, ll, lon_min, lon_max, lat_min, lat_max, nmcmc) {

  # Setup variables
  nindv <- dim(y)[1]
  nsite <- dim(y)[2]
  nyear <- dim(y)[3]
  nreps <- dim(y)[4]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, 1:3, max)
  obs01 <- ifelse(apply(y, c(1,3), sum) > 0, 1, 0)
  beta_phi_update <- rep(0, nmcmc)
  beta_gamma_update <- rep(0, nmcmc)

  psi_save <- rep(0, nmcmc)
  kappa_save <- rep(0, nmcmc)
  beta_upsilon_save <- matrix(0, ncovs+1, nmcmc)
  beta_phi_save <- matrix(0, ncovs+3, nmcmc)
  beta_gamma_save <- matrix(0, ncovs+3, nmcmc)
  mu_gene_save <- rep(0, nmcmc)
  sigma_gene_save <- rep(0, nmcmc)
  alpha_save <- matrix(0, npcvs, nmcmc)
  g_save <- matrix(0, nindv, 2000)
  z_save <- array(, dim=c(nindv, nyear, 2000))

  # Priors
  log_kappa_mean <- 0
  log_kappa_sd <- 2
  beta_upsilon_mean <- 0
  beta_upsilon_sd <- 2
  beta_phi_mean <- 0
  beta_phi_sd <- 2
  beta_gamma_mean <- 0
  beta_gamma_sd <- 2
  mu_gene_mean <- 0
  mu_gene_sd <- 2
  sigma_gene_mean <- 2 # mean of sigma_gene square
  sigma_gene_sd <- 10 # standard deviation of sigma_gene square
  alpha_mean <- 0
  alpha_sd <- 2

  # Starting values
  psi <- 0.5
  kappa <- 1
  beta_upsilon <- rep(0, ncovs+1)
  beta_phi   <- rep(0, ncovs+3)
  beta_gamma <- rep(0, ncovs+3)
  mu_gene <- 0
  sigma_gene <- 2
  sigma_gene_a <- (sigma_gene_mean ^ 2) / (sigma_gene_sd ^ 2) + 2
  sigma_gene_b <- 1 / (sigma_gene_mean * (sigma_gene_a - 1))
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
  phi   <- matrix(0.5, nindv, nyear-1)
  gamma <- matrix(0.5, nindv, nyear-1)
  g <- q
  g[which(is.na(g))] <- rnorm(length(which(is.na(g))), mean(g,na.rm=T), sd(g,na.rm=T))
  p <- array(0.5, dim=c(nsite, nyear, nreps))
  zi <- obs01
  for (v in 1:nindv) {
    if (sum(zi[v,]) > 0) {
      zi[v,min(which(zi[v,] == 1)):max(which(zi[v,] == 1))] <- 1
    }
  } # v
  z <- zi
  e <- matrix(0, nindv, nyear)
  e[,1] <- z[,1]
  for (t in 2:nyear) {
    e[,t] <- ifelse(rowSums(z[,1:t])==0, 0, 1)
  } # t

  # Tuning factors
  kappa_tune <- 0.2
  beta_upsilon_tune <- 0.06
  beta_phi_tune <- 0.2
  beta_gamma_tune <- 0.4
  alpha_tune <- 0.035
  ss_tune <- 1
  g_tune <- 0.1

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
    mhK <- exp(mh1 - mh2)
    if (mhK > runif(1)) {
      kappa <- kappa_star
      upsilon<- upsilon_star
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
    beta_phi_star <- rnorm(ncovs+3, beta_phi, beta_phi_tune)
    phi_star <- matrix(0, nindv, nyear-1)
    for (v in 1:nindv) {
      for (t in 2:nyear) {
        phi_star[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g[v] - mu_gene) / sigma_gene) %*% beta_phi_star)
      } # t
    } # v
    prob1 <- prob2 <- rep(0, nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[t] <- sum(dbinom(z[which(z[,t]==1),t+1], 1, phi_star[which(z[,t]==1),t], log=TRUE))
      prob2[t] <- sum(dbinom(z[which(z[,t]==1),t+1], 1, phi     [which(z[,t]==1),t], log=TRUE))
    } # t
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_phi_star, beta_phi_mean, beta_phi_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta_phi     , beta_phi_mean, beta_phi_sd, log=TRUE))
    mhP <- exp(mh1 - mh2)
    if (mhP > runif(1)) {
      beta_phi <- beta_phi_star
      phi <- phi_star
      beta_phi_update[k] <- 1
    }
    if (k > 1000) {
      if (sum(beta_phi_update[(k-9):k]) == 0) {
        beta_phi <- runif(ncovs+3, min=apply(beta_phi_save[,(k-99):k],1,min), max=apply(beta_phi_save[,(k-99):k],1,max))
        phi <- matrix(0, nindv, nyear-1)
        for (v in 1:nindv) {
          for (t in 2:nyear) {
            phi[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g[v] - mu_gene) / sigma_gene) %*% beta_phi)
          } # t
        } # v
        beta_phi_update[k] <- 1
      }
    }

    ### Sample beta_gamma
    beta_gamma_star <- rnorm(ncovs+3, beta_gamma, beta_gamma_tune)
    gamma_star <- matrix(0, nindv, nyear-1)
    for (v in 1:nindv) {
      for (t in 2:nyear) {
        gamma_star[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g[v] - mu_gene) / sigma_gene) %*% beta_gamma_star)
      } # t
    } # v
    prob1 <- prob2 <- rep(0, nyear-1)
    for (t in 1:(nyear-1)) {
      prob1[t] <- sum(dbinom(z[which(e[,t]==0),t+1], 1, gamma_star[which(e[,t]==0),t], log=TRUE))
      prob2[t] <- sum(dbinom(z[which(e[,t]==0),t+1], 1, gamma     [which(e[,t]==0),t], log=TRUE))
    } # t
    mh1 <- sum(prob1) + 
           sum(dnorm(beta_gamma_star, beta_gamma_mean, beta_gamma_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(beta_gamma     , beta_gamma_mean, beta_gamma_sd, log=TRUE))
    mhP <- exp(mh1 - mh2)
    if (mhP > runif(1)) {
      beta_gamma <- beta_gamma_star
      gamma <- gamma_star
      beta_gamma_update[k] <- 1
    }
    if (k > 1000) {
      if (sum(beta_gamma_update[(k-9):k]) == 0) {
        beta_gamma <- runif(ncovs+3, min=apply(beta_gamma_save[,(k-99):k],1,min), max=apply(beta_gamma_save[,(k-99):k],1,max))
        gamma <- matrix(0, nindv, nyear-1)
        for (v in 1:nindv) {
          for (t in 2:nyear) {
            gamma[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g[v] - mu_gene) / sigma_gene) %*% beta_gamma)
         } # t
        } # v
        beta_gamma_update[k] <- 1
      }
    }

    ### Sample sigma_gene
    temp_a <- 1 / (sum((q - mu_gene) ^ 2, na.rm=T) / 2 + 1 / sigma_gene_b)
    temp_b <- length(which(!(is.na(q)))) / 2 + sigma_gene_a
    sigma_gene <- sqrt(1 / rgamma(1, shape=temp_b, scale=temp_a))

    ### Sample mu_gene
    temp_sd <- sqrt(1 / (length(which(!(is.na(q)))) / (sigma_gene ^ 2) + 1 / 4))
    temp_mean <- (temp_sd ^ 2) * (sum(q, na.rm=T) / (sigma_gene ^ 2))
    mu_gene <- rnorm(1, temp_mean, temp_sd)

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

    ### Sample g
    g_star <- q
    g_star[which(is.na(q))] <- rnorm(length(which(is.na(q))), g[which(is.na(q))], g_tune)
    phi_star <- gamma_star <- matrix(0, nindv, nyear-1)
    for (v in 1:nindv) {
      if (is.na(q[v])) {
        prob1 <- prob2 <- rep(0, nyear-1)
        for (t in 2:nyear) {
          phi_star  [v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g_star[v] - mu_gene) / sigma_gene) %*% beta_phi  )
          gamma_star[v,t-1] <- inv.logit(c(1, sum(upsilon[v,,t-1] * colSums(upsilon[-v,,t-1] * z[-v,t-1])) / sum(upsilon[v,,t-1]), colSums(upsilon[v,,t-1] * x[,t-1,]) / sum(upsilon[v,,t-1]), (g_star[v] - mu_gene) / sigma_gene) %*% beta_gamma)
          if (z[v,t-1] == 1) {
            prob1[t-1] <- dbinom(z[v,t], 1, phi_star[v,t-1], log=TRUE)
            prob2[t-1] <- dbinom(z[v,t], 1, phi     [v,t-1], log=TRUE)
          } else if (e[v,t-1] == 0) {
            prob1[t-1] <- dbinom(z[v,t], 1, gamma_star[v,t-1], log=TRUE)
            prob2[t-1] <- dbinom(z[v,t], 1, gamma     [v,t-1], log=TRUE)
          }
        } # t
        mh1 <- sum(prob1) + 
               dnorm(g_star[v], mu_gene, sigma_gene, log=TRUE)
        mh2 <- sum(prob2) + 
               dnorm(g     [v], mu_gene, sigma_gene, log=TRUE)
        mhg <- exp(mh1 - mh2)
        if (mhg > runif(1)) {
          g[v] <- g_star[v]
          phi[v,] <- phi_star[v,]
          gamma[v,] <- gamma_star[v,]
        }
      }
    } # v

    ### Sample z
    for (v in 1:nindv) {
      if (zi[v,1] == 0) {
        prod_temp <- psi * prod(exp(-1 * upsilon[v,,1] * p[,1,] * o))
        psi_temp <- prod_temp / (prod_temp + 1 - psi)
        psi_temp <- min(c(1, psi_temp))
        z[v,1] <- rbinom(1, 1, psi_temp)
      }
      e[v,1] <- z[v,1]
      for (t in 2:nyear) {
        if (e[v,t-1] == 1 & z[v,t-1] == 0 & zi[v,t] == 0) {
          z[v,t] <- 0
        } else if (z[v,t-1] == 1 & zi[v,t] == 0) {
          prod_temp <- phi[v,t-1] * prod(exp(-1 * upsilon[v,,t] * p[,t,] * o))
          psi_temp <- prod_temp / (prod_temp + 1 - phi[v,t-1])
          psi_temp <- min(c(1, psi_temp))
          z[v,t] <- rbinom(1, 1, psi_temp)
        } else if (e[v,t-1] == 0 & zi[v,t] == 0) {
          prod_temp <- gamma[v,t-1] * prod(exp(-1 * upsilon[v,,t] * p[,t,] * o))
          psi_temp <- prod_temp / (prod_temp + 1 - gamma[v,t-1])
          psi_temp <- min(c(1, psi_temp))
          z[v,t] <- rbinom(1, 1, psi_temp)
        }
        e[v,t] <- ifelse(sum(z[v,1:t]) == 0, 0, 1)
      } # t
    } # v

    ### Save samples
    psi_save[k] <- psi
    kappa_save[k] <- kappa
    beta_upsilon_save[,k] <- beta_upsilon
    beta_phi_save[,k] <- beta_phi
    beta_gamma_save[,k] <- beta_gamma
    mu_gene_save[k] <- mu_gene
    sigma_gene_save[k] <- sigma_gene
    alpha_save[,k] <- alpha
    if (nmcmc > 2000) {
      if (k > nmcmc-2000) {
        g_save[,k-nmcmc+2000] <- g
        z_save[,,k-nmcmc+2000] <- z
      }
    }
  } # k

  ### Write output
  list(psi_save=psi_save, 
       kappa_save=kappa_save, beta_upsilon_save=beta_upsilon_save, 
       beta_phi_save=beta_phi_save, beta_gamma_save=beta_gamma_save, 
       mu_gene_save=mu_gene_save, sigma_gene_save=sigma_gene_save, 
       alpha_save=alpha_save, 
       g_save=g_save, z_save=z_save)

} # gene_scr_mcmc

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
  gene_scr_mcmc(y=y, q=q, x=x, w=w, o=o, ll=ll, 
                lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, 
                nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='5.3 gene dynamic scr_output.RData')

#==============
# Plot results
#==============
pdf(file='5.3.1 gene dynamic scr_chains.pdf', width=10, height=8)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta_upsilon <- c(expression(beta[0]^""["["*upsilon*"]"]), expression(beta["bamboo"]^""["["*upsilon*"]"]), 
                       expression(beta["temperature"]^""["["*upsilon*"]"]))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["density"]^""["["*phi*"]"]), 
                   expression(beta["bamboo"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]), 
                   expression(beta["gene"]^""["["*phi*"]"]))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta["density"]^""["["*gamma*"]"]), 
                     expression(beta["bamboo"]^""["["*gamma*"]"]), expression(beta["temperature"]^""["["*gamma*"]"]), 
                     expression(beta["gene"]^""["["*gamma*"]"]))
ylab_alpha <- c(expression(alpha[1]), expression(alpha[2]), expression(alpha[3]), expression(alpha[4]))

par(mfrow=c(7,3))
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
text(x=nmcmc*0.42, y=psi + yint * 5, pos=4,  cex=1.2, 
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
text(x=nmcmc*0.42, y=kappa + yint * 1.25, pos=4,  cex=1.2, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

for (i in 1:(ncovs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_upsilon_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_upsilon[i] - yint * 2
  ymax <- beta_upsilon[i] + yint * 2
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_upsilon[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_upsilon[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_upsilon[i] + yint * 1.25, pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+3)) {
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
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_phi[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_phi[i] + yint * 5, pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+3)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_gamma_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma[i] - yint * c(16,8,8,8,16)[i]
  ymax <- beta_gamma[i] + yint * c(16,8,8,8,16)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*c(8,4,4,4,8)[i]), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=beta_gamma[i] + yint * c(10,5,5,5,10)[i], pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$mu_gene_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- mu_gene - yint * 8
ymax <- mu_gene + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=mu_gene, col='grey16', lwd=1.5)
title(main=expression(mu^"[gene]"), cex.main=2, line=1.2)
text(x=nmcmc*0.42, y=mu_gene + yint * 5, pos=4,  cex=1.2, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

tt <- matrix(, nmcmc, chain)
for (j in 1:chain) {
  tt[,j] <- out[[j]]$sigma_gene_save
} # j
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- sigma_gene - yint * 8
ymax <- sigma_gene + yint * 8
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.1, las=2)
box()
for (j in 1:chain) {
  lines(tt[,j], col=cols[j])
} # j
abline(h=sigma_gene, col='grey16', lwd=1.5)
title(main=expression(sigma^"[gene]"), cex.main=2, line=1.2)
text(x=nmcmc*0.42, y=sigma_gene + yint * 5, pos=4,  cex=1.2, 
     labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))

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
  axis(2, at=round(seq(ymin, ymax, yint*1), digits=1), cex.axis=1.1, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.42, y=alpha[i] + yint * 1.25, pos=4,  cex=1.2, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=1, outer=T)

dev.off()

pdf(file='5.3.2 gene dynamic scr_gene.pdf', width=10, height=8)

library(vioplot)

g_post <- cbind(out[[1]]$g_save, out[[2]]$g_save, out[[3]]$g_save)
g_med <- apply(g_post, 1, median)

par(mfrow=c(1,1))
par(mar=c(5.5,6,1,1))

plot(1, xlim=c(0.5,nyear+0.5), ylim=c(0.5,4.5), type='n', axes=F, xlab='', ylab='')
for (t in 1:nyear) {
  z_post <- cbind(out[[1]]$z_save[,t,], out[[2]]$z_save[,t,], out[[3]]$z_save[,t,])
  z_med <- apply(z_post, 1, median)
  vioplot(g[which(z[,t] == 1)], at=t, side='left', col='royalblue', rectCol='grey80', lineCol=NA, colMed=NA, border=NA, add=T)
  vioplot(g_med[which(z_med == 1)], at=t, side='right', col='lightcoral', rectCol='grey20', lineCol=NA, colMed=NA, border=NA, add=T)
  points(median(g[which(z[,t] == 1)]) ~ t, pch=16, col='yellow', cex=1.6)
  points(median(g_med[which(z_med == 1)]) ~ t, pch=16, col='tomato', cex=1.2)
} # v

axis(1, at=seq(1,nyear,1), cex.axis=1.8)
axis(2, at=seq(0,6,1), cex.axis=1.8, las=2)
box()
title(xlab='Year', cex.lab=2.5, line=4)
title(ylab='Genetic Score', cex.lab=2.5, line=3.8)

vioplot(rnorm(1000,4,0.1), at=2, side='left', col='royalblue', rectCol='grey80', lineCol=NA, colMed='yellow', pchMed=16, border=NA, cex=1.6, add=T)
vioplot(rnorm(1000,4,0.1), at=4, side='right', col='lightcoral', rectCol='grey20', lineCol=NA, colMed='tomato', pchMed=16, border=NA, cex=1.2, add=T)
text(x=2.1, y=4, pos=4, label='True', cex=1.8)
text(x=4.5, y=4, pos=4, label='Estimated', cex=1.8)

dev.off()


