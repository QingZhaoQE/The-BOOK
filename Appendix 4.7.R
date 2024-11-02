#========================================================================================================
# Spatial integrated population model (SIPM) with gene-class-specific reproduction, survival and movement
# code for simulating data, defining MCMC algorithm, implementing the model, and creating figures
# written by Qing Zhao, 2023 in Colorado
#========================================================================================================

setwd('c:/Zhao/RESEARCH/C. model/a. SDNM/b. The BOOK/Chapter 4_Integrated population models/')

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions
library(LaplacesDemon) # for generating random numbers from categorical distributions (rcat)

set.seed(8)

# Basic values
nsite <- 30 # number of sites
nyear <- 10 # number of years
nreps <- 5  # number of within-season replicates
ngene <- 6  # number of classes of genetic variation
ncovs <- 2  # number of environmental covariates
npcvs <- 1  # number of observational covariates

beta0 <- c(3.4, 1, -0.8)                          # intercept and slopes for initial abundance
beta_gamma <- c(-1.2, -0.3, 0.4, -0.2, 0.4, -0.3) # intercept and slopes for reproduction
beta_phi   <- c( 1.0, -0.5, 0.6, -0.3, 0.6, -0.4) # intercept and slopes for apparent survival
beta_kappa <- c(-0.8, -0.6, 0.6, 0.8, -0.4)       # intercept and slopes for movement
alpha_det <- c( 1  , -0.4)                        # intercept and slopes for detection probability
alpha_cap <- c(-1.4, -0.8)                        # intercept and slopes for capture probability

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
N <- array(, dim=c(nsite, nyear, ngene)) # abundance for each site, year and gene class
lambda0 <- exp(cbind(1,x[,1,]) %*% beta0) # expectation of initial abundance
for (i in 1:nsite) {
  N[i,1,] <- rpois(ngene, lambda0[i]/ngene) # initial abundance
} # i

lon <- cov$lon[1:nsite] # longitude or easting
lat <- cov$lat[1:nsite] # latitude or northing
d <- matrix(, nsite, nsite)  # distance matrix
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1

nindv <- 15000 # maximum number of individuals that ever exist
z <- matrix(0, nindv, nyear) # individual survival-location history
gene <- numeric(nindv) # gene class for each individual
for (i in 1:nsite) {
  for (g in 1:ngene) {
    if (i == 1 & g == 1) {
      tz <- rep(i, N[i,1,g])
      tg <- rep(g, N[i,1,g])
    } else { 
      tz <- c(tz, rep(i, N[i,1,g]))
      tg <- c(tg, rep(g, N[i,1,g]))
    }
  } # g
} # i
z[1:length(tz),1] <- tz
gene[1:length(tg)] <- tg

gamma <- array(, dim=c(nsite, nyear-1, ngene)) # reproduction rate
phi   <- array(, dim=c(nsite, nyear-1, ngene)) # survival probability
kappa <- array(, dim=c(nsite, nyear-1, ngene)) # distance effect on movement
eta   <- array(, dim=c(nsite, nsite, nyear-1, ngene)) # unstandardized movement rate
theta <- array(, dim=c(nsite, nsite, nyear-1, ngene)) #   standardized movement rate
R     <- array(, dim=c(nsite, nyear-1, ngene)) # number of reproduced individuals

for (t in 2:nyear) {
  for (i in 1:nsite) {
    for (g in 1:ngene) {
      x_temp1 <- (sum(N[i,t-1,])-lambda0[i])/lambda0[i]
      x_temp2 <- x[i,t,1]
      x_temp3 <- x[i,t,2]
      x_temp4 <- (g - 3.5) / 1.7 # standardize gene class with hypothesized mean (3.5) and SD (1.7)
      x_temp5 <- x_temp2 * x_temp4
      x_temp <- c(1, x_temp1, x_temp2, x_temp3, x_temp4, x_temp5)
      gamma[i,t-1,g] <- exp(x_temp %*% beta_gamma)
      phi[i,t-1,g] <- inv.logit(x_temp %*% beta_phi)
      kappa[i,t-1,g] <- exp(x_temp[-4] %*% beta_kappa)
      eta[i,,t-1,g] <- exp(-1 * kappa[i,t-1,g] * d[i,])
      theta[i,,t-1,g] <- eta[i,,t-1,g] / sum(eta[i,,t-1,g])
      R[i,t-1,g] <- rpois(1, N[i,t-1,g] * gamma[i,t-1,g])
      if (i == 1 & g == 1) {
        tz <- rep(i, R[i,t-1,g])
        tg <- rep(g, R[i,t-1,g])
      } else { 
        tz <- c(tz, rep(i, R[i,t-1,g]))
        tg <- c(tg, rep(g, R[i,t-1,g]))
      }
    } # g
  } # i
  z[max(which(rowMeans(z,na.rm=T)>0))+c(1:length(tz)), t-1] <- tz
  gene[max(which(gene>0))+c(1:length(tg))] <- tg

  for (v in 1:nindv) {
    if (z[v,t-1] %in% 1:nsite) {
      z[v,t] <- rcat(1, c(phi[z[v,t-1],t-1,gene[v]] * theta[z[v,t-1],,t-1,gene[v]], 1 - phi[z[v,t-1],t-1,gene[v]]))
    }
  } # v
  z[which(!(z %in% 1:nsite))] <- 0

  for (i in 1:nsite) {
    for (g in 1:ngene) {
      N[i,t,g] <- length(which(z[,t] == i & gene == g))
    } # g
  } # i
} # t

exit <- which(rowMeans(z, na.rm=T) > 0) # indicating if an individual ever exists
z <- z[exit,]      # only keep individuals that ever exist
gene <- gene[exit] # only keep individuals that ever exist

# Simulate count data
w <- array(rnorm(nsite * nyear * nreps * npcvs, 0, 1), dim=c(nsite, nyear, nreps, npcvs)) # observational covariates

pdet <- array(, dim=c(nsite, nyear, nreps)) # detection probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    pdet[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_det)
  } # t
} # i

y <- array(, dim=c(nsite, nyear, nreps, ngene)) # count data of breeders
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (g in 1:ngene) {
      y[i,t,,g] <- rbinom(nreps, N[i,t,g], pdet[i,t,])
    } # g
  } # t
} # i

r <- array(, dim=c(nsite, nyear-1, nreps, ngene)) # count data of offspring
for (i in 1:nsite) {
  for (t in 1:(nyear-1)) {
    for (g in 1:ngene) {
      r[i,t,,g] <- rbinom(nreps, R[i,t,g], pdet[i,t,])
    } # g
  } # t
} # i

# simulate capture-recapture data
pcap_per_visit <- array(, dim=c(nsite, nyear, nreps)) # per visit capture probability
for (i in 1:nsite) {
  for (t in 1:nyear) {
    pcap_per_visit[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_cap)
  } # t
} # i

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
  f[v] <- min(which(c[v,] > 0))
} # v
c <- c[which(f %in% c(1:(nyear-1))),]

gene_cap <- gene[which(if_cap==1)]
gene_cap <- gene_cap[which(f %in% c(1:(nyear-1)))]

# Extend individual capture history
for (v in 1:dim(c)[1]) {
  t1 <- c[v,]
  t1_gene <- gene_cap[v]
  ycap <- which(t1 %in% c(1:nsite))
  ncap <- length(ycap)
  if (ncap == 1) {
    t2 <- t1
    t2_gene <- t1_gene 
  } else {
    t2 <- matrix(, ncap, nyear)
    for (k in 1:(ncap-1)) {
      t2[k, ycap[k]:ycap[k+1]] <- t1[ycap[k]:ycap[k+1]]
    } # k
    t2[ncap, ycap[ncap]] <- t1[ycap[ncap]]
    t2_gene <- rep(t1_gene, ncap)
  }
  if (v == 1) {
    c_ext <- t2
    gene_ext <- t2_gene
  } else {
    c_ext <- rbind(c_ext, t2)
    gene_ext <- c(gene_ext, t2_gene)
  }
} # v
c_ext[which(!(c_ext %in% 1:nsite))] <- 0

# Converting individual capture history to m-array format
f_ext <- numeric(dim(c_ext)[1])
for (v in 1:length(f_ext)) {
  f_ext[v] <- min(which(c_ext[v,] %in% c(1:nsite)))
} # v

nm <- matrix(, ngene, nsite)
for (k in 1:ngene) {
  for (i in 1:nsite) {
    nm[k,i] <- paste('site ', i, ', gene ', k, sep='')
  } # i
} # k
nm <- as.vector(nm)

m <- array(0, # truncated m-array of capture-recapture (i.e., next-year recapture only)
           dim=c(nsite*ngene, nsite+1, nyear-1), 
           dimnames=list(nm, c(paste('site',1:nsite,sep=' '),'not seen'), paste('year',1:(nyear-1),sep=' ')))
for (t in 1:(nyear-1)) {
  t1c <- c_ext[which(f_ext == t),]
  t1g <- gene_ext[which(f_ext == t)]
  for (i in 1:nsite) {
    for (g in 1:ngene) {
      nband <- length(which(t1c[,t] == i & t1g == g))
      if (nband > 0) {
        t2 <- matrix(t1c[which(t1c[,t] == i & t1g == g),], nrow=nband, ncol=nyear)
        for (j in 1:nsite) {
          m[ngene*(i-1)+g,j,t] <- length(which(t2[,t+1] == j))
        } # j
        m[ngene*(i-1)+g,nsite+1,t] <- nband - sum(m[ngene*(i-1)+g,1:nsite,t])
      }
    } # g
  } # i
} # t

#=======================
# Define MCMC algorithm
#=======================
# Function to calculate the cell probability of truncated m-array
cell_prob_fun <- function(nsite, ngene, nyear, phi, theta, pcap) {
  cell_prob <- array(0, dim=c(nsite*ngene, nsite+1, nyear-1))
  for (t in 1:(nyear-1)) {
    for (i in 1:nsite) {
      for (g in 1:ngene) {
        cell_prob[ngene*(i-1)+g,1:nsite,t] <- phi[i,t,g] * theta[i,,t,g] * pcap[,t]
        cell_prob[ngene*(i-1)+g,nsite+1,t] <- 1 - sum(cell_prob[ngene*(i-1)+g,1:nsite,t])
      } # g
    } # i
  } # t

  return(cell_prob)
} # cell_prob_fun

gene_sipm_mcmc <- function(y, r, m, x, w, d, nmcmc) {

  # Setup variables
  nsite <- dim(y)[1]
  nyear <- dim(y)[2]
  nreps <- dim(y)[3]
  ngene <- dim(y)[4]
  ncovs <- dim(x)[3]
  npcvs <- dim(w)[4]
  ymax <- apply(y, c(1,2,4), max)
  xgene <- (1:ngene - 3.5) / 1.7

  beta0_save <- matrix(0, ncovs+1, nmcmc)
  beta_gamma_save <- matrix(0, ncovs+4, nmcmc)
  beta_phi_save <- matrix(0, ncovs+4, nmcmc)
  beta_kappa_save <- matrix(0, ncovs+3, nmcmc)
  alpha_det_save <- matrix(0, npcvs+1, nmcmc)
  alpha_cap_save <- matrix(0, npcvs+1, nmcmc)
  N_save <- array(0, dim=c(nsite, nyear, ngene, 500))

  # Prior
  beta0_mean <- rep(0, ncovs+1)
  beta0_sd <- 2
  beta_gamma_mean <- rep(0, ncovs+4)
  beta_gamma_sd <- 2
  beta_phi_mean <- rep(0, ncovs+4)
  beta_phi_sd <- 2
  beta_kappa_mean <- rep(0, ncovs+3)
  beta_kappa_sd <- 2
  alpha_det_mean <- rep(0, npcvs+1)
  alpha_det_sd <- 2
  alpha_cap_mean <- rep(0, npcvs+1)
  alpha_cap_sd <- 2

  # Starting values
  beta0 <- rep(0, ncovs+1)
  beta_gamma <- rep(0, ncovs+4)
  beta_phi   <- rep(0, ncovs+4)
  beta_kappa <- rep(0, ncovs+3)
  alpha_det <- rep(0, npcvs+1)
  alpha_cap <- rep(0, npcvs+1)
  lambda0 <- rep(1, nsite)
  gamma <- array(1, dim=c(nsite, nyear-1, ngene))
  phi <- array(0.5, dim=c(nsite, nyear-1, ngene))
  kappa <- array(1, dim=c(nsite, nyear-1, ngene))
  eta <- theta <- array(, dim=c(nsite, nsite, nyear-1, ngene))
  for (i in 1:nsite) {
    for (t in 1:(nyear-1)) {
      for (g in 1:ngene) {
        eta[i,,t,g] <- exp(-1 * kappa[i,t,g] * d[i,])
        theta[i,,t,g] <- eta[i,,t,g] / sum(eta[i,,t,g])
      } # g
    } # t
  } # i
  N <- ymax + 1
  pdet <- array(0.5, dim=c(nsite, nyear, nreps))
  pcap_per_visit <- array(0.5, dim=c(nsite, nyear-1, nreps))
  pcap <- matrix(, nsite, nyear-1)
  for (i in 1:nsite) {
    for (t in 1:(nyear-1)) {
      pcap[i,t] <- 1 - prod(1 - pcap_per_visit[i,t,])
    } # t
  } # i
  cell_prob <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi, theta=theta, pcap=pcap)

  # Tuning factor
  beta0_tune <- 0.08
  beta_gamma_tune <- 0.02
  beta_phi_tune <- 0.08
  beta_kappa_tune <- 0.04
  alpha_det_tune <- 0.025
  alpha_cap_tune <- 0.05
  N_tune <- 1

  # Begin MCMC loop
  for (k in 1:nmcmc) {
    ### Sample N
    N_star <- array(rpois(nsite*nyear*ngene, N + N_tune), dim=c(nsite, nyear, ngene))
    prob1 <- prob2 <- array(, dim=c(nsite, nyear, nreps, ngene))
    for (j in 1:nreps) {
      for (g in 1:ngene) {
        prob1[,,j,g] <- dbinom(y[,,j,g], N_star[,,g], pdet[,,j], log=T)
        prob2[,,j,g] <- dbinom(y[,,j,g], N     [,,g], pdet[,,j], log=T)
      } # g
    } # j
    Nexp <- array(, dim=c(nsite, nyear, ngene))
    Nexp[,1,] <- lambda0 / ngene
    for (t in 2:nyear) {
      for (g in 1:ngene) {
        Nexp[,t,g] <- (N[,t-1,g] * (1 + gamma[,t-1,g]) * phi[,t-1,g]) %*% theta[,,t-1,g] + 1e-6
      } # g
    } # t
    mh1 <- apply(prob1, c(1,2,4), sum) + 
           dpois(N_star, Nexp, log=T) + 
           dpois(N, N_star+N_tune, log=T)
    mh2 <- apply(prob2, c(1,2,4), sum) + 
           dpois(N     , Nexp, log=T) + 
           dpois(N_star, N+N_tune, log=T)
    mhN <- exp(mh1 - mh2)
    Nkeep <- ((mhN > runif(nsite*nyear*ngene)) & (N_star >= ymax))
    N[Nkeep] <- N_star[Nkeep]

    ### Sample beta0
    beta0_star <- rnorm(ncovs+1, beta0, beta0_tune)
    lambda0_star <- exp(cbind(1,x[,1,]) %*% beta0_star)
    mh1 <- sum(dpois(N[,1,], lambda0_star / ngene, log=T)) + 
           sum(dnorm(beta0_star, beta0_mean, beta0_sd, log=T))
    mh2 <- sum(dpois(N[,1,], lambda0      / ngene, log=T)) + 
           sum(dnorm(beta0     , beta0_mean, beta0_sd, log=T))
    mhB0 <- exp(mh1 - mh2)
    if (mhB0 > runif(1)) {
      beta0 <- beta0_star
      lambda0 <- lambda0_star
    }

    # Sample beta_gamma
    beta_gamma_star <- rnorm(ncovs+4, beta_gamma, beta_gamma_tune)
    gamma <- gamma_star <- array(, dim=c(nsite, nyear-1, ngene))
    for (t in 2:nyear) {
      for (g in 1:ngene) {
        x_temp1 <- (rowSums(N[,t-1,]) - lambda0) / lambda0
        x_temp2 <- x[,t,1]
        x_temp3 <- x[,t,2]
        x_temp4 <- (g - 3.5) / 1.7
        x_temp5 <- x_temp2 * x_temp4
        x_temp <- cbind(1, x_temp1, x_temp2, x_temp3, x_temp4, x_temp5)
        gamma     [,t-1,g] <- exp(x_temp %*% beta_gamma     )
        gamma_star[,t-1,g] <- exp(x_temp %*% beta_gamma_star)
      } # g
    } # t
    Nexp <- Nexp_star <- array(, dim=c(nsite, nyear-1, ngene))
    for (t in 1:(nyear-1)) {
      for (g in 1:ngene) {
        Nexp     [,t,g] <- (N[,t,g] * (1 + gamma     [,t,g]) * phi[,t,g]) %*% theta[,,t,g] + 1e-6
        Nexp_star[,t,g] <- (N[,t,g] * (1 + gamma_star[,t,g]) * phi[,t,g]) %*% theta[,,t,g] + 1e-6
      } # g
    } # t
    prob1 <- prob2 <- numeric(nreps)
    for (j in 1:nreps) {
      prob1[j] <- sum(dbinom(r[,,j,], r[,,j,] + y[,-nyear,j,], gamma_star / (gamma_star + 1), log=T))
      prob2[j] <- sum(dbinom(r[,,j,], r[,,j,] + y[,-nyear,j,], gamma      / (gamma      + 1), log=T))
    } # j
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1,], Nexp_star, log=T)) + 
           sum(dnorm(beta_gamma_star, beta_gamma_mean, beta_gamma_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1,], Nexp     , log=T)) + 
           sum(dnorm(beta_gamma     , beta_gamma_mean, beta_gamma_sd, log=T))
    mhBG <- exp(mh1 - mh2)
    if (mhBG > runif(1)) {
      beta_gamma <- beta_gamma_star
      gamma <- gamma_star
    }

    # Sample beta_phi
    beta_phi_star <- rnorm(ncovs+4, beta_phi, beta_phi_tune)
    phi <- phi_star <- array(, dim=c(nsite, nyear-1, ngene))
    for (t in 2:nyear) {
      for (g in 1:ngene) {
        x_temp1 <- (rowSums(N[,t-1,]) - lambda0) / lambda0
        x_temp2 <- x[,t,1]
        x_temp3 <- x[,t,2]
        x_temp4 <- (g - 3.5) / 1.7
        x_temp5 <- x_temp2 * x_temp4
        x_temp <- cbind(1, x_temp1, x_temp2, x_temp3, x_temp4, x_temp5)
        phi     [,t-1,g] <- inv.logit(x_temp %*% beta_phi     )
        phi_star[,t-1,g] <- inv.logit(x_temp %*% beta_phi_star)
      } # g
    } # t
    cell_prob      <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi     , theta=theta, pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi_star, theta=theta, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite*ngene, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        for (g in 1:ngene) {
          prob1[ngene*(i-1)+g,t] <- dmultinom(m[ngene*(i-1)+g,,t], size=sum(m[ngene*(i-1)+g,,t]), prob=cell_prob_star[ngene*(i-1)+g,,t], log=T)
          prob2[ngene*(i-1)+g,t] <- dmultinom(m[ngene*(i-1)+g,,t], size=sum(m[ngene*(i-1)+g,,t]), prob=cell_prob     [ngene*(i-1)+g,,t], log=T)
        } # g
      } # t
    } # i
    Nexp <- Nexp_star <- array(, dim=c(nsite, nyear-1, ngene))
    for (t in 1:(nyear-1)) {
      for (g in 1:ngene) {
        Nexp     [,t,g] <- (N[,t,g] * (1 + gamma[,t,g]) * phi     [,t,g]) %*% theta[,,t,g] + 1e-6
        Nexp_star[,t,g] <- (N[,t,g] * (1 + gamma[,t,g]) * phi_star[,t,g]) %*% theta[,,t,g] + 1e-6
      } # g
    } # t
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1,], Nexp_star, log=T)) + 
           sum(dnorm(beta_phi_star, beta_phi_mean, beta_phi_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1,], Nexp     , log=T)) + 
           sum(dnorm(beta_phi     , beta_phi_mean, beta_phi_sd, log=T))
    mhBO <- exp(mh1 - mh2)
    if (mhBO > runif(1)) {
      beta_phi <- beta_phi_star
      phi <- phi_star
    }

    # Sample beta_kappa
    beta_kappa_star <- rnorm(ncovs+3, beta_kappa, beta_kappa_tune)
    kappa <- kappa_star <- array(, dim=c(nsite, nyear-1, ngene))
    eta <- theta <- eta_star <- theta_star <- array(, dim=c(nsite, nsite, nyear-1, ngene))
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        for (g in 1:ngene) {
          x_temp1 <- (sum(N[i,t,]) - lambda0[i]) / lambda0[i]
          x_temp2 <- x[i,t+1,1]
          x_temp3 <- (g - 3.5) / 1.7
          x_temp4 <- x_temp2 * x_temp3
          x_temp <- c(1, x_temp1, x_temp2, x_temp3, x_temp4)
          kappa     [i,t,g] <- exp(x_temp %*% beta_kappa     )
          kappa_star[i,t,g] <- exp(x_temp %*% beta_kappa_star)
          eta     [i,,t,g] <- exp(-1 * kappa     [i,t,g] * d[i,])
          eta_star[i,,t,g] <- exp(-1 * kappa_star[i,t,g] * d[i,])
          theta     [i,,t,g] <- eta     [i,,t,g] / sum(eta     [i,,t,g])
          theta_star[i,,t,g] <- eta_star[i,,t,g] / sum(eta_star[i,,t,g])
        } # g
      } # t
    } # i
    cell_prob      <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi, theta=theta     , pcap=pcap)
    cell_prob_star <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi, theta=theta_star, pcap=pcap)
    prob1 <- prob2 <- matrix(, nsite*ngene, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        for (g in 1:ngene) {
          prob1[ngene*(i-1)+g,t] <- dmultinom(m[ngene*(i-1)+g,,t], size=sum(m[ngene*(i-1)+g,,t]), prob=cell_prob_star[ngene*(i-1)+g,,t], log=T)
          prob2[ngene*(i-1)+g,t] <- dmultinom(m[ngene*(i-1)+g,,t], size=sum(m[ngene*(i-1)+g,,t]), prob=cell_prob     [ngene*(i-1)+g,,t], log=T)
        } # g
      } # t
    } # i
    Nexp <- Nexp_star <- array(, dim=c(nsite, nyear-1, ngene))
    for (t in 1:(nyear-1)) {
      for (g in 1:ngene) {
        Nexp     [,t,g] <- (N[,t,g] * (1 + gamma[,t,g]) * phi[,t,g]) %*% theta     [,,t,g] + 1e-6
        Nexp_star[,t,g] <- (N[,t,g] * (1 + gamma[,t,g]) * phi[,t,g]) %*% theta_star[,,t,g] + 1e-6
      } # g
    } # t
    mh1 <- sum(prob1) + 
           sum(dpois(N[,-1,], Nexp_star, log=T)) + 
           sum(dnorm(beta_kappa_star, beta_kappa_mean, beta_kappa_sd, log=T))
    mh2 <- sum(prob2) + 
           sum(dpois(N[,-1,], Nexp     , log=T)) + 
           sum(dnorm(beta_kappa     , beta_kappa_mean, beta_kappa_sd, log=T))
    mhBK <- exp(mh1 - mh2)
    if (mhBK > runif(1)) {
      beta_kappa <- beta_kappa_star
      kappa <- kappa_star
      eta <- eta_star
      theta <- theta_star
    }

    ### Sample alpha_det
    alpha_det_star <- rnorm(npcvs+1, alpha_det, alpha_det_tune)
    pdet_star <- array(, dim=c(nsite, nyear, nreps))
    for (i in 1:nsite) {
      for (t in 1:nyear) {
        pdet_star[i,t,] <- inv.logit(cbind(1,w[i,t,,]) %*% alpha_det_star)
      } # t
    } # i
    prob1 <- prob2 <- array(, dim=c(nsite, nyear, nreps, ngene))
    for (j in 1:nreps) {
      for (g in 1:ngene) {
        prob1[,,j,g] <- dbinom(y[,,j,g], N[,,g], pdet_star[,,j], log=T)
        prob2[,,j,g] <- dbinom(y[,,j,g], N[,,g], pdet     [,,j], log=T)
      } # g
    } # j
    mh1 <- sum(prob1) + 
           sum(dnorm(alpha_det_star, alpha_det_mean, alpha_det_sd, log=TRUE))
    mh2 <- sum(prob2) + 
           sum(dnorm(alpha_det     , alpha_det_mean, alpha_det_sd, log=TRUE))
    mhAD <- exp(mh1 - mh2)
    if (mhAD > runif(1)) {
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
    cell_prob      <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi, theta=theta, pcap=pcap     )
    cell_prob_star <- cell_prob_fun(nsite=nsite, ngene=ngene, nyear=nyear, phi=phi, theta=theta, pcap=pcap_star)
    prob1 <- prob2 <- matrix(, nsite*ngene, nyear-1)
    for (i in 1:nsite) {
      for (t in 1:(nyear-1)) {
        for (g in 1:ngene) {
          prob1[ngene*(i-1)+g,t] <- dmultinom(m[ngene*(i-1)+g,,t], size=sum(m[ngene*(i-1)+g,,t]), prob=cell_prob_star[ngene*(i-1)+g,,t], log=T)
          prob2[ngene*(i-1)+g,t] <- dmultinom(m[ngene*(i-1)+g,,t], size=sum(m[ngene*(i-1)+g,,t]), prob=cell_prob     [ngene*(i-1)+g,,t], log=T)
        } # g
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
    beta_gamma_save[,k] <- beta_gamma
    beta_phi_save[,k] <- beta_phi
    beta_kappa_save[,k] <- beta_kappa
    alpha_det_save[,k] <- alpha_det
    alpha_cap_save[,k] <- alpha_cap
    if (nmcmc > 500) {
      if (k > nmcmc-500) {
        N_save[,,,k-nmcmc+500] <- N
      }
    }
  } # k

  list(beta0_save=beta0_save, 
       beta_gamma_save=beta_gamma_save, beta_phi_save=beta_phi_save, 
       beta_kappa_save=beta_kappa_save, 
       alpha_det_save=alpha_det_save, alpha_cap_save=alpha_cap_save, 
       N_save=N_save)
} # gene_sipm_mcmc

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
  gene_sipm_mcmc(y=y, r=r, m=m, x=x, w=w, d=d, nmcmc=nmcmc)
} # i
stopImplicitCluster() # clean up the cluster
end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time

output <- list(out=out, time_taken=time_taken)
save(output, file='4.3.4 gene sipm_output.RData')

#==============
# Plot results
#==============
pdf(file='4.3.4.1 gene sipm_chains.pdf', width=10, height=10)

library(rstan) # for calculating rhat statistics

cols <- c('gold', 'tomato', 'royalblue')
ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta["reef"]^"[0]"), expression(beta["temperature"]^"[0]"))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta["density"]^""["["*gamma*"]"]), 
                     expression(beta["reef"]^""["["*gamma*"]"]), expression(beta["temperature"]^""["["*gamma*"]"]), 
                     expression(beta["gene"]^""["["*gamma*"]"]), expression(beta["reef "%*%" gene"]^""["["*gamma*"]"]))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta["density"]^""["["*phi*"]"]), 
                   expression(beta["reef"]^""["["*phi*"]"]), expression(beta["temperature"]^""["["*phi*"]"]), 
                   expression(beta["gene"]^""["["*phi*"]"]), expression(beta["reef "%*%" gene"]^""["["*phi*"]"]))
ylab_beta_kappa <- c(expression(beta[0]^""["["*kappa*"]"]), expression(beta["density"]^""["["*kappa*"]"]), 
                     expression(beta["reef"]^""["["*kappa*"]"]), expression(beta["gene"]^""["["*kappa*"]"]), 
                     expression(beta["reef "%*%" gene"]^""["["*kappa*"]"]))
ylab_alpha_det <- c(expression(alpha[0]^"[det]"), expression(alpha[1]^"[det]"), expression(alpha[2]^"[det]"))
ylab_alpha_cap <- c(expression(alpha[0]^"[cap]"), expression(alpha[1]^"[cap]"), expression(alpha[2]^"[cap]"))

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
  text(x=nmcmc*0.4, y=beta0[i]+yint*8*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+4)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_gamma_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma[i] - yint * 4
  ymax <- beta_gamma[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*2), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_gamma[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_gamma[i]+yint*4*0.6, pos=4,  cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+4)) {
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
  text(x=nmcmc*0.4, y=beta_phi[i]+yint*8*0.6, pos=4, cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(ncovs+3)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$beta_kappa_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_kappa[i] - yint * 8
  ymax <- beta_kappa[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=beta_kappa[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_kappa[i], cex.main=2, line=1.2)
  text(x=nmcmc*0.4, y=beta_kappa[i]+yint*8*0.6, pos=4, cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

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
  if (i == 1) {
    axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  } else {
    axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  }
  axis(2, at=round(seq(ymin, ymax, yint*c(2,1)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha_det[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha_det[i], cex.main=2, line=1.5)
  text(x=nmcmc*0.4, y=alpha_det[i]+yint*c(4,2)[i]*0.6, pos=4, cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

for (i in 1:(npcvs+1)) {
  tt <- matrix(, nmcmc, chain)
  for (j in 1:chain) {
    tt[,j] <- out[[j]]$alpha_cap_save[i,]
  } # j
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc,]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha_cap[i] - yint * c(8,8,16)[i]
  ymax <- alpha_cap[i] + yint * c(8,8,16)[i]
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), cex.axis=1.1)
  axis(2, at=round(seq(ymin, ymax, yint*c(4,4,8)[i]), digits=1), cex.axis=1.2, las=2)
  box()
  for (j in 1:chain) {
    lines(tt[,j], col=cols[j])
  } # j
  abline(h=alpha_cap[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha_cap[i], cex.main=2, line=1.5)
  text(x=nmcmc*0.4, y=alpha_cap[i]+yint*c(8,8,16)[i]*0.6, pos=4, cex=1.6, 
       labels=paste(c('Rhat = ', format(round(Rhat(tt[(nmcmc*0.2+1):nmcmc]), digits=2), nsmall=2)), collapse=''))
} # i

title(xlab='Iteration', cex.lab=3, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=3, line=0.8, outer=T)

dev.off()

pdf(file='4.3.4.2 gene sipm_N.pdf', width=9, height=8)

library(plotrix) # for making pie charts at selected positions

N_post <- array(, dim=c(nsite, nyear, ngene, 1500))
N_post[,,,   1: 500] <- out[[1]]$N
N_post[,,, 501:1000] <- out[[2]]$N
N_post[,,,1001:1500] <- out[[3]]$N
N_med <- apply(N_post, 1:3, median)
years_plot <- c(1, 8)

par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(4,5,5,2))

for (t in 1:2) {
  plot(1, xlim=c(-11, 11), ylim=c(-8.8, 8.8), type='n', axes=F, xlab='', ylab='')
  for (i in 1:nsite) {
    tt <- N[i,years_plot[t],]
    if (sum(tt) > 0) {
      prop <- tt / sum(tt)
      size <- sum(N[i,years_plot[t],]) ^ 0.67 / 14
      cols <- colorRampPalette(c('white','royalblue'))(ngene+1)[-1]
      floating.pie(xpos=cov$lon[i], ypos=cov$lat[i], x=prop, radius=size, col=cols, border=cols[1], lwd=1)
      draw.circle(x=cov$lon[i], y=cov$lat[i], radius=size, nv=100, border='grey18', col=NA, lwd=1.2)
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
    tt <- N_med[i,years_plot[t],]
    if (sum(tt) > 0) {
      prop <- tt / sum(tt)
      size <- sum(N_med[i,years_plot[t],]) ^ 0.67 / 14
      cols <- colorRampPalette(c('white','lightcoral'))(ngene+1)[-1]
      floating.pie(xpos=cov$lon[i], ypos=cov$lat[i], x=prop, radius=size, col=cols, border=cols[1], lwd=1)
      draw.circle(x=cov$lon[i], y=cov$lat[i], radius=size, nv=100, border='grey18', col=NA, lwd=1.2)
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
title(main=expression('Proportion of Gene Class'), line=3.2, cex.main=3, outer=T)

dev.off()


