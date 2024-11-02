#========================================================================
# Code for generating and visualizing habitat patch size and temperature
# written by Qing Zhao, 2023 in Colorado
#========================================================================

#setwd('')

#=====================
# Simulate covariates
#=====================
library(boot) # for using inv.logit functions

set.seed(1)

nsite <- 360 # number of sites
nyear <- 20  # number of years

lon <- runif(nsite, -10, 10) # longitude or easting
lat <- runif(nsite, - 8,  8) # latitude or northing
d <- matrix(, nsite, nsite) # distance between sites
for (i1 in 1:nsite) {
  for (i2 in 1:nsite) {
    d[i1,i2] <- sqrt((lon[i1] - lon[i2]) ^ 2 + (lat[i1] - lat[i2]) ^ 2)
  } # i2
} # i1
diag(d) <- NA

size <- matrix(, nsite, nyear) # habitat patch size
size_mean <- numeric(nsite)
for (i in 1:nsite) {
  size_mean[i] <- exp(rnorm(n=1, mean=log(min(d[i,],na.rm=T) / 2), sd=0.05))
  for (t in 1:nyear) {
    size[i,t] <- exp(rnorm(n=1, mean=log(size_mean[i]), sd=0.05))
  } # t
} # i

temp <- matrix(, nsite, nyear) # temperature
temp_mean <- inv.logit(lat / 30 + lon / 20) * 10
for (t in 1:nyear) {
  temp[,t] <- (temp_mean + sin(t / 1.5) * 2) / 0.6 - 4
} # t

# Save simulated covariate data
cov <- list(lon=lon, lat=lat, size=size, temp=temp)
save(cov, file='covariates.RData')

#==========
# Graphing
#==========
pdf(file='1.4 covariates.pdf', width=12, height=10)

ngrid <- 10000
x <- rep(seq(-11,11,length.out=sqrt(ngrid)), times=sqrt(ngrid))
y <- rep(seq(- 9, 9,length.out=sqrt(ngrid)), each =sqrt(ngrid))
temp_area_mean <- numeric(ngrid)
temp_area <- matrix(, ngrid, nyear)
for (i in 1:ngrid) {
  temp_area_mean[i] <- inv.logit(x[i] / 30 + y[i] / 20) * 10
  for (t in 1:nyear) {
    temp_area[i,t] <- (temp_area_mean[i] + sin(t / 1.5) * 2) / 0.6 - 4
  } # t
} # i

nscale <- 100
scale <- round((temp_area - min(temp_area)) / (max(temp_area) - min(temp_area)) * (nscale - 1)) + 1
xsize <- (max(x) - min(x)) / (sqrt(ngrid)-1) / 2
ysize <- (max(y) - min(y)) / (sqrt(ngrid)-1) / 2

laymat <- matrix(
  c(rep(rep( 1: 2, each=5), times=2), 
    rep(rep( 3: 7, each=2), times=3), 
    rep(rep( 8:12, each=2), times=3), 
    rep(rep(13:17, each=2), times=3), 
    rep(rep(18:22, each=2), times=3)), 
  nrow=14, ncol=10, byrow=TRUE)
layout(mat=laymat)
par(mar=c(1,0,0,0))
par(oma=c(6,7,0,1))

plot(1, type='n', axes=F, xlab='', ylab='')
points(x=seq(0.8,1.2,length.out=5), y=rep(0.92,5), pch=21, bg='white', col='navy', lwd=2, cex=seq(1,9,2)/2)
text(x=0.8, y=0.86, labels='small',      pos=1, cex=2.4)
text(x=1.2, y=0.86, labels='large',      pos=1, cex=2.4)
text(x=1.0, y=1.05, labels='Habitat patch', pos=3, cex=3.0)

plot(1, type='n', axes=F, xlab='', ylab='')
nscale <- 100
for (i in 1:nscale) {
  rect(seq(0.8,1.2,length.out=nscale+1)[i],   0.88, 
       seq(0.8,1.2,length.out=nscale+1)[i+1], 0.96, 
       col=rev(heat.colors(nscale))[i], border=NA)
} # i
text(x=0.8, y=0.86, labels='low',         pos=1, cex=2.4)
text(x=1.2, y=0.86, labels='high',        pos=1, cex=2.4)
text(x=1.0, y=1.05, labels='Temperature', pos=3, cex=3.0)

par(mar=c(0,0,0,0))
for (t in 1:nyear) {
  plot(1, xlim=c(-11,11), ylim=c(-9,12), type='n', axes=F, xlab='', ylab='')
  for (i in 1:ngrid) {
    rect(x[i]-xsize,y[i]-ysize,x[i]+xsize,y[i]+ysize,col=rev(heat.colors(nscale))[scale[i,t]],border=NA)
  } # i
  points(lon, lat, pch=21, cex=size[,t]*3, col='navy', lwd=0.5)

  if (t %% 5 == 1) {
    axis(2, at=seq(-08,08,length.out=5), cex.axis=1.8, las=2)
  }
  if (t > 15) {
    axis(1, at=seq(-10,10,length.out=5), cex.axis=1.8)
  }
  box()
  text(x=0, y=11, labels=paste('Year', t, sep=' '), cex=2.4)
  if (t == 11) {
    axis(2, at=12, labels='Northing', line=3.2, cex.axis=4, tick=F)
  }
  if (t == 18) {
    axis(1, at= 0, labels='Easting',  line=3.4, cex.axis=4, tick=F)
  }
} # t

dev.off()


