########## Aug 18 Freezed MixFit 4 built on:
########## Aug 10 Freezed MixFit3 Modified by Qinghua ############

rm(list=ls())
library(foreign)
library(lfe)
library(HI) # implements rejection sampling
setwd("c:/WB/EU/!Research")


########## control values set by user ########
nFV<-7
nFEps<-5
muQuantile<-0.1
nIter<-5000


### vector y,p and scalor n and g should be adjusted when calling the following density function
fv.blok.restricted <- function(niter) {
  sig.sq <- sig^2
  dnormMatrix <- outer(1:n,1:g,FUN=function(j,i) dnorm(y[j],m=mu[i],s=sqrt(sig.sq[i]+v.eps[j])))
  for(k in 1:niter) {
    tau<-dnormMatrix %*% diag(p)	# E-step
    llhm <- log(tau)
    tau<-tau/rowSums(tau)
    llh <- sum(tau * llhm)
    p<-colSums(tau)/n;        		# M-step
  }
  list(p=p,mu=mu,sig=sig,llh=llh) 
}

feps.blok.restricted <- function(niter) {
  gs.fv <- Vectorize(function(x,mui,sig2i) sum(dnorm(x,m=mui+mu.v.restricted,s=sqrt(sig2i+sig.v.restricted^2))*p.v.restricted))
  sig.sq <- sig^2
  dnormMatrix <- outer(1:n,1:g,FUN=function(j,i) gs.fv(y[j],mu[i],sig.sq[i]))
  for(k in 1:niter) {
    tau<-dnormMatrix %*% diag(p)	# E-step
    llhm<-log(tau)
    tau<-tau/rowSums(tau)
    llh <- sum(tau * llhm)
    p<-colSums(tau)/n;          	# M-step
  }
  list(p=p,mu=mu,sig=sig,llh=llh)
}


######################## main program ##########################
#
#------------------------------ Read Philippine data 2009 ----------------------
# Missing data was handled in Stata (data is free of missing values for variables included in the model)
#
data <- read.dta("phi2009.dta")
clusterSize<-as.vector(table(data$id_bgy))
#write(clusterSize,file="clusterSize.txt")

# Administrative data are numeric in dataset. Change type into "factors". (To facilitate fixed-effects regression.) "regn" is already a "factor".
data$prov <- as.factor(data$prov)
data$id_mun <- as.factor(data$id_mun)
data$id_bgy <- as.factor(data$id_bgy)

# Run linear regression without FE
lmreg <- lm(lny~1+urban+hh_logsize+hd_female+hd_logage+hd_married+hd_unemp+hd_olf+hd_edu4+hd_edu5+hd_edu6+min1_emp+min1_edu56+pch_5+pch_516+peld_60+dw_basic+urb_logsize+urb_hd_fem+urb_hd_unemp+urb_hd_olf+urb_hd_edu5+urb_hd_edu6+urb_min_emp+urb_min_edu56+urb_pch5+urb_pch516+urb_peld60+urb_dw_basic,data)
eps <- resid(lmreg)
#write(eps,file="eps_lm.txt");

# --------------------------- Barangay level -----------------------
# Get bgy means
ubar <- aggregate(eps,by=list(data$id_bgy),FUN=mean)$x # don't ask; check documentation of aggregate() function

# With bgy FE (to get H3-estimate of var[eps])
fereg_bgy <- felm(lny~1+hh_logsize+hd_female+hd_logage+hd_married+hd_unemp+hd_olf+hd_edu4+hd_edu5+hd_edu6+min1_emp+min1_edu56+pch_5+pch_516+peld_60+dw_basic+urb_logsize+urb_hd_fem+urb_hd_unemp+urb_hd_olf+urb_hd_edu5+urb_hd_edu6+urb_min_emp+urb_min_edu56+urb_pch5+urb_pch516+urb_peld60+urb_dw_basic+G(id_bgy),data=data)
eps_felm <- summary(fereg_bgy)$residuals
#write(eps_felm,file="eps_felm.txt");

v.eps.bgy <- summary(fereg_bgy)$rse^2    # v.eps.bgy<-sqrt((t(eps_felm) %*% eps_felm)/(38369-26-3015))

# Get var(eps)/k_a
v.eps <- v.eps.bgy / as.numeric(table(data$id_bgy))
#
# ---------------------------- Fit the mixture-distribution for F_v ----------------------------
set.seed(153)
y <- ubar
n <- length(y) # should match m (i.e. number of clusters)
g <- nFV # we expect to need a larger number of components for the restricted approach

############### Restricted approach (estimating mixing probabilities only) ################
murange <- quantile(y,prob=c(muQuantile,1-muQuantile))
mu <- seq(murange[1],murange[2],by=(murange[2]-murange[1])/(g-1))
sig<-rep(var(y)/g,g)^0.5
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)
p0<-p      # keep it for later use in F-ebar

# check the initial likelyhood 
fv.blok.restricted(1)

# Apply EM-algorithm (note that 5000 iterations needed for convergence)
psi.v.restricted <- fv.blok.restricted(nIter)

# Keep mixture distribution for Fv using parameter estimates
g.v.restricted<-g
p.v.restricted=psi.v.restricted$p
mu.v.restricted=psi.v.restricted$mu
sig.v.restricted=psi.v.restricted$sig

fv.dens.restricted <- function(x) as.numeric(sapply(x, function(t) p.v.restricted %*% dnorm(t,m=mu.v.restricted,s=sig.v.restricted)))

#
# -------------------------- Fit the mixture-distribution for F_ebar. --------------------------------
#  This trunk can be ignored if not do directly estimate f(Va|e_bar_a)=f(e_bar_a|Va)*f(Va)/f(e_bar_a)
#  Done by set v.eps to zero and recompute fv_blok_restricted
#  
v.eps<-0*v.eps  # we will not estimate the distribution of ebar = v + epsbar, which is data that we observe, so we set the variance of the "measurement error" to zero
p<-p0

# check the initial likelyhood 
fv.blok.restricted(1)

# Apply EM-algorithm (note that 5000 iterations needed for convergence)
psi.ebar.restricted <- fv.blok.restricted(5000);
psi.ebar.restricted             # psi contains final parameter estimates

# construct mixture distribution for Fv using parameter estimates
g.ebar.restricted<-g
p.ebar.restricted=psi.ebar.restricted$p
mu.ebar.restricted=psi.ebar.restricted$mu
sig.ebar.restricted=psi.ebar.restricted$sig

febar.dens.restricted <- function(x) as.numeric(sapply(x, function(t) p.ebar.restricted %*% dnorm(t,m=mu.ebar.restricted,s=sig.ebar.restricted)))

#
# ---------------------------------- Fit the mixture-distribution for F_eps ----------------------------
#
# using the household data
set.seed(153)
y <- eps
n <- length(y)      # should match total number of households
g <- nFEps          # we expect to need a larger number of components for the restricted approach

############### Restricted approach (estimating mixing probabilities only) ################
murange <- quantile(y,prob=c(muQuantile,1-muQuantile))
mu <- seq(murange[1],murange[2],by=(murange[2]-murange[1])/(g-1))
sig<-rep(var(y)/g,g)^0.5
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)

# check the initial likelyhood 
feps.blok.restricted(1)

# Apply EM-algorithm (note that 1000 iterations needed for convergence)
psi.eps.restricted <- feps.blok.restricted(300);psi.eps.restricted # psi contains final parameter estimates

# construct mixture distribution for F_eps using parameter estimates
g.eps.restricted<-g
p.eps.restricted=psi.ebar.restricted$p
mu.eps.restricted=psi.ebar.restricted$mu
sig.eps.restricted=psi.ebar.restricted$sig

feps.dens.restricted <- function(x) as.numeric(sapply(x, function(t) p.eps.restricted %*% dnorm(t,m=mu.eps.restricted,s=sig.eps.restricted)))


### ------------------ sample from the estimated unconditional distribution f(v) ---------------------

# fv.dens.restricted (unconditional density for v)

# use "rejection sampling" to sample from the unconditional distribution
my.domain.fun <- function(x) x>-1 & x < 1
my.ldens <- function(x) log(fv.dens.restricted(x))

# sample 5000 observations
v <- arms(runif(1,-1,1), my.ldens, my.domain.fun, 5000)


### ------------------ now derive and sample from the conditional distribution f(v|ebar) ---------------------

# need to set benchmark k for ebar_a
k <- 10

# febar.dens.restricted (unconditional density for ebar)
# fv.dens.restricted (unconditional density for v)
# condition on u=ebar

# derive conditional density function: f (v | ebar)
# note that: f(v_a | ebar_a) = f(ebar | v) f(v) / f(ebar) ~ f(ebar | v) f(v) -- eq. (46) in saenote [f(ebar) ends up in the integration constant]
# note that: f(ebar | v) is a normal pdf with mean=v and variance=var(eps)/k
fv.con.dens <- function(x,u) dnorm(u,m=x,s=(v.eps.bgy/k)^0.5)*fv.dens.restricted(x)/febar.dens.restricted(u)

# use "rejection sampling" to sample from the conditional distribution
my.domain.fun <- function(x) x>-1 & x < 1
my.con.ldens <- function(x,u) log(dnorm(u,m=x,s=(v.eps.bgy/k)^0.5))+log(fv.dens.restricted(x))

# set value of ubar on which we will be conditioning
ebar<-0.09
fv.con.dens1 <- function(x) fv.con.dens(x,ebar)

# sample 5000 observations
cv <- arms(runif(1,-1,1), function(x,u) my.con.ldens(x,u), function(x,u) my.domain.fun(x), 5000, u=ebar)


### ------------------ now derive f(v|ebar) as a normal-mixture (see Proposition 17) ---------------------

# need to set benchmark k for ebar_a
k <- 10

# get s2_a = var[eps]/k
s2<-v.eps.bgy/k

# set value of ebar for test/example
ebar<-0.09
fv.con.dens1 <- function(x) fv.con.dens(x,ebar)

# compute mixture parameters for f(v|ebar) using Proposition 17
alpha<-sig.v.restricted^2/(s2+sig.v.restricted^2)
sig.cv.restricted<-(1/s2+1/(sig.v.restricted^2))^(-1/2)
mu.cv.restricted<-alpha*ebar+(1-alpha)*mu.v.restricted
delta<-alpha*ebar^2+(1-alpha)*mu.v.restricted^2-mu.cv.restricted^2
cstar<-p.v.restricted*exp(-delta/(2*sig.cv.restricted^2))/((s2+sig.v.restricted^2)^(1/2))
p.cv.restricted<-cstar/sum(cstar)

fv.con.dens.prop17 <- function(x) {  
  as.numeric(sapply(x, function(t) p.cv.restricted %*% dnorm(t,m=mu.cv.restricted,s=sig.cv.restricted)))
}

# sample 5000 observations from the normal-mixture f(v|ebar)
n.sample<-5000
rm.comp<-rmultinom(n.sample,size=1,prob=p.cv.restricted)

# for some reason rnorm(n.sample,m=mu.cv.restricted,s=sig.cv.restricted) does not work
rn.comp<-matrix(0,nrow=g.v.restricted,ncol=n.sample)
for(j in 1:g.v.restricted) {
  rn.comp[j,]<-rnorm(n.sample,m=mu.cv.restricted[j],s=sig.cv.restricted[j])
}
cv.sample<-colSums(rn.comp*rm.comp)

# compare histogram to analytic density functions
hist(cv.sample,prob=T,n=50)
plot(fv.con.dens.prop17,xlim=c(-0.7,0.7),col="red",add=T)
plot(fv.con.dens1,xlim=c(-0.7,0.7),col="blue",add=T)

### PERFECT MATCH! NOTE HOW CONDITIONAL DENSITY DERIVED AS A NORMAL-MIXTURE NICELY MATCHES THE CONDITIONAL DENSITY OBTAINED USING EQ. (48)
### THE 5000 OBSERVATIONS SAMPLED FROM THE NORMAL-MIXTURE ALSO PERFECTLY MATCH THE ANALYTIC MIXTURE DENSITY
### NOTE THAT THE SMALL DISREPANCY BETWEEN THE CONDIOTONAL DENSITY OBTAINED USING PROP. 17 AND OBTAINED USING EQ. (48)
### IS BECAUSE THE LATTER DENSITY IS SUBJECT TO ESTIMATION ERROR IN f(ebar) WHILE THE FORMER IS NOT