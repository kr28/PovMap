rm(list=ls())
library(foreign)
#library(lfe)
# library(HI) # implements rejection sampling
library(psych) # for matrix functions
library(MASS)
# setwd("c:/Roy/My WB reports/EU project/!Research")
#setwd("c:/WB/EU/!Research")


########## control values set by user ########
nFV<-3
nFEps<-3
muQuantile<-0.1
nIter<-1300


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
#------------------------------ Read survey data ----------------------
# Missing data was handled in Stata (data is free of missing values for variables included in the model)
#
## data <- read.dta("fakesurvey.dta")
## data$cluster <- data$id_mun
data <- read.dta("data1000.dta",convert.factors=F)
data$cluster <- data$cid
clusterSize<-as.vector(table(data$cluster))
names(data)[4] <- "lny"
#write(clusterSize,file="clusterSize.txt")

# Administrative data are numeric in dataset. Change type into "factors". (To facilitate fixed-effects regression.) "regn" is already a "factor".
## data$prov <- as.factor(data$prov)
#data$id_mun <- as.factor(data$id_mun)
#data$id_bgy <- as.factor(data$id_bgy)

# Run linear regression without FE
#reg.beta <- lm(lny~1+hh_logsize+min1_edu56+hd_female+hd_married+hd_edu4+hd_edu5+hd_edu6+hd_olf+nchild_516+pch_5+dw_basic,data)
## reg.beta <- lm(lny~1+dw_basic+hd_edu4+hd_edu5+hd_edu6+hd_female+hd_married+hd_olf+hh_logsize+min1_edu56+nchild_516+pch_5,data)
reg.beta <- lm(lny~1+x,data=data)
r_hat <- resid(reg.beta)
n_hlds <- length(r_hat)
#write(r_hat,file="r_hat_lm.txt");

# store independent variables in X matrix
#X<-as.matrix(subset(data,select=c(hh_logsize,min1_edu56,hd_female,hd_married,hd_edu4,hd_edu5,hd_edu6,hd_olf,nchild_516,pch_5,dw_basic)))
X<-as.matrix(subset(data,select=c(x)))
X.ols<-matrix(0,dim(X)[1],dim(X)[2]+1)
X.ols[,2:dim(X.ols)[2]]<-X
X.ols[,1]<-matrix(rep(1,dim(X.ols)[1]))

# number of independent variables (variables in X excluding the constant)
d0 <- dim(X)[2]

# and including the constant
d <- d0+1
d.ols <- dim(X.ols)[2]

# --------------------------- cluster level -----------------------
# Get cluster means
ubar <- aggregate(r_hat,by=list(data$cluster),FUN=mean)$x
n_area <- length(ubar)

# H3 variance of epsilon: within estimator
ka <- as.numeric(table(data$cluster))
Dka <- diag(ka)
y <- data$lny
ybar <- aggregate(data$lny,by=list(data$cluster),FUN=mean)$x
a.X <- aggregate(X,by=list(data$cluster),FUN=mean)
Xbar <- as.matrix(a.X[,2:dim(a.X)[2]])
XTX <- (t(X)%*%X)-(t(Xbar)%*%Dka%*%Xbar)
XTy <- (t(X)%*%y)-(t(Xbar)%*%Dka%*%ybar)
beta.within <- ginv(XTX)%*%XTy
eps_hat <- (y-rep(ybar,ka))-(X%*%beta.within-rep(Xbar%*%beta.within,ka))
sigma2.eps <- sum(eps_hat*eps_hat)/(n_hlds-(d-1)-n_area)

v.eps <- sigma2.eps / ka

# H3 variance of v: Corollary 3
a.X.ols <- aggregate(X.ols,by=list(data$cluster),FUN=mean)
UX.ols <- as.matrix(a.X.ols[,2:dim(a.X.ols)[2]])
UX.ols <- Dka%*%UX.ols
trupxu <- tr(ginv(t(X.ols)%*%X.ols)%*%t(UX.ols)%*%UX.ols)
sigma2.v <- (n_hlds-d.ols)*((sum(r_hat*r_hat)/(n_hlds-d.ols))-sigma2.eps)/(n_hlds-trupxu)

#
# obtain GLS estimate of beta
#
xtax <- matrix(0,dim(X.ols)[2],dim(X.ols)[2])
xtay <- matrix(0,dim(X.ols)[2],1)
j <- 1
for(i in 1:n_area) {
   psi <- sigma2.v/(ka[i]*sigma2.v+sigma2.eps)
   #V.inverse <- diag(ka[i])
   V.inverse <- (1/sigma2.eps)*(diag(ka[i])-psi*matrix(1,ka[i],ka[i]))
   xtax <- xtax + t(X.ols[j:(j+ka[i]-1),])%*%V.inverse%*%X.ols[j:(j+ka[i]-1),]
   xtay <- xtay + t(X.ols[j:(j+ka[i]-1),])%*%V.inverse%*%y[j:(j+ka[i]-1)]
   j <- j+ka[i]
}
beta.gls <- ginv(xtax)%*%xtay

#---------------- decompose total error into cluster and household error ---------------------
gamma_a <- sigma2.v/(sigma2.v+(sigma2.eps/ka))
v_tilde <- gamma_a * ubar - mean(gamma_a * ubar);
eps_tilde <- r_hat - rep(v_tilde,ka) - mean(r_hat - rep(v_tilde,ka));

eps_hat.within <- eps_hat;
eps_hat <- (sqrt(sigma2.eps)/sd(eps_tilde))*eps_tilde;
v_hat <- (sqrt(sigma2.v)/sd(v_tilde))*v_tilde;

#--------------------------------- alpha model option -----------------------------------------
A <- 1.05*max(eps_hat*eps_hat)
h <- log(eps_hat*eps_hat/(A-eps_hat*eps_hat))
reg.alpha <- lm(h~1+x,data)
sigma2.alpha <- sum(reg.alpha$residuals*reg.alpha$residuals)/reg.alpha$df.residual
B <- exp(reg.alpha$fitted.values)
sigma2.eps.ah <- A*B/(1+B) + 0.5*sigma2.alpha*(A*B)*(1-B)/((1+B)^3)

# impose restriction that the average conditional variance equals the unconditional variance
sigma2.eps.ah <- sigma2.eps*sigma2.eps.ah/mean(sigma2.eps.ah)

# obtain GLS estimate of beta (with alpha model)
xtax <- matrix(0,dim(X.ols)[2],dim(X.ols)[2])
xtay <- matrix(0,dim(X.ols)[2],1)
j <- 1
for(i in 1:n_area) {
   Deps.inverse <- 1/sigma2.eps.ah[j:(j+ka[i]-1)]
   psi <- sigma2.v/(1+sigma2.v*sum(Deps.inverse))
   V.inverse <- diag(Deps.inverse) - psi*(Deps.inverse%*%t(Deps.inverse))
   xtax <- xtax + t(X.ols[j:(j+ka[i]-1),])%*%V.inverse%*%X.ols[j:(j+ka[i]-1),]
   xtay <- xtay + t(X.ols[j:(j+ka[i]-1),])%*%V.inverse%*%y[j:(j+ka[i]-1)]
   j <- j+ka[i]
}
beta.gls.alpha <- ginv(xtax)%*%xtay

# EB with alpha model
mu_EB <- rep(0,n_area)
var_EB <- rep(0,n_area)
j <- 1
for(i in 1:n_area) {
   mu_EB[i] <- sigma2.v*(sum(r_hat[j:(j+ka[i]-1)]/sigma2.eps.ah[j:(j+ka[i]-1)]))/(1+sigma2.v*sum(1/sigma2.eps.ah[j:(j+ka[i]-1)]))
   var_EB[i] <- sigma2.v*(1-(sigma2.v*sum(1/sigma2.eps.ah[j:(j+ka[i]-1)]))/(1+sigma2.v*sum(1/sigma2.eps.ah[j:(j+ka[i]-1)])))
   j <- j+ka[i]
}


#
# ---------------------------- Fit the mixture-distribution for F_v ----------------------------
set.seed(153)
y <- ubar
n <- length(y) # should match m (i.e. number of clusters)
g <- nFV # we expect to need a larger number of components for the restricted approach
g <- 5

h <- (g-1)/2 # number of components left (right) from center

############### Restricted approach (estimating mixing probabilities only) ################
murange <- quantile(y,prob=c(muQuantile,1-muQuantile))
mu <- rep(0,g)
mu[1:h] <- seq(1,1/h,by=-1/h)*murange[1]
mu[h+1] <- 0
mu[(h+2):g] <- seq(1/h,1,by=1/h)*murange[2]
sig <- rep(sigma2.v/g,g)^0.5
sig[h+1] <- sigma2.v^0.5 # variance of center component is set to the variance of the data
# sig[h+1] <- sig 
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)
p0<-p      # keep it for later use in F-ebar

# check the initial likelyhood 
fv.blok.restricted(1)

# Apply EM-algorithm
psi.v.restricted <- fv.blok.restricted(nIter)

# Keep mixture distribution for Fv using parameter estimates
g.v.restricted<-g
p.v.restricted=psi.v.restricted$p
mu.v.restricted=psi.v.restricted$mu
sig.v.restricted=psi.v.restricted$sig

fv.dens.restricted <- function(x) as.numeric(sapply(x, function(t) p.v.restricted %*% dnorm(t,m=mu.v.restricted,s=sig.v.restricted)))

hist(y,prob=T,n=50,main="histogram of cluster mean residuals and eta densities")
plot(fv.dens.restricted,xlim=c(-1.0,1.0),add=T)
plot(function(x) fv.nm(x,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),
                            s.v=c(0.2,0.2,0.2)),-1,1,add=T,col="red")

#
# ---------------------------------- Fit the mixture-distribution for F_eps ----------------------------
#
# using the household data
set.seed(153)
g <- 3 #
g <- 5 
# select the first household from each cluster
y <- rep(0,n_area)
j <- 10
for(i in 1:n_area) {
   y[i] <- r_hat[j]
   j <- j+ka[i]
}
n <- length(y)      # should match total number of clusters
#g <- 5
h <- (g-1)/2 # number of components left (right) from center
y <- r_hat
############### Restricted approach (estimating mixing probabilities only) ################
murange <- quantile(y,prob=c(muQuantile,1-muQuantile))
mu <- rep(0,g)
mu[1:h] <- seq(1,1/h,by=-1/h)*murange[1]
mu[h+1] <- 0
mu[(h+2):g] <- seq(1/h,1,by=1/h)*murange[2]
sig <- rep(sigma2.eps/g,g)^0.5
sig[h+1] <- sigma2.eps^0.5 # variance of center component is set to the variance of the data
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)

# check the initial likelyhood 
feps.blok.restricted(1)

# Apply EM-algorithm
psi.eps.restricted <- feps.blok.restricted(nIter)

# construct mixture distribution for F_eps using parameter estimates
g.eps.restricted<-g
p.eps.restricted=psi.eps.restricted$p
mu.eps.restricted=psi.eps.restricted$mu
sig.eps.restricted=psi.eps.restricted$sig

feps.dens.restricted <- function(x) as.numeric(sapply(x, function(t) p.eps.restricted %*% dnorm(t,m=mu.eps.restricted,s=sig.eps.restricted)))

plot(feps.dens.restricted,xlim=c(-1.7,1.7),col="red",add=T)

hist(y,prob=T,n=50,main="histogram of residuals and epsilon densities",ylim=c(0,0.5))
plot(feps.dens.restricted,xlim=c(-3.0,3.0),add=T)
plot(function(x) fv.nm(x,p=c(0.2,0.5,0.3),m.v=c(-1.55,0,0.77),
                            s.v=sqrt(c(0.45,0.45,0.45))),-3,3,add=T,col="red")

plot(0,0,main="eta densities",xlab="area error",ylab="density",xlim=c(-1,1),ylim=c(0,3),type="n")
plot(function(x) fv.nm(x,p=psi.v.restricted$p,m.v=psi.v.restricted$mu,
                            s.v=psi.v.restricted$sig),-1,1,add=T,col="red")
plot(function(x) fv.cond.nm(x,mea=-0.3,p=psi.v.restricted$p,m.v=psi.v.restricted$mu,
                            s.v=psi.v.restricted$sig,
                            s2.eps=sigma2.eps,k=30),-1,1,add=T)
plot(function(x) fv.cond.nm(x,mea=0.0,p=psi.v.restricted$p,m.v=psi.v.restricted$mu,
                            s.v=psi.v.restricted$sig,
                            s2.eps=sigma2.eps,k=30),-1,1,col="blue",add=T)
plot(function(x) fv.cond.nm(x,mea=0.3,p=psi.v.restricted$p,m.v=psi.v.restricted$mu,
                            s.v=psi.v.restricted$sig,
                            s2.eps=sigma2.eps,k=30),-1,1,col="green",add=T)

plot(0,0,main="eta densities",xlab="area error",ylab="density",xlim=c(-1,1),ylim=c(0,1.5),type="n")
plot(function(x) fv.nm(x,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),
                            s.v=c(0.2,0.2,0.2)),-1,1,add=T)
plot(function(x) fv.n(x,m=0,s=sqrt(0.1)),-1,1,add=T,col="red")
lines(c(-0.5,-0.5),c(0,1.5),col="blue",add=T)

#########################
class(rmultinom(1,1,c(0.3,0.4,0.4))[,1])
# sample from v-mixture

draw.one <- function(pvec,muvec,sigvec){ # vectors should be conformal
  choice <- (rmultinom(1,size=1,prob=pvec)[,1])
  mu <- muvec %*% choice
  sig <- sigvec %*% choice
  rnorm(1,mu,sig)
}

vdraw <-  sapply(1:length(ubar),function(i) draw.one(p.v.restricted,mu.v.restricted,sig.v.restricted))

epsdraw <- sapply(1:length(data$lny),function(i) draw.one(p.eps.restricted,mu.eps.restricted,sig.eps.restricted))

#########################
### ------------------ sample from the estimated unconditional distribution f(v) ---------------------

# fv.dens.restricted (unconditional density for v)

# use "rejection sampling" to sample from the unconditional distribution
my.domain.fun <- function(x) x>-1 & x < 1
my.ldens <- function(x) log(fv.dens.restricted(x))

# sample 5000 observations
v <- arms(runif(1,-1,1), my.ldens, my.domain.fun, 5000)


### ------------------ now derive and sample from the conditional distribution f(v|ebar) obtained using Bayes rule ---------------------

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
ebar <- 0.5
fv.con.dens1 <- function(x) fv.con.dens(x,ebar)

# sample 5000 observations
cv <- arms(runif(1,-1,1), function(x,u) my.con.ldens(x,u), function(x,u) my.domain.fun(x), 5000, u=ebar)


### ------------------ now derive f(v|ebar) as a normal-mixture using Proposition 17 ---------------------

# need to set benchmark k for ebar_a
k <- 20

# get s2_a = var[eps]/k
s2 <- sigma2.eps/k

# set value of ebar for test/example
ebar <- 0
#fv.con.dens1 <- function(x) fv.con.dens(x,ebar)

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

fv.cond.nm <- function(et,mea,p.v,m.v,s.v,s2.eps,k){ # mea: mean residual in area a
  s2 <- s2.eps / k
  alpha <-s.v^2/(s2+s.v^2)
  s.cv <- (1/s2+1/s.v^2)^(-1/2)
  m.cv <- alpha*mea + (1-alpha)*m.v
  delta <- alpha*mea^2 + (1-alpha)*m.v^2-m.cv^2
  cstar <- p.v*exp(-delta/(2*s.cv^2))/(s2+s.v^2)^(1/2)
  p.cv<-cstar/sum(cstar)
  as.numeric(sapply(et, function(t) p.cv %*% dnorm(t,m=m.cv,s=s.cv)))
}

fv.cond.n <- function(et,mea,s.v,s2.eps,k){ # mea: mean residual in area a
  s2 <- s2.eps / k
  gamma <-s.v^2/(s2+s.v^2)
  s.cv <- (1/s2+1/s.v^2)^(-1/2)
  m.cv <- gamma * mea 
  as.numeric(sapply(et, function(t) dnorm(t,m=m.cv,s=s.cv)))
}

fv.nm <- function(et,p.v,m.v,s.v) {
  as.numeric(sapply(et, function(t) p.v %*% dnorm(t,m=m.v,s=s.v)))
}

fv.n <- function(et,m,s) dnorm(et,m,s)

parms.cond.nm <- function(mea,p.v,m.v,s.v,s2.eps,k){ # mea: mean residual in area a
  s2 <- s2.eps / k
  alpha <-s.v^2/(s2+s.v^2)
  s.cv <- (1/s2+1/s.v^2)^(-1/2)
  m.cv <- alpha*mea + (1-alpha)*m.v
  delta <- alpha*mea^2 + (1-alpha)*m.v^2-m.cv^2
  cstar <- p.v*exp(-delta/(2*s.cv^2))/(s2+s.v^2)^(1/2)
  p.cv<-cstar/sum(cstar)
  ## as.numeric(sapply(et, function(t) p.cv %*% dnorm(t,m=m.cv,s=s.cv)))
  list(p=p.cv,m=m.cv,s=s.cv)
}

parms.cond.n <- function(mea,s.v,s2.eps,k){ # mea: mean residual in area a
  s2 <- s2.eps / k
  gamma <-s.v^2/(s2+s.v^2)
  s.cv <- (1/s2+1/s.v^2)^(-1/2)
  m.cv <- gamma * mea 
  ## as.numeric(sapply(et, function(t) dnorm(t,m=m.cv,s=s.cv)))
  list(m=m.cv,s=s.cv)
}

parms.cond.nm(x,mea=0.0,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),s.v=c(0.2,0.2,0.2),s2.eps=0.9,k=30)

quartz()
plot.new()
plot(0,0,main="conditional densities",xlab="area error",ylab="density",xlim=c(-1,1),ylim=c(0,3),type="n")
plot(function(x) fv.cond.nm(x,mea=0.0,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),
                            s.v=c(0.2,0.2,0.2),s2.eps=0.9,k=30),-1,1,add=T)
plot(function(x) fv.nm(x,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),
                            s.v=c(0.2,0.2,0.2)),-1,1,add=T,col="red")
plot(function(x) fv.cond.n(x,mea=0.2,s.v=sqrt(0.1),s2.eps=0.9,k=30),-0.4,0.7,add=T,col="blue")

integrate(function(x) x*fv.cond.nm(x,mea=0.2,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),s.v=c(0.2,0.2,0.2),s2.eps=0.9,k=30),-0.4,0.7)
integrate(function(x) x*fv.cond.n(x,mea=0.2,s.v=sqrt(0.1),s2.eps=0.9,k=30),-0.4,0.7)

integrate(function(x) fv.nm(x,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),s.v=c(0.2,0.2,0.2)),-1,1)

plot(0,0,main="eta densities",xlab="area error",ylab="density",xlim=c(-1,1),ylim=c(0,3),type="n")
plot(function(x) fv.nm(x,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),
                            s.v=c(0.2,0.2,0.2)),-1,1,add=T,col="red")
plot(function(x) fv.nm(x,p=p.v.restricted,m.v=mu.v.restricted,
                            s.v=sig.v.restricted),-1,1,add=T)


plot(0,0,main="eps densities",xlab="individual error",ylab="density",xlim=c(-3,3),ylim=c(0,0.6),type="n")
plot(function(x) fv.nm(x,p=c(0.2,0.5,0.3),m.v=c(-1.55,0,0.77),s.v=sqrt(c(0.45,0.45,0.45))),-3,3,add=T)
plot(function(x) dnorm(x,m=0,s=sqrt(0.9)),-3,3,add=T,col="red")

## example: area (single eta) with 5000 observations and survey mean residual mea

# 1. unconditional normal eps and eta
povline <- 4.58 
set.seed(101)
x.1 <- rnorm(5000,m=2)
eps.1 <- rnorm(5000,m=0,s=sqrt(0.9))

pov.1 <- function(et) {
  sum(x.1+eps.1+3+et < povline) / length(x.1)
}
pov.1(0) # 0.3848

# 2. unconditional mixed normal eps and eta, same x
peps <- c(0.2,0.5,0.3)
meps <- c(-1.55,0,0.77)
seps <- sqrt(c(0.45,0.45,0.45))
eps.2 <-  sapply(1:5000, function(i) draw.one(peps,meps,seps))
pov.2 <- function(et) {
  sum(x.1+eps.2+3+et < povline) / length(x.1)
}
pov.2(0) # 0.3966

# 3. conditional normal eta, same x, same eps
# conditional expectation of eta, given mea = 0

integrate(function(x) x*fv.cond.n(x,mea=0,s.v=sqrt(0.1),s2.eps=0.9,k=30),-1,1) # 0
pov.1(0) # 0.3848

#4. conditional mixed normal eta, same x, same eps
et.4 <- integrate(function(x) x*fv.cond.nm(x,mea=0,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),s.v=c(0.2,0.2,0.2),s2.eps=0.9,k=30),-1,1) $value # -0.02349776
pov.2(et.4) # 0.4028

# 1. simulation unconditional normal eta
eta.1 <- rnorm(1000,m=0,s=sqrt(0.1))
sim.1 <- sapply(eta.1,pov.1)
# hist(sim.1,nclass=50)
s1 <- c(mean=mean(sim.1),sd=sd(sim.1),hlow=quantile(sim.1,0.025),hhigh=quantile(sim.1,0.975))
  ##      mean          sd   hlow.2.5% hhigh.97.5% 
  ## 0.3865520   0.0796706   0.2328000   0.5524250 

#2. simulation unconditional mixed normal eta and eps
eta.2 <- sapply(1:1000,function (i) draw.one(c(0.4,0.4,0.2),c(-0.22,0,0.44),c(0.2,0.2,0.2)))
sim.2 <- sapply(eta.2,pov.2)
#hist(sim.2,nclass=50)
s2 <- c(mean=mean(sim.2),sd=sd(sim.2),hlow=quantile(sim.2,0.025),hhigh=quantile(sim.2,0.975))
 ##       mean          sd   hlow.2.5% hhigh.97.5% 
 ## 0.40056560  0.07454577  0.25338500  0.53988500 

# 3. simulation conditional normal eta, same x, same eps
# conditional distribution of eta, given mea = 0

sol.3 <- parms.cond.n(0,sqrt(0.1),0.9,30)
eta.3 <- rnorm(1000,sol.3$m,sol.3$s)
sim.3 <- sapply(eta.3,pov.1)
#hist(sim.3,nclass=50)
s3 <- c(mean=mean(sim.3),sd=sd(sim.3),hlow=quantile(sim.3,0.025),hhigh=quantile(sim.3,0.4))
  ##      mean          sd   hlow.2.5% hhigh.97.5% 
  ## 0.3880064   0.0390752   0.3113600   0.4654300 

#4. simulation conditional mixed normal eta, same x, same eps
sol.4 <- parms.cond.nm(0,p=c(0.4,0.4,0.2),m.v=c(-0.22,0,0.44),s.v=c(0.2,0.2,0.2),s2.eps=0.9,k=30)
eta.4 <- sapply(1:1000, function(i) draw.one(sol.4$p,sol.4$m,sol.4$s))
sim.4 <- sapply(eta.4,pov.2)
#hist(sim.4,nclass=50)
s4 <- c(mean=mean(sim.4),sd=sd(sim.4),hlow=quantile(sim.4,0.025),hhigh=quantile(sim.4,0.975))
 ##       mean          sd   hlow.2.5% hhigh.97.5% 
 ## 0.40473380  0.03865883  0.32959000  0.47840000 

s1.4 <- data.frame(matrix(c(s1,s2,s3,s4),ncol=4,byrow=T))
names(s1.4) <- c("mean","sd","q_0.025","q_0.975")
s1.4$simulation <- c("normal-unconditional", "mixed.nrm-unc","normal-cond","mixed.nrm-cond")
s1.4
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
#hist(cv,prob=T,n=50)
plot(fv.con.dens.prop17,xlim=c(-0.7,0.7),col="blue",add=T)
plot(fv.con.dens1,xlim=c(-0.7,0.7),col="blue",add=T)

### PERFECT MATCH! NOTE HOW CONDITIONAL DENSITY DERIVED AS A NORMAL-MIXTURE NICELY MATCHES THE CONDITIONAL DENSITY OBTAINED USING EQ. (48)
### THE 5000 OBSERVATIONS SAMPLED FROM THE NORMAL-MIXTURE ALSO PERFECTLY MATCH THE ANALYTIC MIXTURE DENSITY
### NOTE THAT THE SMALL DISREPANCY BETWEEN THE CONDIOTONAL DENSITY OBTAINED USING PROP. 17 AND OBTAINED USING EQ. (48)
### IS BECAUSE THE LATTER DENSITY IS SUBJECT TO ESTIMATION ERROR IN f(ebar) WHILE THE FORMER IS NOT
