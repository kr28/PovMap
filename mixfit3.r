# Fitting mixture-distribution to Philippine data (2003, 2006 and 2009)

library(HI) # implements rejection sampling

fv.blok.restricted <- function(niter=1000,inner.iter=10,printen=F) {
  #
  # arguments
  #   niter is the number of iterations to perform; more advanced is to stop upon convergence
  #   printen is a flag; if set to T it prints results for each iteration
  #
  # The many 'as.numeric' wrappers in the code are to get rid of matrix or list structures
  # I don't know if they are all really necessary but this appears to work.
  #

    
  for(k in 1:niter) {

                                        # E-step
    sig.sq <- sig^2
    tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*dnorm(y[j],m=mu[i],s=sqrt(sig.sq[i]+v.eps[j])))
    llhm <- log(tau)
    tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
    tau <- tau/tau.rs
    # warning:
    # this works because of the way R stores arrays and the dangerous 'recycling rule'
    # safer would be to define tau.rs by something like
    # tau.rs <- matrix(rep(tau %*% rep(1,g),each=g),ncol=g,byrow=T)

    llh <- sum(tau * llhm)
    progr.llh <<- rbind(progr.llh,c(p,mu,sig,llh)) # note the double arrow: this changes global variable progr.llh
    

                                        # M-step

    s.vec.n <- rep(1,n)
    
    # for p there is a closed-form solution
    p <- as.numeric(s.vec.n %*% tau) / n

    # note that mu and sigma are kept fixed in this case, only p is being estimated/updated
  }
  c(p=p,mu=mu,sig=sig) # value of blok function (put all parameters in psi vector)
}

feps.blok.restricted <- function(niter=1000,inner.iter=10,printen=F) {
  #
  # arguments
  #   niter is the number of iterations to perform; more advanced is to stop upon convergence
  #   printen is a flag; if set to T it prints results for each iteration
  #
  # The many 'as.numeric' wrappers in the code are to get rid of matrix or list structures
  # I don't know if they are all really necessary but this appears to work.
  #

  gs.fv <- Vectorize(function(x,m.gs,s.gs) sum(dnorm(x,m=m.gs+mu.v.restricted,s=sqrt(s.gs^2+sig.v.restricted^2))*p.v.restricted))
    
  for(k in 1:niter) {

                                        # E-step
    sig.sq <- sig^2
    tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*gs.fv(y[j],mu[i],sig[i]))
    llhm <- log(tau)
    tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
    tau <- tau/tau.rs
    # warning:
    # this works because of the way R stores arrays and the dangerous 'recycling rule'
    # safer would be to define tau.rs by something like
    # tau.rs <- matrix(rep(tau %*% rep(1,g),each=g),ncol=g,byrow=T)

    llh <- sum(tau * llhm)
    progr.llh <<- rbind(progr.llh,c(p,mu,sig,llh)) # note the double arrow: this changes global variable progr.llh
    

                                        # M-step

    s.vec.n <- rep(1,n)
    
    # for p there is a closed-form solution
    p <- as.numeric(s.vec.n %*% tau) / n

    # note that mu and sigma are kept fixed in this case, only p is being estimated/updated
  }
  c(p=p,mu=mu,sig=sig) # value of blok function (put all parameters in psi vector)
}

fv.blok.general <- function(niter=1000,inner.iter=10,printen=F) {
  #
  # arguments
  #   niter is the number of iterations to perform; more advanced is to stop upon convergence
  #   printen is a flag; if set to T it prints results for each iteration
  #
  # The many 'as.numeric' wrappers in the code are to get rid of matrix or list structures
  # I don't know if they are all really necessary but this appears to work.
  #

  #
  # In the iteration the FOCs for component mu and sig2 do not allow a closed solution
  # if the variance components are not equal.
  # Note that within sig.fun.g I am using lots of variables defined outside its context (ugly but fast).
  # It is therefore necessary to (re)define the function within the blok function. Otherwise it would
  # always refer to global variables defined outside the blok context and not respond to iterated changes.
  #
  
  #specify function that solves for mu_g_t+1 and s2_g_t+1 given alpha_g and mu_g_t

  mu.fun.g <- function(alpha,i) {
    as.numeric(sum(tau[,i]*y/alpha)/sum(tau[,i]/alpha))
  }

  sig.fun.g <- function(alpha,mu.t,i) {
    as.numeric(sum(tau[,i]*((y-mu.t)^2)/alpha^2)/sum(tau[,i]/alpha))
  }
  
  for(k in 1:niter) {
    
    # print(sig)

                                        # E-step
    sig.sq <- sig^2
    tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*dnorm(y[j],m=mu[i],s=sqrt(sig.sq[i]+v.eps[j])))
    llhm <- log(tau)
    tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
    tau <- tau/tau.rs
    # warning:
    # this works because of the way R stores arrays and the dangerous 'recycling rule'
    # safer would be to define tau.rs by something like
    # tau.rs <- matrix(rep(tau %*% rep(1,g),each=g),ncol=g,byrow=T)

    llh <- sum(tau * llhm)
    progr.llh <<- rbind(progr.llh,c(p,mu,sig,llh)) # note the double arrow: this changes global variable progr.llh
    
    s2ji <- outer(1:n,1:g,FUN=function(j,i) sig.sq[i]+v.eps[j]) # allowing for heteroskedastic observations within groups

                                        # M-step

    s.vec.n <- rep(1,n)
    
    # for p there are closed-form solutions
    p <- as.numeric(s.vec.n %*% tau) / n

    #mu <- as.numeric(t(tau/s2ji) %*% y / (t(tau/s2ji)%*% s.vec.n ))
    #y.min.mu.sq <- outer(1:n,1:g,FUN=function(j,i) (y[j]-mu[i])^2)
    
    # for mu_g and s2_g FOCs amount to two equations for each component that need to be solved simultaneously
    for(i in 1:g) {
      #print("component")
	#print(i)
	mu.t=mu[i]
	alpha=1+(v.eps/sig.sq[i])
	for(j in 1:inner.iter) {
	  mu.t<-mu.fun.g(alpha,i)
	  #print("mu")
	  #print(mu.t)
	  s2.t<-sig.fun.g(alpha,mu.t,i)
	  #print("s2")
	  #print(s2.t)
	  alpha<-1+(v.eps/s2.t)
	  #print("alpha")
	  #print(alpha[1])
      }
	mu[i]<-mu.t
	sig.sq[i]<-s2.t
    }
    sig <- sqrt(sig.sq)
    if(printen) {
      print(c(p=p,mu=mu,sig=sig)) # put more diagnostics here if necessary
    }
  }
  c(p=p,mu=mu,sig=sig) # value of blok function
}

#
# Read Philippine data 2009
#
# Missing data was handled in Stata (data is free of missing values for variables included in the model)
library(foreign)
library(lfe)
data <- read.dta("c:/WB/EU/!Research/phi2009.dta")

# Administrative data are numeric in dataset. Change type into "factors". (To facilitate fixed-effects regression.) "regn" is already a "factor".
data$prov <- as.factor(data$prov)
data$id_mun <- as.factor(data$id_mun)
data$id_bgy <- as.factor(data$id_bgy)

# Run linear regression without FE
lmreg <- lm(lny~1+urban+hh_logsize+hd_female+hd_logage+hd_married+hd_unemp+hd_olf+hd_edu4+hd_edu5+hd_edu6+min1_emp+min1_edu56+pch_5+pch_516+peld_60+dw_basic+urb_logsize+urb_hd_fem+urb_hd_unemp+urb_hd_olf+urb_hd_edu5+urb_hd_edu6+urb_min_emp+urb_min_edu56+urb_pch5+urb_pch516+urb_peld60+urb_dw_basic,data)
eps <- resid(lmreg)

# ----------------------- Barangay level -----------------------

# Get bgy means
#epsreg_bgy <- felm(eps~0+G(id_bgy),data=data)
#ubar <- getfe(epsreg_bgy,se=FALSE)$effect
### CHECK if aggregating "eps" yields the same result
ubar <- aggregate(eps,by=list(data$id_bgy),FUN=mean)$x # don't ask; check documentation of aggregate() function
#hist(ubar,breaks=50)

# With bgy FE (to get H3-estimate of var[eps])
fereg_bgy <- felm(lny~1+urban+hh_logsize+hd_female+hd_logage+hd_married+hd_unemp+hd_olf+hd_edu4+hd_edu5+hd_edu6+min1_emp+min1_edu56+pch_5+pch_516+peld_60+dw_basic+urb_logsize+urb_hd_fem+urb_hd_unemp+urb_hd_olf+urb_hd_edu5+urb_hd_edu6+urb_min_emp+urb_min_edu56+urb_pch5+urb_pch516+urb_peld60+urb_dw_basic+G(id_bgy),data=data)
v.eps.bgy <- summary(fereg_bgy)$rse^2 # Assuming that the appropriate degrees of freedom correction is applied

# Get var(eps)/k_a
v.eps <- v.eps.bgy / as.numeric(table(data$id_bgy))

#
# ---------------------------------- Fit the mixture-distribution for F_v ----------------------------
#

# Initialize data
y <- ubar
set.seed(153)
n <- length(y) # should match m (i.e. number of clusters)

######## General approach (estimating component parameters jointly with mixing probabilities) ########

# Initialize parameters
g <- 4 # (4, "153" 2500, iter=1400) (5, "4925" 3500) 4 = best

mu <- rnorm(g,s=var(y)^0.5)
sig2 <- rnorm(g,s=var(y)^0.5)^2
sig <- sqrt(sig2)
p <- rep(1/g,g) # start with equal weights

# calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*dnorm(y[j],m=mu[i],s=sqrt(sig[i]^2+v.eps[j])))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <- c(p,mu,sig,llh) # collect steps in EM algorithm in progr.llh
p.c <- paste("p",1:g,sep="")
mu.c <- paste("mu",1:g,sep="")
sig.c <- paste("sig",1:g,sep="")
names(progr.llh) <- c(p.c,mu.c,sig.c,"llh")
options(width=190) # adjust width of output lines if you want (especially when using 4 or 5 components)

# Apply EM-algorithm (note that 1400 iterations needed for convergence)
psi.v.general <- fv.blok.general(1400);psi.v.general # psi contains final parameter estimates

p.v.general <- psi.v.general[1:g]
mu.v.general <- psi.v.general[(g+1):(2*g)]
sig.v.general <- psi.v.general[(2*g+1):(3*g)]

# Calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p.v.general[i]*dnorm(y[j],m=mu.v.general[i],s=sqrt(sig.v.general[i]^2+v.eps[j])))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <-  rbind(progr.llh,c(p.v.general,mu.v.general,sig.v.general,llh)) # add last solution to progr.llh data frame
psi.general
plot(progr.llh[,3*g+1],type="l") # plot progress of likelihood to see if it is nondecreasing

lines(sort(progr.llh[,3*g+1]),col="red") # plot llh progress in increasing order
psi.v.general <- progr.llh[order(progr.llh[,3*g+1]),][dim(progr.llh)[1],1:(3*g)] # llh maximizing parameters
p.v.general <- psi.v.general[1:g]
mu.v.general <- psi.v.general[(g+1):(2*g)]
sig.general <- psi.v.general[(2*g+1):(3*g)]

# construct mixture distribution for Fv using parameter estimates
g.v.general<-g
fv.dens.general <- function(x) as.numeric(sapply(x, function(t) psi.v.general[1:g.v.general] %*% dnorm(t,m=psi.v.general[(g.v.general+1):(2*g.v.general)],s=psi.v.general[(2*g.v.general+1):(3*g.v.general)])))

# compare to normal distribution with same overall mean and variance
#m.v.dens <- psi.v.general[1:g.v.general] %*% psi.v.general[(g.v.general+1):(2*g.v.general)]
#s.v.dens <- sqrt(psi.v.general[1:g.v.general] %*% (psi.v.general[(g.v.general+1):(2*g.v.general)]^2 + psi.v.general[(2*g.v.general+1):(3*g.v.general)]^2))


############### Restricted approach (estimating mixing probabilities only) ################

g <- 7 # we expect to need a larger number of components for the restricted approach

murange <- quantile(y,prob=c(0.1,0.9))
mu <- seq(murange[1],murange[2],by=(murange[2]-murange[1])/(g-1))
sig<-rep(var(y)/g,g)^0.5
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)
# alternatively set equal probabilities
#p <- rep(1/g,g)

# calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*dnorm(y[j],m=mu[i],s=sqrt(sig[i]^2+v.eps[j])))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <- c(p,mu,sig,llh) # collect steps in EM algorithm in progr.llh
p.c <- paste("p",1:g,sep="")
mu.c <- paste("mu",1:g,sep="")
sig.c <- paste("sig",1:g,sep="")
names(progr.llh) <- c(p.c,mu.c,sig.c,"llh")
options(width=190) # adjust width of output lines if you want (especially when using 4 or 5 components)

# Apply EM-algorithm (note that 5000 iterations needed for convergence)
psi.v.restricted <- fv.blok.restricted(5000);psi.v.restricted # psi contains final parameter estimates
 
p.v.restricted <- psi.v.restricted[1:g]
mu.v.restricted <- psi.v.restricted[(g+1):(2*g)]
sig.v.restricted <- psi.v.restricted[(2*g+1):(3*g)]

# Calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p.v.restricted[i]*dnorm(y[j],m=mu.v.restricted[i],s=sqrt(sig.v.restricted[i]^2+v.eps[j])))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <-  rbind(progr.llh,c(p.v.restricted,mu.v.restricted,sig.v.restricted,llh)) # add last solution to progr.llh data frame
plot(progr.llh[,3*g+1],type="l") # plot progress of likelihood to see if it is nondecreasing

# construct mixture distribution for Fv using parameter estimates
g.v.restricted<-g
fv.dens.restricted <- function(x) as.numeric(sapply(x, function(t) psi.v.restricted[1:g.v.restricted] %*% dnorm(t,m=psi.v.restricted[(g.v.restricted+1):(2*g.v.restricted)],s=psi.v.restricted[(2*g.v.restricted+1):(3*g.v.restricted)])))

# compare to normal distribution with same overall mean and variance
m.v.dens <- psi.v.restricted[1:g.v.restricted] %*% psi.v.restricted[(g.v.restricted+1):(2*g.v.restricted)]
s.v.dens <- sqrt(psi.v.restricted[1:g.v.restricted] %*% (psi.v.restricted[(g.v.restricted+1):(2*g.v.restricted)]^2 + psi.v.restricted[(2*g.v.restricted+1):(3*g.v.restricted)]^2))

#
# Compare plots of f_v between restricted and general approach

# left and right points for plots
plot.low <- quantile(y,0.01)
plot.high <- quantile(y,0.99)

plot(fv.dens.general,xlim=c(plot.low,plot.high),col="red")
plot(fv.dens.restricted,xlim=c(plot.low,plot.high),col="black",add=T)
plot(function(x) dnorm(x,m=m.v.dens,s=s.v.dens),xlim=c(plot.low,plot.high),add=T,col="blue")

#
# ---------------------------------- Fit the mixture-distribution for F_ebar ----------------------------
#

# Initialize data
y <- ubar
set.seed(153)
n <- length(y) # should match m (i.e. number of clusters)
v.eps<-0*v.eps # we will not estimate the distribution of ebar = v + epsbar, which is data that we observe, so we set the variance of the "measurement error" to zero

############### Restricted approach (estimating mixing probabilities only) ################

g <- 7 # we expect to need a larger number of components for the restricted approach

murange <- quantile(y,prob=c(0.1,0.9))
mu <- seq(murange[1],murange[2],by=(murange[2]-murange[1])/(g-1))
sig<-rep(var(y)/g,g)^0.5
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)

# calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*dnorm(y[j],m=mu[i],s=sqrt(sig[i]^2+v.eps[j])))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <- c(p,mu,sig,llh) # collect steps in EM algorithm in progr.llh
p.c <- paste("p",1:g,sep="")
mu.c <- paste("mu",1:g,sep="")
sig.c <- paste("sig",1:g,sep="")
names(progr.llh) <- c(p.c,mu.c,sig.c,"llh")
options(width=190) # adjust width of output lines if you want (especially when using 4 or 5 components)

# Apply EM-algorithm (note that 5000 iterations needed for convergence)
psi.ebar.restricted <- fv.blok.restricted(5000);psi.ebar.restricted # psi contains final parameter estimates
 
p.ebar.restricted <- psi.ebar.restricted[1:g]
mu.ebar.restricted <- psi.ebar.restricted[(g+1):(2*g)]
sig.ebar.restricted <- psi.ebar.restricted[(2*g+1):(3*g)]

# Calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p.ebar.restricted[i]*dnorm(y[j],m=mu.ebar.restricted[i],s=sqrt(sig.ebar.restricted[i]^2+v.eps[j])))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <-  rbind(progr.llh,c(p.ebar.restricted,mu.ebar.restricted,sig.ebar.restricted,llh)) # add last solution to progr.llh data frame
plot(progr.llh[,3*g+1],type="l") # plot progress of likelihood to see if it is nondecreasing

# construct mixture distribution for Fv using parameter estimates
g.ebar.restricted<-g
febar.dens.restricted <- function(x) as.numeric(sapply(x, function(t) psi.ebar.restricted[1:g.ebar.restricted] %*% dnorm(t,m=psi.ebar.restricted[(g.ebar.restricted+1):(2*g.ebar.restricted)],s=psi.ebar.restricted[(2*g.ebar.restricted+1):(3*g.ebar.restricted)])))

#
# ---------------------------------- Fit the mixture-distribution for F_eps ----------------------------
#

# Initialize data
y <- eps
set.seed(153)
n <- length(y) # should match total number of households

gs.fv <- Vectorize(function(x,m.gs,s.gs) sum(dnorm(x,m=m.gs+mu.v.restricted,s=sqrt(s.gs^2+sig.v.restricted^2))*p.v.restricted))

############### Restricted approach (estimating mixing probabilities only) ################

g <- 5 # we expect to need a larger number of components for the restricted approach

murange <- quantile(y,prob=c(0.1,0.9))
mu <- seq(murange[1],murange[2],by=(murange[2]-murange[1])/(g-1))
sig<-rep(var(y)/g,g)^0.5
p <- dnorm(mu ,m=0,s=var(y)^0.5)
p <- p/sum(p)
# alternatively set equal probabilities
#p <- rep(1/g,g)

# calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p[i]*gs.fv(y[j],mu[i],sig[i]))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <- c(p,mu,sig,llh) # collect steps in EM algorithm in progr.llh
p.c <- paste("p",1:g,sep="")
mu.c <- paste("mu",1:g,sep="")
sig.c <- paste("sig",1:g,sep="")
names(progr.llh) <- c(p.c,mu.c,sig.c,"llh")
options(width=190) # adjust width of output lines if you want (especially when using 4 or 5 components)

# Apply EM-algorithm (note that 1000 iterations needed for convergence)
psi.eps.restricted <- feps.blok.restricted(300);psi.eps.restricted # psi contains final parameter estimates
 
p.eps.restricted <- psi.eps.restricted[1:g]
mu.eps.restricted <- psi.eps.restricted[(g+1):(2*g)]
sig.eps.restricted <- psi.eps.restricted[(2*g+1):(3*g)]

# Calculate log likelihood
tau <- outer(1:n,1:g,FUN=function(j,i) p.eps.restricted[i]*gs.fv(y[j],mu.eps.restricted[i],sig.eps.restricted[i]))
tau.rs <- as.numeric(tau %*% rep(1,g)) # row sums
llhm <- log(tau)
tau <- tau/tau.rs
llh <- sum(tau * llhm);llh
progr.llh <-  rbind(progr.llh,c(p.eps.restricted,mu.eps.restricted,sig.eps.restricted,llh)) # add last solution to progr.llh data frame
plot(progr.llh[,3*g+1],type="l") # plot progress of likelihood to see if it is nondecreasing

# construct mixture distribution for F_eps using parameter estimates
g.eps.restricted<-g
feps.dens.restricted <- function(x) as.numeric(sapply(x, function(t) psi.eps.restricted[1:g.eps.restricted] %*% dnorm(t,m=psi.eps.restricted[(g.eps.restricted+1):(2*g.eps.restricted)],s=psi.eps.restricted[(2*g.eps.restricted+1):(3*g.eps.restricted)])))

# compare to normal distribution with same overall mean and variance
m.eps.dens <- psi.eps.restricted[1:g.eps.restricted] %*% psi.eps.restricted[(g.eps.restricted+1):(2*g.eps.restricted)]
s.eps.dens <- sqrt(psi.eps.restricted[1:g.eps.restricted] %*% (psi.eps.restricted[(g.eps.restricted+1):(2*g.eps.restricted)]^2 + psi.eps.restricted[(2*g.eps.restricted+1):(3*g.eps.restricted)]^2))

# Compare plot of f_eps to normal with identical moments

# left and right points for plots
plot.low <- quantile(y,0.01)
plot.high <- quantile(y,0.99)

plot(feps.dens.restricted,xlim=c(plot.low,plot.high),col="black")
plot(function(x) dnorm(x,m=m.eps.dens,s=s.eps.dens),xlim=c(plot.low,plot.high),add=T,col="blue")



### ------------------ sample from the estimated unconditional distribution f(v) ---------------------

# fv.dens.restricted (unconditional density for v)

# use "rejection sampling" to sample from the unconditional distribution
my.domain.fun <- function(x) x>-1 & x < 1
my.ldens <- function(x) log(fv.dens.restricted(x))

# sample 5000 observations
#v <- matrix(1,ncol=5000) # 5000 perhaps a bit much
#v[1,] <- arms(runif(1,-1,1), function(x) my.ldens(x), function(x) my.domain.fun(x), 5000)
v <- arms(runif(1,-1,1), function(x) my.ldens(x), function(x) my.domain.fun(x), 5000)

# compare histogram of sampled data to analytic density
hist(v,prob=T,n=50)
plot(fv.dens.restricted,xlim=c(-0.7,0.7),col="red",add=T)

###
### the same approach may be adopted to sample from the unconditional distribution f(eps)
###


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
#cv <- matrix(1,ncol=5000) # 5000 perhaps a bit much
#cv[1,] <- arms(runif(1,-1,1), function(x,u) my.con.ldens(x,u), function(x,u) my.domain.fun(x), 5000, u=ebar)
cv <- arms(runif(1,-1,1), function(x,u) my.con.ldens(x,u), function(x,u) my.domain.fun(x), 5000, u=ebar)

# compare histogram of sampled data to analytic density
hist(cv,prob=T,n=50)
plot(fv.con.dens1,xlim=c(-0.7,0.7),col="red",add=T)
plot(fv.dens.restricted,xlim=c(-0.7,0.7),col="blue",add=T)



### ------------------ evaluate conditional expectation ------------------------- (not needed for PovMap)

# set intervals for deconvolution
etarange <- c(-0.7,0.7)
intervals <- seq(etarange[1],etarange[2],0.1);intervals # limit points of intervals
n.intervals <- length(intervals)-1 

# my.ldens denotes log of f(ebar_a | v_a) f(v_a), where x = v_a, and u = ebar_a
my.domain.fun <- function(x) x>-1 & x < 1
my.con.ldens <- function(x,u) log(dnorm(u,m=x,s=(v.eps.bgy/k)^0.5))+log(fv.dens.restricted(x))

mcv <- matrix(nrow=n.intervals+1,ncol=5000) # 5000 perhaps a bit much
for(i in 1:(n.intervals+1)) {
  mcv[i,] <- arms(runif(1,-1,1), function(x,u) my.con.ldens(x,u), function(x,u) my.domain.fun(x), 5000, u=intervals[i])
}

# the benchmark line represents the conditional expectation under normal distribution
plot(intervals,rowMeans(mcv),ylim=c(-0.5,0.5)) # conditional expectation with deconvoluted v-distribution
vsigma2<-(psi.v.restricted[1:g.v.restricted] %*% (psi.v.restricted[(g.v.restricted+1):(2*g.v.restricted)]^2 + psi.v.restricted[(2*g.v.restricted+1):(3*g.v.restricted)]^2))
abline(a=0,b=1-(v.eps.bgy/k)/(vsigma2+(v.eps.bgy/k)))

