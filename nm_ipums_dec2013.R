
rm(list=ls())
library(foreign)
#library(lfe)
## library(HI) # implements rejection sampling
library(psych) # for matrix functions
library(MASS)
##setwd("c:/WB/NM/Empirical Study")
## setwd("c:/WB/EU/!Research")
setwd("/Users/Chris/Documents/PovMap/USA IPUMS 2010 voor NM studie Folder")

#### FUNCTIONS USED

reset.pars.udist <- function(){
  mu.old <<- mu.ini <- seq(from=-(ncomp.uest/2.5)*sqrt(v.u),to=(ncomp.uest/2.5)*sqrt(v.u),length=ncomp.uest)
  ppi.ini <<- dnorm(seq(from=-(ncomp.uest/2.5),to=(ncomp.uest/2.5),length=ncomp.uest))
  ppi.ini <<- ppi.old <<- (ppi.ini <<- ppi.ini / sum(ppi.ini))
  sig.old <<- sig.ini <- rep(sqrt(v.u-ppi.ini%*%mu.ini^2),ncomp.uest)
  sig2.ini <<- sig2.old <- sig.ini^2
}

## i refers to area, k to component

piphifun <- function(t,y){
  return(ppi.old * dnorm(t,y-mu.old,sig.old))
}
  
taufun <- function(t,y){
  return(piphifun(t,y) / sum(piphifun(t,y)))
}

condens <- function(t,y){
  omeg <- ppi.old*dnorm(y,mu.old,sqrt(v.eps+sig2.old))
  omeg <- omeg / sum(omeg)
  alf <- v.eps/(v.eps+sig2.old)
  return(omeg%*%dnorm(t,alf*(y-mu.old),sqrt(alf)*sig.old))
}

run.prog.udist <- function(niter=20){
  s2ky <- mky <- pky <- rep(0,length(sig.old))%o%rep(0,length(y.sample))
  for(iter in 1:niter){
    cat("iteration: ",iter,"\n")
    ## y.sample <- ybar               
    s2ky <- mky <- pky <- rep(0,length(sig.old))%o%rep(0,length(y.sample))
                                        #
    for(i in 1:length(y.sample)){
      fun2 <- function(t) condens(t,y.sample[i])
      for(k in 1:length(sig.old)){
        fun1 <- function(t) {
          r <- taufun(t,y.sample[i])[k]
          if(is.nan(r)) return(0) else return(r)
        }
        fun <- Vectorize(function(t) fun1(t)*fun2(t))
        pky[k,i] <- integrate(fun,-5,5)$value
        fun1 <- function(t) {
          r <- (y.sample[i]-t) * taufun(t,y.sample[i])[k]
          if(is.nan(r)) return(0) else return(r)
        }
        fun <- Vectorize(function(t) fun1(t)*fun2(t))
        mky[k,i] <- integrate(fun,-5,5)$value
        fun1 <- function(t) {
          r <- (y.sample[i]-t-mu.old[k])^2* taufun(t,y.sample[i])[k]
          if(is.nan(r)) return(0) else return(r)
        }
        fun <- Vectorize(function(t) fun1(t)*fun2(t))
        s2ky[k,i] <- integrate(fun,-5,5)$value    
      }
    }
                                        #
    ## plot(y.sample,pky[5,])
                                        #
    s2ky.p <- mky.p <- pky.p <- rep(0,length(sig.old))
    for(i in 1:length(sig.old)){
      pky.p[i] <- predict(interpSpline(y.sample,pky[i,]),ybar)$y %*% rep(1,length(ybar))
      mky.p[i] <- predict(interpSpline(y.sample,mky[i,]),ybar)$y %*% rep(1,length(ybar))
      s2ky.p[i] <- predict(interpSpline(y.sample,s2ky[i,]),ybar)$y %*% rep(1,length(ybar))
    }
                                        #
    ## check!
    ##    print(c(A=length(ybar),spkyp=sum(pky.p)))
    ppi.alt <- pky.p/length(ybar)
    mu.alt <- mky.p/pky.p ## length(ybar)
    s2.alt <- s2ky.p/pky.p ## length(ybar)
                                        #
    print(data.frame(ppi.ini=ppi.ini,ppi.old=ppi.old,ppi.new=ppi.alt))
    ppi.old <<- ppi.alt
    print(data.frame(mu.ini=mu.ini,mu.old=mu.old,mu.new=mu.alt))
    mu.old <<- mu.alt
    print(data.frame(sig.ini=sig.ini,sig.old=sig.old,sig.new=sqrt(s2.alt)))
    sig.old <<- sqrt(s2.alt)
    sig2.old <<- s2.alt
  }
  return(list(ppi.new=ppi.alt,mu.new=mu.alt,sig.new=sqrt(s2.alt),sig2.new=s2.alt))
} ## end run.prog.udist

reset.pars.eps <- function(){
  nu.old <<- nu.ini <- seq(from=-(ncomp.epsest/2.5)*sqrt(v.eps),to=(ncomp.epsest/2.5)*sqrt(v.eps),length=ncomp.epsest)
  lam.ini <<- dnorm(seq(from=-(ncomp.epsest/2.5),to=(ncomp.epsest/2.5),length=ncomp.epsest))
  lam.ini <- lam.old <<- (lam.ini <<- lam.ini / sum(lam.ini))
  tau.old <<- tau.ini <- rep(sqrt(v.eps-lam.ini%*%nu.ini^2),ncomp.epsest)
  tau2.ini <<- tau2.old <- tau.ini^2
}

piphifun.eps <- function(t,y){
  return(lam.old * dnorm(t,y-nu.old,tau.old))
}
  
taufun.eps <- function(t,y){
  return(piphifun.eps(t,y) / sum(piphifun.eps(t,y)))
}

condens.u <- function(t,y){ ## this is the approximated conditional distribution of u (or eta), given y (i.e. ebar)
  ## note that this function does not change during the iterations
  omeg <- ppi.u*dnorm(y,mu.u,sqrt(v.epsbar+sig2.u))
  omeg <- omeg / sum(omeg)
  alf <- sig2.u/(v.epsbar+sig2.u)
  return(omeg%*%dnorm(t,alf*y+(1-alf)*mu.u,sqrt(alf)*sqrt(v.epsbar)))
}


run.prog.epsdist <- function(niter=21){
  s2kyy <- mkyy <- pkyy <- array(0,dim=c(length(tau.old),length(ybar.sample),length(ybar.sample)),dimnames=c("cmp","yi","ybar"))
  for(iter in 1:niter){
    cat("iteration: ",iter,"\n")
## y.sample <- ybar               
    s2kyy <- mkyy <- pkyy <- array(0,dim=c(length(tau.old),length(y.sample),length(y.sample)),dimnames=c("cmp","yi","ybar"))
    for(i1 in 1:n.ybar.sample){
      fun2 <- function(t) condens.u(t,ybar.sample[i1]) ## in principle fun2 could be taken out of the for(iter...)-loop
      for(i2 in 1:n.y.sample){
        for(k in 1:length(tau.old)){
          fun1 <- function(t) {
            r <- taufun.eps(t,y.sample[i2])[k]
            if(is.nan(r)) return(0) else return(r)
          }
          fun <- Vectorize(function(t) fun1(t)*fun2(t))
          pkyy[k,i2,i1] <- integrate(fun,-5,5)$value
          fun1 <- function(t) {
            r <- (y.sample[i2]-t) * taufun.eps(t,y.sample[i2])[k]
            if(is.nan(r)) return(0) else return(r)
          }
          fun <- Vectorize(function(t) fun1(t)*fun2(t))
          mkyy[k,i2,i1] <- integrate(fun,-5,5)$value
          fun1 <- function(t) {
            r <- (y.sample[i2]-t-nu.old[k])^2* taufun.eps(t,y.sample[i2])[k]
            if(is.nan(r)) return(0) else return(r)
          }
          fun <- Vectorize(function(t) fun1(t)*fun2(t))
          s2kyy[k,i2,i1] <- integrate(fun,-5,5)$value    
        }
      }
    }
## plot(y.sample,pky[5,])
  
    s2ky.p <- mky.p <- pky.p <- matrix(0,nrow=length(tau.old),ncol=length(ybar))
    for(a in unique(ydata$area)){
      yvals <- ydata$y[ydata$area==a]
      ybar.a <- ybardata$ybar[ybardata$area==a]
      for(k in 1:length(tau.old)){
        pky.p[k,a] <- sum(bicubic(y.sample,ybar.sample,pkyy[k,1:n.y.sample,1:n.ybar.sample],yvals,rep(ybar.a,length(yvals)))$z) 
        mky.p[k,a] <- sum(bicubic(y.sample,ybar.sample,mkyy[k,1:n.y.sample,1:n.ybar.sample],yvals,rep(ybar.a,length(yvals)))$z) 
        s2ky.p[k,a] <- sum(bicubic(y.sample,ybar.sample,s2kyy[k,1:n.y.sample,1:n.ybar.sample],yvals,rep(ybar.a,length(yvals)))$z) 
      }
    }

    lam.alt <- pky.p%*%rep(1,length(ybar))
    nu.alt <- mky.p%*%rep(1,length(ybar)) / lam.alt
    t2.alt <- s2ky.p%*%rep(1,length(ybar)) / lam.alt
    lam.alt <- lam.alt / sum(lam.alt)
  
## check!
## print(c(A=length(ydata$y),spkyp=sum(pky.p)))

    print(data.frame(lam.ini=lam.ini,lam.old=lam.old,lam.new=lam.alt))
    lam.old <<- lam.alt
    print(data.frame(nu.ini=nu.ini,nu.old=nu.old,nu.new=nu.alt))
    nu.old <<- nu.alt
    print(data.frame(tau.ini=tau.ini,tau.old=tau.old,tau.new=sqrt(t2.alt)))
    tau.old <<- sqrt(t2.alt)
    tau2.old <<- t2.alt
  }
  return(list(lam.new=lam.alt,nu.new=nu.alt,tau.new=sqrt(t2.alt),tau2.new=t2.alt))
} # end run.prog.epsdist

### poverty: head-count
est.pov.ind <- function(x,z=pov.line,ppi=ppi.new,mu=mu.new,sig=sig.new,
                    lam=lam.new,nu=nu.new,tau=tau.new,beta=beta.gls){
  ## x does not include the intercept
  yhat <- beta[1] + x %*% beta[-1]
  ppilam <- as.numeric(ppi %o% lam)
  munu <- as.numeric(outer(mu,nu,FUN="+"))
  sigtau <- as.numeric(sqrt(outer(sig^2,tau^2,"+")))
  return(as.numeric(ppilam %*% pnorm(z-yhat,munu,sigtau)))
}

est.pov.ind.normal <- function(x,z=pov.line,sig2.eps=v.eps,sig2.u=v.u,beta=beta.gls){
  ## x does not include the intercept
  yhat <- beta[1] + x %*% beta[-1]
  return(pnorm(z-yhat,0,sqrt(sig2.eps+sig2.u)))
}

c.est.pov.ind <- function(x,ar,z=pov.line,ppi=ppi.new,mu=mu.new,sig=sig.new,
                    lam=lam.new,nu=nu.new,tau=tau.new,beta=beta.gls){
  ## x does not include the constant term
  ## ar is area indicator
  y <- ybardata$y[ybardata$area==ar] ## or ybardata$ybar[ybardata$area==ar]?
  omeg <- ppi*dnorm(y,mu,sqrt(v.epsbar+sig^2))
  omeg <- omeg / sum(omeg)
  alf <- sig^2/(v.epsbar+sig^2)
  cmu.u <- alf*y+(1-alf)*mu
  csig2.u <- v.epsbar+sig^2
  csig.u <- sqrt(csig2.u)
  return(est.pov.ind(x,pov.line,omeg,cmu.u,csig.u,lam,nu,tau,beta))
##  return(list(omeg=omeg,cmu=cmu.u,csig=csig.u,csig2=csig2.u,exp.pov=est.pov.ind(x,pov.line,omeg,cmu.u,csig.u,lam,nu,tau,beta)))
}

c.est.pov.ind.normal <- function(x,ar,beta=beta.gls){
  ## x does not include the constant term
  ## ar is area indicator
  return(c.est.pov.ind(x,ar,pov.line,1,0,sqrt(v.u),1,0,sqrt(v.eps),beta))
}

### inequality: MLD
est.mld.ind <- function(x,ppi=ppi.new,mu=mu.new,sig=sig.new,
                    lam=lam.new,nu=nu.new,tau=tau.new,beta=beta.gls){
  ## x does not include the intercept
  yhat <- beta[1] + x %*% beta[-1]
  #ppilam <- as.numeric(ppi %o% lam)
  #munu <- as.numeric(outer(mu,nu,FUN="+"))
  #sigtau2 <- as.numeric(outer(sig^2,tau^2,"+"))
  fbar <- as.numeric(lam %*% exp((nu+0.5*(tau^2))))
  mbar <- as.numeric(lam %*% nu)
  return(c(fbar*exp(yhat),yhat+mbar))
}

est.mld.ind.normal <- function(x,sig2.eps=v.eps,sig2.u=v.u,beta=beta.gls){
  ## x does not include the intercept
  yhat <- beta[1] + x %*% beta[-1]
  return(c(exp(0.5*(sig2.eps))*exp(yhat),yhat))
}

c.est.mld.ind <- function(x,ar,ppi=ppi.new,mu=mu.new,sig=sig.new,
                    lam=lam.new,nu=nu.new,tau=tau.new,beta=beta.gls){
  ## x does not include the constant term
  ## ar is area indicator
  y <- ybardata$y[ybardata$area==ar] ## or ybardata$ybar[ybardata$area==ar]?
  omeg <- ppi*dnorm(y,mu,sqrt(v.epsbar+sig^2))
  omeg <- omeg / sum(omeg)
  alf <- sig^2/(v.epsbar+sig^2)
  cmu.u <- alf*y+(1-alf)*mu
  csig2.u <- v.epsbar+sig^2
  csig.u <- sqrt(csig2.u)
  return(est.mld.ind(x,omeg,cmu.u,csig.u,lam,nu,tau,beta))
}

c.est.mld.ind.normal <- function(x,ar,beta=beta.gls){
  ## x does not include the constant term
  ## ar is area indicator
  return(c.est.mld.ind(x,ar,1,0,sqrt(v.u),1,0,sqrt(v.eps),beta))
}

msdiff <- function(v1,v2){ ## mean square difference
  mean((v1-v2)^2)
}

all.functions <- ls()

######################## main program ##########################
#
#------------------------------ Read survey data ----------------------
# (data is free of missing values for variables included in the model)
#
############# step one: data, model, etc.
#
sur.data <- read.dta("ipums2010_survey.dta")
##sur.data <- read.csv("/Users/Chris/Documents/PovMap/sim3_survey_ah.csv",header=T)
##sur.data.u <- read.csv("sim_dagum_survey_a_rho005_shape05.csv",header=T)
##names(sur.data)[1] <- "hhid"
##sur.data.u <- sur.data.u[,2:3]

sur.data$cluster <- sur.data$id_county
clusterSize<-as.vector(table(sur.data$cluster)) ## assumed here: equal cluster size throughout
#write(clusterSize,file="clusterSize.txt")

sur.data$y.ah <- sur.data$lnycap
reg.beta <- lm(y.ah~metro_suburb+metro_none+hh_size_1+hh_size_2+hh_size_3+share_child15+hh_hedu2+hh_employed+hd_female+hd_logage+hd_ledu+hd_hedu1+hd_hedu2+hd_employed+hd_hisp+hd_black,data=sur.data)
summary(reg.beta)
r_hat <- resid(reg.beta)
n_hlds <- length(r_hat)

##
## Henderson III
##

X.ols <- model.matrix(reg.beta)
X <- as.matrix(X.ols[,-1])

# number of independent variables (variables in X excluding the constant)
d0 <- dim(X)[2]

# and including the constant
d <- d0+1
d.ols <- dim(X.ols)[2]


# --------------------------- cluster level -----------------------
# Get cluster means
ubar <- aggregate(r_hat,by=list(sur.data$cluster),FUN=mean)$x
n_area <- length(ubar)

# H3 variance of epsilon: within estimator
ka <- as.numeric(table(sur.data$cluster))
##y <- data$lny
y <- sur.data$y.ah
Dka <- diag(ka)
ybar <- aggregate(sur.data$y,by=list(sur.data$cluster),FUN=mean)$x
a.X <- aggregate(X,by=list(sur.data$cluster),FUN=mean)
Xbar <- as.matrix(a.X[,2:dim(a.X)[2]])
XTX <- (t(X)%*%X)-(t(Xbar)%*%Dka%*%Xbar)
XTy <- (t(X)%*%y)-(t(Xbar)%*%Dka%*%ybar)
beta.within <- ginv(XTX)%*%XTy
eps_hat <- (y-rep(ybar,ka))-(X%*%beta.within-rep(Xbar%*%beta.within,ka))
sigma2.eps <- sum(eps_hat*eps_hat)/(n_hlds-(d-1)-n_area)

v.eps <- sigma2.eps / ka  ## warning: v.eps has different meaning further down!

# H3 variance of v: Corollary 3
a.X.ols <- aggregate(X.ols,by=list(sur.data$cluster),FUN=mean)
UX.ols <- as.matrix(a.X.ols[,2:dim(a.X.ols)[2]])
UX.ols <- Dka%*%UX.ols
trupxu <- tr(ginv(t(X.ols)%*%X.ols)%*%t(UX.ols)%*%UX.ols)
sigma2.v <- (n_hlds-d.ols)*((sum(r_hat*r_hat)/(n_hlds-d.ols))-sigma2.eps)/(n_hlds-trupxu)

#estimate of rho
rho <- sigma2.v/(sigma2.v+sigma2.eps)

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

res.gls <- y-X.ols%*%beta.gls

##
## switch to names of variables of the other script files
##

ydata <- data.frame(area=sur.data$id_county,hh=sur.data$serial,y=res.gls)
##ydata <- data.frame(area=sur.data$id_county,hh=sur.data$serial,y=res.gls,ytrue=y-X.ols%*%c(0,1))
## ydata <- cbind(ydata,model.matrix(reg.beta))
## ydata$"(Intercept)" <- NULL

##u.expand <- rep(sur.data.u$u.a,ka)
##ydata$u <- u.expand
##ydata$eps <- ydata$ytrue-u.expand

## clean up
rm(list=setdiff(ls(),c(all.functions,"all.functions","ydata","sigma2.eps","sigma2.v","reg.beta","beta.gls","sur.data"))) ## now only ydata, sigma2.eps etc. remain
detach(package:foreign)
detach(package:psych)
detach(package:MASS)

save.image(file="nmdec2013.RData")

####################### U DISTRIBUTION ################

library(splines)

#
# ---------------------------- Fit the mixture-distribution for F_v ----------------------------

## replace by new procedure
##
## ebar moet gemiddelde 0 hebben, en variantie volgens H3
##

ybardata <- aggregate(ydata[,c("y")],by=list(ydata$area),FUN=mean)
##ybardata <- aggregate(ydata[,c("u","eps","y","ytrue")],by=list(ydata$area),FUN=mean)
names(ybardata)[1] <- "area"
names(ybardata)[2] <- "y"

## centering the y:
ybardata$ybar <- ybardata$y-mean(ybardata$y)

## PARAMETER
##
## warning: v.eps has different meaning than earlier defined
##
v.eps <- sigma2.eps/(length(ydata$y)/length(ybardata$y)) ## expectation of eps is set to zero when estimating ubar distribution

## initial mu of u distribution; var u will come from Henderson-3 or similar in practice

## PARAMETER
##
v.u <- sigma2.v
## rm(list=c("sigma2.eps","sigma2.v"))

## PARAMETER
##
## ncomp.uest is the number of components in the distribution of u
##
## ncomp.uest <- 5
## ncomp.uest <- 2
#ncomp.uest <- 3
ncomp.uest <- 3
mu.old <- mu.ini <- seq(from=-(ncomp.uest/2.5)*sqrt(v.u),to=(ncomp.uest/2.5)*sqrt(v.u),length=ncomp.uest)
ppi.ini <- dnorm(seq(from=-(ncomp.uest/2.5),to=(ncomp.uest/2.5),length=ncomp.uest))
ppi.ini <- ppi.old <- (ppi.ini <- ppi.ini / sum(ppi.ini))
sig.old <- sig.ini <- rep(sqrt(v.u-ppi.ini%*%mu.ini^2),ncomp.uest)
sig2.ini <- sig2.old <- sig.ini^2

reset.pars.udist()
mu.ini
ppi.ini
sig.ini
sig2.ini

##
## using spline interpolation
##

## PARAMETER
##
## grid.size is the number of y-points for which the expectation is computed
##

grid.size <- 50

y.sample <- do.call("seq",as.list(c(range(ybardata$ybar),length.out=grid.size)))

ybar <- (ybardata$ybar)

reset.pars.udist()
(outcomes <- run.prog.udist(60))  ## 100 is the number of iterations
ppi.new <- outcomes$ppi.new
mu.new <- outcomes$mu.new
sig.new <- outcomes$sig.new
sig2.new <- outcomes$sig2.new

dens.uhat <- Vectorize(function(u) as.numeric(ppi.new%*% dnorm(u,mu.new,sig.new)))

## check expectation

## Hier gaat iets fout: ybardata$u bestaat niet.
## Moet dit zijn ybardata$y of ybardata$ybar? In dit geval maakt het niet uit.
## ...omdat survey cluster size constant is????
## Volgens comment line hieronder moet het ybardata$y zijn

c(mean.u=mean(ybardata$y),var.u=v.u,mean.uhat=ppi.new%*%mu.new,var.uhat=ppi.new%*%(mu.new^2+sig2.new))


## check densities
## NOTE: ydata$u would denote the true area error data! (in simulations)

ybar.rescaled <- ybar*(sigma2.v^0.5)/sd(ybar)

curve(dens.uhat,min(ybar.rescaled),max(ybar.rescaled),col="red")
## curve(dens.u,-3,3,add=TRUE,col="red")
curve(dnorm(x,mean=0,sd=sigma2.v^0.5),min(ybar.rescaled),max(ybar.rescaled),add=TRUE,col="blue")
lines(density(ybar.rescaled,from=min(ybar.rescaled),to=max(ybar.rescaled))) ## since i don't have the true density


## save these on pdf
pdf(file="plot1.pdf")
curve(dens.uhat,min(ybar),max(ybar),col="red")
lines(density(ybar,n=200,from=min(ybar),to=max(ybar))) ## since i don't have the true density
dev.off()

## clean up

rm(list=setdiff(ls(),c(all.functions,"all.functions","mu.new","ppi.new","sig.new","sig2.new","v.eps","v.u","ybardata","ydata","sigma2.eps","sigma2.v","reg.beta","beta.gls","sur.data","sur.data.u")))

################################### EPS DISTRIBUTION #####################
###### eps distribution ####

library(akima)

##
## Warning: another redefinition of v.eps
##
## PARAMETER
v.epsbar <- v.eps
## PARAMETER
## v.eps <- v.epsbar * length(ydata$y)/length(ybardata$y) ## sig
v.eps <- sigma2.eps

## distribution of u:
ppi.u <- ppi.new
mu.u <- mu.new
sig2.u <- sig2.new
sig.u <- sig.new

## PARAMETER
## choose number of components for eps distribution
## ncomp.epsest <- 4
ncomp.epsest <- 2

## PARAMETER
## 2.5 below can be changed but appears to work fine
##
nu.old <<- nu.ini <- seq(from=-(ncomp.epsest/2.5)*sqrt(v.eps),to=(ncomp.epsest/2.5)*sqrt(v.eps),length=ncomp.epsest)
lam.ini <<- dnorm(seq(from=-(ncomp.epsest/2.5),to=(ncomp.epsest/2.5),length=ncomp.epsest))
lam.ini <- lam.old <<- (lam.ini <- lam.ini / sum(lam.ini))
tau.old <<- tau.ini <- rep(sqrt(v.eps-lam.ini%*%nu.ini^2),ncomp.epsest)
tau2.ini <<- tau2.old <- tau.ini^2

reset.pars.eps()
nu.ini
lam.ini
tau.ini
tau2.ini


## prepare for iterations
## using bicubic spline interpolation provided by package akima
##
## PARAMETER
n.ybar.sample <- 12
## PARAMETER
n.y.sample <- 12
ybar.sample <- do.call("seq",as.list(c(range(ybardata$ybar),length.out=n.ybar.sample)))
y.sample <- do.call("seq",as.list(c(range(ydata$y),length.out=n.y.sample)))
ybar <- (ybardata$ybar)

reset.pars.eps()
(outcomes <- run.prog.epsdist(100)) ## again, use more iterations if you want

lam.new <- as.numeric(outcomes$lam.new)
nu.new <- as.numeric(outcomes$nu.new)
tau.new <- as.numeric(outcomes$tau.new)
tau2.new <- as.numeric(outcomes$tau2.new)

dens.epshat <- Vectorize(function(t) as.numeric(lam.new%*% dnorm(t,nu.new,tau.new)))

## check expectation
c(mean.eps=mean(ydata$eps),var.eps=v.eps,mean.epshat=lam.new%*%nu.new,var.epshat=lam.new%*%(nu.new^2+tau2.new))

## check densities
curve(dens.epshat,min(ydata$eps),max(ydata$eps),col="red")

##lines(density(ydata$eps,n=200,from=min(ydata$eps),to=max(ydata$eps))) ## since i don't have the true density

## on pdf:
pdf(file="plot2.pdf")
curve(dens.epshat,min(ydata$eps),max(ydata$eps),col="red")
lines(density(ydata$eps,n=200,from=min(ydata$eps),to=max(ydata$eps))) ## since i don't have the true density
dev.off()

## just to show non-normality of the artificial "data"
qqnorm(ydata$u+ydata$eps,cex=0.1)
qqnorm(ydata$u,cex=0.1)
qqnorm(ydata$eps,cex=0.1)

############################################ SAMPLE IMPUTATION ###############
## prepare for census imputation
##
## you might want to remove some variables:
##
ls()

## PARAMETER
pov.line <- -0.5
## sample poverty rate: around 0.32
## > length(sur.data$y.ah[sur.data$y.ah < -0.35])/dim(sur.data)[1]
## [1] 0.3157333

## check:
## individual probability of poverty
##
## > mean(sapply(sur.data$x.ah,est.pov.ind))
## [1] 0.3130839
## > 

## individual estimated poverty and variance
##
estpov.ind <- sapply(sur.data$x.ah,est.pov.ind)
povvar.ind <- estpov.ind-estpov.ind^2

estpov.ind.normal <- sapply(sur.data$x.ah,est.pov.ind.normal)
povvar.ind.normal <- estpov.ind.normal-estpov.ind.normal^2
  
poor.index <- ifelse(sur.data$y.ah < -0.35,1,0)
## MSE

cat("MSE normal mixture (within sample):    ",mean((estpov.ind-poor.index)^2),"\n")
cat("MSE normal nomixture (within sample):  ",mean((estpov.ind.normal-poor.index)^2),"\n")

cat("Bias normal mixture (within sample):    ",mean((estpov.ind-poor.index)),"\n")
cat("Bias normal nomixture (within sample):  ",mean((estpov.ind.normal-poor.index)),"\n")


## Hardly any difference, but the two estimates are not completely equivalent:
plot(estpov.ind,estpov.ind.normal,cex=0.1)
pdf(file="plot3.pdf")
plot(estpov.ind,estpov.ind.normal,cex=0.1)
dev.off()

cestpov.ind <- sapply(1:(dim(sur.data)[1]),function(n) c.est.pov.ind(sur.data$x.ah[n],sur.data$area.id[n]))

cestpov.ind.normal <- sapply(1:(dim(sur.data)[1]),function(n) c.est.pov.ind.normal(sur.data$x.ah[n],sur.data$area.id[n]))

cat("MSE normal mixture (within sample, conditional):    ",mean((cestpov.ind-poor.index)^2),"\n")
cat("MSE normal nomixture (within sample, conditional):  ",mean((cestpov.ind.normal-poor.index)^2),"\n")
## conditioned on ebar the differences are non-trivial:
plot(cestpov.ind,cestpov.ind.normal,cex=0.1)
cestpov.ind.agg <- aggregate(cestpov.ind,by=list(ydata$area),FUN=mean)$x
cestpov.ind.normal.agg <- aggregate(cestpov.ind.normal,by=list(ydata$area),FUN=mean)$x
plot(cestpov.ind.agg,cestpov.ind.normal.agg,cex=0.1)
plot(cestpov.ind.agg[order(cestpov.ind.agg)],type="l")
lines(cestpov.ind.normal.agg[order(cestpov.ind.agg)],cex=0.1,col="red")
cat(" root mean square difference NM and N:  ",sqrt(mean((cestpov.ind.agg-cestpov.ind.normal.agg)^2)),"\n")
pdf(file="plot4.pdf")
plot(cestpov.ind,cestpov.ind.normal,cex=0.3)
plot(cestpov.ind.agg,cestpov.ind.normal.agg,cex=0.3)
plot(cestpov.ind.agg[order(cestpov.ind.agg)],type="l")
lines(cestpov.ind.normal.agg[order(cestpov.ind.agg)],cex=0.1,col="red")
dev.off()


################ CENSUS IMPUTATION ############
###
### Method 1: read in all data (if you have enough memory)
###
cen.data <- read.csv("sim_dagum_census_ah_rho005_shape05.csv",header=T)
cen.data.u <- read.csv("sim_dagum_census_a_rho005_shape05.csv",header=T)
##cen.data <- read.csv("/Users/Chris/Documents/PovMap/sim3_census_ah.csv",header=T)
cen.data.u <- sur.data.u

pov.line <- 0

###
### EB!
###
census.cestpov.ind <- sapply(1:(dim(cen.data)[1]),function(n) c.est.pov.ind(cen.data$x.ah[n],cen.data$area.id[n]))
census.cestpov.ind.normal <- sapply(1:(dim(cen.data)[1]),function(n) c.est.pov.ind.normal(cen.data$x.ah[n],cen.data$area.id[n]))

pov.true.ind <- ifelse(cen.data$y.ah<pov.line,1,0)
pov.pred <- aggregate(census.cestpov.ind,by=list(cen.data$area.id),FUN=mean)
names(pov.pred) <- c("area","pred.pov")
pov.true <- aggregate(pov.true.ind,by=list(cen.data$area.id),FUN=mean)
names(pov.true) <- c("area","true.pov")
pov.pred.normal <- aggregate(census.cestpov.ind.normal,by=list(cen.data$area.id),FUN=mean)
names(pov.pred.normal) <- c("area","pred.pov.normal")

###
### RMSE!!!
###
cat("POV - MSE MN prediction at area level:  ",msdiff(pov.pred$pred.pov,pov.true$true.pov),"\n")
cat("POV - MSE N prediction at area level:   ",msdiff(pov.pred.normal$pred.pov.normal,pov.true$true.pov),"\n")
cat("POV - RMSE MN prediction at area level:  ",sqrt(msdiff(pov.pred$pred.pov,pov.true$true.pov)),"\n")
cat("POV - RMSE N prediction at area level:   ",sqrt(msdiff(pov.pred.normal$pred.pov.normal,pov.true$true.pov)),"\n")
cat("POV - BIAS MN prediction at area level:  ",mean(pov.pred$pred.pov-pov.true$true.pov),"\n")
cat("POV - BIAS N prediction at area level:   ",mean(pov.pred.normal$pred.pov.normal-pov.true$true.pov),"\n")

### INEQUALITY: MLD
census.cestmld.ind <- sapply(1:(dim(cen.data)[1]),function(n) c.est.mld.ind(cen.data$x.ah[n],cen.data$area.id[n]))
census.cestmld.ind.normal <- sapply(1:(dim(cen.data)[1]),function(n) c.est.mld.ind.normal(cen.data$x.ah[n],cen.data$area.id[n]))
mld.true.ind <- cbind(exp(cen.data$y.ah),cen.data$y.ah)

mld.aggr <- aggregate(t(census.cestmld.ind),by=list(cen.data$area.id),FUN=mean)
mld.aggr.normal <- aggregate(t(census.cestmld.ind.normal),by=list(cen.data$area.id),FUN=mean)
mld.aggr.true <- aggregate(mld.true.ind,by=list(cen.data$area.id),FUN=mean)

mld.pred <- log(mld.aggr[2])-mld.aggr[3]
mld.pred.normal <- log(mld.aggr.normal[2])-mld.aggr.normal[3]
mld.true <- log(mld.aggr.true[2])-mld.aggr.true[3]

cat("MLD - RMSE MN prediction at area level:  ",sqrt(msdiff(mld.pred,mld.true)),"\n")
cat("MLD - RMSE N prediction at area level:   ",sqrt(msdiff(mld.pred.normal,mld.true)),"\n")
cat("MLD - BIAS MN prediction at area level:  ",mean(mld.pred$V1-mld.true$V1),"\n")
cat("MLD - BIAS N prediction at area level:   ",mean(mld.pred.normal$V1-mld.true$V1),"\n")

###
### POVERTY FIGURE!!!
###
inc.pov <- order(pov.true$true.pov)
#plot(pov.true$true.pov,pov.pred$pred.pov,cex=0.1)
#plot(pov.true$true.pov,pov.pred.normal$pred.pov.normal,cex=0.1)
plot(pov.true$true.pov[inc.pov],type="l")
lines(pov.pred.normal$pred.pov.normal[inc.pov],col="blue")
lines(pov.pred$pred.pov[inc.pov],col="red")

###
### NON-EB!
###
census.estpov.ind <- sapply(1:(dim(cen.data)[1]),function(n) est.pov.ind(cen.data$x.ah[n]))
census.estpov.ind.normal <- sapply(1:(dim(cen.data)[1]),function(n) est.pov.ind.normal(cen.data$x.ah[n]))

pov.true.ind <- ifelse(cen.data$y.ah<pov.line,1,0)
pov.pred <- aggregate(census.estpov.ind,by=list(cen.data$area.id),FUN=mean)
names(pov.pred) <- c("area","pred.pov")
pov.true <- aggregate(pov.true.ind,by=list(cen.data$area.id),FUN=mean)
names(pov.true) <- c("area","true.pov")
pov.pred.normal <- aggregate(census.estpov.ind.normal,by=list(cen.data$area.id),FUN=mean)
names(pov.pred.normal) <- c("area","pred.pov.normal")

###
### RMSE!!!
###
cat("MSE MN prediction at area level:  ",msdiff(pov.pred$pred.pov,pov.true$true.pov),"\n")
cat("MSE N prediction at area level:   ",msdiff(pov.pred.normal$pred.pov.normal,pov.true$true.pov),"\n")
cat("RMSE MN prediction at area level:  ",sqrt(msdiff(pov.pred$pred.pov,pov.true$true.pov)),"\n")
cat("RMSE N prediction at area level:   ",sqrt(msdiff(pov.pred.normal$pred.pov.normal,pov.true$true.pov)),"\n")
cat("BIAS MN prediction at area level:  ",mean(pov.pred$pred.pov-pov.true$true.pov),"\n")
cat("BIAS N prediction at area level:   ",mean(pov.pred.normal$pred.pov.normal-pov.true$true.pov),"\n")









t.sm <- loess(pov.pred$pred.pov[inc.pov] ~ I(1:500) )
t2.sm <- loess(pov.pred.normal$pred.pov.normal[inc.pov] ~ I(1:500) )
plot(pov.true$true.pov[inc.pov],type="l")
lines(t.sm$fitted,col="red")
lines(t2.sm$fitted,col="blue")

rm(cen.data)
rm(cen.data.u)


################ CENSUS IMPUTATION - 2 ############
###
### Method 2: read in all data in blocks, slower but needs less resources


## use connection to read data from a large file several lines at a time

cen.data <- file("../Shared with Chris/sim_census_ah.csv","r")
readLines(cen.data,n=1) # skip header line

bloxice <- 100000 ## number of lines to read each time from census file
doorlezen <- TRUE
agg.file <- data.frame(
                       area=numeric(0),
                       pov.true=numeric(0),
                       pov.uc.nm=numeric(0),
                       pov.uc.n=numeric(0),
                       pov.c.nm=numeric(0),
                       pov.c.n=numeric(0),
                       count=numeric(0)
                       )

teller=0
while(doorlezen){
  deeltje <- as.data.frame(scan(cen.data,what=list(num="",area=0,y=0,x=0,poor=0),sep=",",nlines=bloxice))
  if(dim(deeltje)[1] == 0) break
  if(dim(deeltje)[1] < bloxice) doorlezen <- FALSE
  d.estpov.ind <- sapply(deeltje$x,est.pov.ind)
  d.estpov.ind.normal <- sapply(deeltje$x,est.pov.ind.normal)
  d.cestpov.ind <- sapply(1:(dim(deeltje)[1]),function(n) c.est.pov.ind(deeltje$x[n],deeltje$area[n]))
  d.cestpov.ind <- sapply(1:(dim(deeltje)[1]),function(n) c.est.pov.ind(deeltje$x[n],deeltje$area[n]))
  d.cestpov.ind.normal <- sapply(1:(dim(deeltje)[1]),function(n) c.est.pov.ind.normal(deeltje$x[n],deeltje$area[n]))
  d.truepov <- ifelse(deeltje$y < pov.line,1,0)
  a.tmp <- data.frame(area=deeltje$area,pov.true=d.truepov,
                      pov.uc.nm=d.estpov.ind,pov.uc.n=d.estpov.ind.normal,
                      pov.c.nm=d.cestpov.ind,pov.c.n=d.cestpov.ind.normal)
  a.count <- aggregate(rep(1,dim(a.tmp)[1]),by=list(a.tmp$area),FUN=sum)
  a.tmp <- aggregate(a.tmp[,c("pov.true","pov.uc.nm","pov.uc.n","pov.c.nm","pov.c.n")],
                     by=list(a.tmp$area),FUN=mean)
  a.tmp$count <- a.count$x
  names(a.tmp)[1] <- "area"
  ##  names(a.tmp)[dim(a.tmp)[2]] <- "count"
  agg.file <- rbind(agg.file,a.tmp)
  cat("teller:  ",teller <- teller+1,"\n")
}
#  
close(cen.data)  
#
## finalize agg.file for double area entries
#

for(ar in unique(agg.file$area)){
  arpart <- agg.file[agg.file$area==ar,]
  rn <- row.names(arpart)
  if(dim(arpart)[1]>1){
    povests <- as.numeric(arpart$count %*% as.matrix(arpart[,2:6]))/sum(arpart$count)
    agg.file[rn[1],] <- c(ar,povests,sum(arpart$count))
    agg.file[rn[-1],] <- c(-1,0,0,0,0,0,0)
  }
}
agg.file <- agg.file[agg.file$area != -1,]

pdf(file="plot5.pdf")
inc.pov <- order(agg.file$pov.true)
plot(agg.file$pov.true,agg.file$pov.c.nm,cex=0.3)
plot(agg.file$pov.true[inc.pov],type="l")
lines(agg.file$pov.c.nm[inc.pov],col="red")
lines(agg.file$pov.c.n[inc.pov],col="blue")
dev.off()

plot(pov.true$true.pov[inc.pov],type="l")
lines(pov.pred$pred.pov[inc.pov],col="red")
lines(pov.pred.normal$pred.pov.normal[inc.pov],col="blue")



cat("MSE MN prediction at area level:  ",msdiff(agg.file$pov.true,agg.file$pov.c.nm),"\n")
cat("MSE N prediction at area level:   ",msdiff(agg.file$pov.true,agg.file$pov.c.n),"\n")


