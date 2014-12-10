### to get quick model ###
library(foreign)

b.data <- read.csv("~/Documents/Povmap/Brazil/mg_survey_2000_nonzero.csv",header=T)

summary(reg.beta <- lm(lny ~single +female_kids+ age07_1+ age07_2plus+ age714_1+ age714_2plus+ lnage+ edu_1+ edu_2+ edu_3+ edu_4+ edu_5,data=b.data)) 
## summary(reg.beta <- lm(lny~ single+ female_kids+ age07_1+ age07_2plus+ age714_1+ age714_2plus+ lnage+ edu_1 +edu_2 +edu_3 +edu_4 +edu_5 +dw_pavimen +dw_sewage_1 +asset_auto +asset_pc+ asset_tv+ asset_microwave+ asset_washing+ asset_radio+ asset_video,data=b.data))

library(nlme)
reg.lme <- lme(formula(reg.beta),random=~1|mun,data=b.data)
(v.eps.lme <- as.numeric(exp(attr(reg.lme$apVar,"Pars"))[2]^2))
(v.u.lme <- as.numeric(exp(attr(reg.lme$apVar,"Pars"))[1]^2))

### FE step ###
###
library(psych) # for matrix functions
library(MASS)

## save.image(file="~/Documents/Povmap/Brazil/temp.RData")

sur.data <- b.data
rm(b.data)
sur.data$cluster <- sur.data$mun
clusterSize<-as.vector(table(sur.data$cluster)) ## assumed here: equal cluster size
r_hat <- resid(reg.beta)
n_hlds <- length(r_hat)


X.ols <- model.matrix(reg.beta)
## take out mu_logpop (is constant within municipio)
##
X.ols <- X.ols[,-29]
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
y <- sur.data$lny
Dka <- diag(ka)
ybar <- aggregate(y,by=list(sur.data$cluster),FUN=mean)$x
a.X <- aggregate(X,by=list(sur.data$cluster),FUN=mean)
Xbar <- as.matrix(a.X[,2:dim(a.X)[2]])
XTX <- (t(X)%*%X)-(t(Xbar)%*%Dka%*%Xbar)
XTy <- (t(X)%*%y)-(t(Xbar)%*%Dka%*%ybar)
beta.within <- ginv(XTX)%*%XTy
eps_hat <- (y-rep(ybar,ka))-(X%*%beta.within-rep(Xbar%*%beta.within,ka))
sigma2.eps <- sum(eps_hat*eps_hat)/(n_hlds-(d-1)-n_area)
### compare with lme
##
## > sigma2.eps
## [1] 0.7871393
##

v.eps <- sigma2.eps / ka  ## warning: v.eps has different meaning further down!

# H3 variance of v: Corollary 3
a.X.ols <- aggregate(X.ols,by=list(sur.data$cluster),FUN=mean)
UX.ols <- as.matrix(a.X.ols[,2:dim(a.X.ols)[2]])
UX.ols <- Dka%*%UX.ols
trupxu <- tr(ginv(t(X.ols)%*%X.ols)%*%t(UX.ols)%*%UX.ols)
sigma2.v <- (n_hlds-d.ols)*((sum(r_hat*r_hat)/(n_hlds-d.ols))-sigma2.eps)/(n_hlds-trupxu)
## > sigma2.v
## [1] 0.003695134  ## very low compared to lme estimate


#
# obtain GLS estimate of beta
#

# reset X.ols
#
X.ols <- model.matrix(reg.beta)

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

## save.image(file="~/Documents/Povmap/Brazil/temp.RData")

ydata <- data.frame(area=sur.data$cluster,hh=rownames(sur.data),y=res.gls)
## ydata <- cbind(ydata,model.matrix(reg.beta))
## ydata$"(Intercept)" <- NULL

#### FUNCTIONS USED

library(deconvo,lib.loc="deconvo")

## save.image(file="~/Documents/Povmap/Brazil/temp.RData")

## clean up
rm(list=setdiff(ls(),c("ydata","sigma2.eps","sigma2.v","reg.beta","beta.gls","sur.data"))) ## now only ydata, sigma2.eps etc. remain

detach(package:foreign)
detach(package:psych)
detach(package:plm)
detach(package:MASS)
detach(package:sandwich)
detach(package:zoo)
detach(package:Formula)
detach(package:bdsmatrix)
detach(package:nlme)

library(splines)

#
# ---------------------------- Fit the mixture-distribution for F_v ----------------------------

## replace by new procedure
##
## ebar moet gemiddelde 0 hebben, en variantie volgens H3
##

ybardata <- aggregate(ydata[,c("y")],by=list(ydata$area),FUN=mean)
names(ybardata) <- c("area","y")

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
ncomp.uest <- 3
# ncomp.uest <- 2
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
(outcomes <- run.prog.udist(5))  ## 50 is the number of iterations
ppi.new <- outcomes$ppi.new
mu.new <- outcomes$mu.new
sig.new <- outcomes$sig.new
sig2.new <- outcomes$sig2.new

dens.uhat <- Vectorize(function(u) as.numeric(ppi.new%*% dnorm(u,mu.new,sig.new)))
dens.uhat.norm <- Vectorize(function(u) as.numeric(dnorm(u,0,sqrt(sigma2.v))))

## check expectation
c(mean.uhat=ppi.new%*%mu.new,var.uhat=ppi.new%*%(mu.new^2+sig2.new))

## check densities

curve(dens.uhat,min(ybar),max(ybar),col="red")
curve(dens.uhat.norm,min(ybar),max(ybar),col="blue",add=T)

## curve(dens.u,-3,3,add=TRUE,col="red")
lines(density(ybar,n=200,from=min(ybar),to=max(ybar))) ## since i don't have the true density

pdf(file="plotzo1.pdf")
curve(dens.uhat,min(ybar),max(ybar),col="red")
curve(dens.uhat.norm,min(ybar),max(ybar),col="blue",add=T)
lines(density(ybar,n=200,from=min(ybar),to=max(ybar))) ## since i don't have the
dev.off()
## ## save these on pdf
## pdf(file="plot1.pdf")
## curve(dens.uhat,min(ybar),max(ybar),col="red")
## lines(density(ybar,n=200,from=min(ybar),to=max(ybar))) ## since i don't have the true density
## dev.off()

## clean up

## save.image(file="~/Documents/Povmap/Brazil/temp.RData")

rm(list=c(setdiff(ls(),c("mu.new","ppi.new","sig.new","sig2.new","v.eps","v.u","ybardata","ydata","sigma2.eps","sigma2.v","reg.beta","beta.gls","sur.data","sur.data.u"))))

################################### EPS DISTRIBUTION #####################

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
## ncomp.epsest <- 3

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
(outcomes <- run.prog.epsdist(5,-2.4,2.4)) ## again, use more iterations if you want

lam.new <- as.numeric(outcomes$lam.new)
nu.new <- as.numeric(outcomes$nu.new)
tau.new <- as.numeric(outcomes$tau.new)
tau2.new <- as.numeric(outcomes$tau2.new)

dens.epshat <- Vectorize(function(t) as.numeric(lam.new%*% dnorm(t,nu.new,tau.new)))
dens.epshat.norm <- Vectorize(function(t) as.numeric(dnorm(t,0,sqrt(sigma2.eps))))

## check expectation
c(mean.epshat=lam.new%*%nu.new,var.epshat=lam.new%*%(nu.new^2+tau2.new))

## check densities
curve(dens.epshat,min(ydata$y),max(ydata$y),col="red")
curve(dens.epshat.norm,min(ydata$y),max(ydata$y),col="blue",add=T)
lines(density(ydata$y,n=200,from=min(ydata$y),to=max(ydata$y))) ## since i don't have the true density

pdf(file="plot2ma.pdf")
curve(dens.epshat,min(ydata$y),max(ydata$y),col="red")
curve(dens.epshat.norm,min(ydata$y),max(ydata$y),col="blue",add=T)
lines(density(ydata$y,n=200,from=min(ydata$y),to=max(ydata$y))) ## since i don't have the true density
dev.off()

## ## on pdf:
## pdf(file="plot2.pdf")
## curve(dens.epshat,min(ydata$eps),max(ydata$eps),col="red")
## lines(density(ydata$eps,n=200,from=min(ydata$eps),to=max(ydata$eps))) ## since i don't have the true density
## dev.off()


##############
## prepare for census imputation
##
## you might want to remove some variables:
##

## PARAMETER
pov.line <- 4.601
X.ols <- model.matrix(reg.beta)
X.vars <- colnames(X.ols)[-1]

################ CENSUS IMPUTATION ############
###
### Method 1: read in all data (if you have enough memory)
###
## cen.data <- read.csv("../Shared with Chris/sim_census_ah.csv",header=T)
## cen.data.u <- read.csv("../Shared with Chris/sim_census_a.csv",header=T)
cen.data <- read.csv("/Users/Chris/Documents/PovMap/Brazil/mg_2000_nonzero.csv",header=T)

## remove zeroes
## cen.data <- cen.data[cen.data$ypc > 0,]
nhh.cen <- dim(cen.data)[1]
X.vars <- colnames(X.ols)[-1]

X.mat <- as.matrix(cen.data[,X.vars])
yhat <- beta.gls[1]+X.mat%*%beta.gls[-1]

## individual, unconditional

estpov.ind <- sapply(1:nhh.cen, est.pov.ind)  #NM
estpov.ind.normal <- sapply(1:nhh.cen,est.pov.ind.normal) #N
poor.index <- ifelse(cen.data$lny < pov.line,1,0)

## MSE, individual, unconditional

cat("MSE normal mixture:    ",mean((estpov.ind-poor.index)^2),"\n")
cat("MSE normal:            ",mean((estpov.ind.normal-poor.index)^2),"\n")
cat("bias normal mixture:   ",mean((estpov.ind-poor.index)),"\n")
cat("bias normal:           ",mean((estpov.ind.normal-poor.index)),"\n")

plot(estpov.ind,estpov.ind.normal,cex=0.1)

## aggregate, unconditional

poverty.agg <- aggregate(poor.index,by=list(cen.data$mun),FUN=mean)$x

estpov.agg <- aggregate(estpov.ind,by=list(cen.data$mun),FUN=mean)$x
estpov.agg.normal <- aggregate(estpov.ind.normal,by=list(cen.data$mun),FUN=mean)$x

plot(estpov.agg,poverty.agg,cex=0.1)
plot(estpov.agg,estpov.agg.normal,cex=0.1)
plot(poverty.agg[order(poverty.agg)],type="l")
lines(estpov.agg.normal[order(poverty.agg)],cex=0.1,col="blue")
lines(estpov.agg[order(poverty.agg)],col="red")

pdf(file="plot6ma.pdf")
plot(poverty.agg[order(poverty.agg)],type="l")
lines(estpov.agg.normal[order(poverty.agg)],cex=0.1,col="blue")
lines(estpov.agg[order(poverty.agg)],col="red")
dev.off()

cat(" root mean square difference NM and N:  ",sqrt(mean((estpov.agg-estpov.agg.normal)^2)),"\n")
cat(" mean difference NM and N:  ",(mean((estpov.agg-estpov.agg.normal))),"\n")

cat("MSE normal mixture (aggregated, unconditional):    ",mean((estpov.agg-poverty.agg)^2),"\n")
cat("MSE normal (aggregated, unconditional):            ",mean((estpov.agg.normal-poverty.agg)^2),"\n")
cat("Bias normal mixture (aggregated, unconditional):   ",mean((estpov.agg-poverty.agg)),"\n")
cat("Bias normal (aggregated, unconditional):           ",mean((estpov.agg.normal-poverty.agg)),"\n")


## individual, conditional

cestpov.ind <- sapply(1:nhh.cen,function(n) c.est.pov.ind(n,cen.data$mun[n]))
cestpov.ind.normal <- sapply(1:nhh.cen ,function(n) c.est.pov.ind.normal(n,cen.data$mun[n]))

cat("MSE normal mixture (conditional):    ",mean((cestpov.ind-poor.index)^2),"\n")
cat("MSE normal (conditional):            ",mean((cestpov.ind.normal-poor.index)^2),"\n")
cat("Bias normal mixture (conditional):   ",mean((cestpov.ind-poor.index)),"\n")
cat("Bias normal (conditional):           ",mean((cestpov.ind.normal-poor.index)),"\n")
## 
plot(cestpov.ind,cestpov.ind.normal,cex=0.1)

## aggregate, conditional

cestpov.agg <- aggregate(cestpov.ind,by=list(cen.data$mun),FUN=mean)$x
cestpov.agg.normal <- aggregate(cestpov.ind.normal,by=list(cen.data$mun),FUN=mean)$x

plot(cestpov.agg,cestpov.agg.normal,cex=0.1)

plot(poverty.agg[order(poverty.agg)],type="l")
lines(cestpov.agg.normal[order(poverty.agg)],cex=0.1,col="blue")
lines(cestpov.agg[order(poverty.agg)],col="red")
pdf("plot3ma.pdf")
plot(poverty.agg[order(poverty.agg)],type="l")
lines(cestpov.agg.normal[order(poverty.agg)],cex=0.1,col="blue")
lines(cestpov.agg[order(poverty.agg)],col="red")
dev.off()


cat(" root mean square difference NM and N:  ",sqrt(mean((cestpov.agg-cestpov.agg.normal)^2)),"\n")
cat(" mean difference NM and N:  ",(mean((cestpov.agg-cestpov.agg.normal))),"\n")

cat("MSE normal mixture (aggregated, conditional):    ",mean((cestpov.agg-poverty.agg)^2),"\n")
cat("MSE normal (aggregated, conditional):            ",mean((cestpov.agg.normal-poverty.agg)^2),"\n")
cat("Bias normal mixture (aggregated, conditional):   ",mean((cestpov.agg-poverty.agg)),"\n")
cat("Bias normal (aggregated, conditional):           ",mean((cestpov.agg.normal-poverty.agg)),"\n")
## 

inc.pov <- order(poverty.agg)

t.sm <- loess(cestpov.agg[inc.pov] ~ I(1:length(cestpov.agg) ))
t2.sm <- loess(cestpov.agg.normal[inc.pov] ~ I(1:length(cestpov.agg)))
plot(poverty.agg[inc.pov],type="l")
lines(t.sm$fitted,col="red")
lines(t2.sm$fitted,col="blue")

pdf("plot4ma.pdf")
plot(poverty.agg[inc.pov],type="l")
lines(t.sm$fitted,col="red")
lines(t2.sm$fitted,col="blue")
dev.off()

rm(cen.data)
rm(cen.data.u)
