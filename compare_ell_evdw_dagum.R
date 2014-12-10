library(VGAM)



#setwd("c:/WB/NM")

#set nested error variance parameters
rho <- 0.3
sigma2.e <- 1

sigma2.u <- rho*sigma2.e
sigma2.eps <- sigma2.e - sigma2.u

sigma2.u
sigma2.eps

#number of areas
A <- 5000

#number of households sampled per area
n.a <- 5

#total sample size
n <- A*n.a

#sample the household errors: eps.ah
shape.eps <- 5
rdag <- rdagum(n, 2, 1, shape.eps)
eps.ah <- log(rdag)
eps.ah <- eps.ah - mean(eps.ah)
eps.ah <- ((sigma2.eps/var(eps.ah))^0.5)*eps.ah
mean(eps.ah)
var(eps.ah)
kdens.eps <- density(eps.ah)
#plot(kdens.eps)

#sample the area errors: u.a
set.seed(123)
shape.u <- 0.2
rdag <- rdagum(A, 2, 1, shape.u)
u.a <- log(rdag)
u.a <- u.a - mean(u.a)
u.a <- ((sigma2.u/var(u.a))^0.5)*u.a
mean(u.a)
var(u.a)
kdens.u <- density(u.a)

#create area identifier
area.id <- rep(1:A,each=n.a)
area.id.a <- rep(1:A,each=1)

#construct the total household error
e.ah <- rep(u.a,each=n.a) + eps.ah
kdens.e <- density(e.ah)
#plot(kdens.e)
#
# implement non-parametric approach by ELL
#

#let us assume that the variance parameters are known
e.a <- aggregate(e.ah,by=list(area.id),FUN=mean)$x

#re-scale e.a to obtain estimate of u.a
u.a.hat <- (e.a - mean(e.a))*sqrt(sigma2.u)/sd(e.a)

#similarly, obtain an estimate for eps.ah
eps.ah.hat <- e.ah - rep(u.a.hat,each=n.a)
eps.ah.hat <- (eps.ah.hat - mean(eps.ah.hat))*sqrt(sigma2.eps)/sd(eps.ah.hat)

#compare histograms of u.a.hat and eps.ah.hat to "true" densities
hist(u.a.hat, breaks=50, freq=F)
lines(kdens.u, add=T, col="Red")

#estimate for eps.ah barely affected
hist(eps.ah.hat, breaks=50, freq=F)
lines(kdens.eps, add=T, col="Red")







