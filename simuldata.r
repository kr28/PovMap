rm(list=ls())
library(foreign)
library(stringr)

# 2. 300 clusters

nclust <- 300
nhhpc <- 10 # max 99

hid <- rep(1:10,nclust)
# hid <- str_pad(hid,width=2,pad="0")
cid <- rep(1:nclust,each=nhhpc)
# cid <- str_pad(cid,width=3,pad="0")
# fid <- paste(cid,hid,sep="")
fid <- cid*100+hid

set.seed(123)
xvec <- rnorm(nclust*nhhpc,m=2,s=sqrt(1))
eps <- rnorm(nclust*nhhpc,m=0,s=sqrt(0.93))
eta1 <-rnorm(nclust,m=0.32,s=sqrt(0.02))
eta2 <- rnorm(nclust,m=-0.24,s=sqrt(0.02))
eta0 <- rnorm(nclust,m=0,s=0.1)
pieta1 <- 0.3/0.7-0.1
pieta2 <- 0.4/0.7-0.1
pieta0 <- 0.2
eta.comp <- rmultinom(nclust,1,c(pieta1,pieta0,pieta2))
eta <- rep(eta1*eta.comp+eta2*(1-eta.comp),each=nhhpc)
y <- 3.0+xvec+eta+eps
# pieta1*0.4+pieta2*-0.3
# var(eta)
# var(eta)/(var(eta)+var(eps))

data300 <- data.frame(id=fid,cid=cid,hid=hid,y=y,x=xvec,eta=eta,eps=eps)
data300[1:24,]
# hist(eta,prob=T)
# plot(function(t) dnorm(t,m=0,s=sd(eta)),-1,1,add=T)
attr(data300,"var.eps") <- 0.93
attr(data300,"var.eta") <- c(0.05,0.05)
attr(data300,"mu.eta") <- c(0.32,-0.24)
attr(data300,"pi.eta") <- c(3.0/7.0,4.0/7.0)

write.dta(data300,file="data300.dta",version=10,convert.factors="string")

# 1. 100 clusters

nclust <- 100
nhhpc <- 10 # max 99

hid <- rep(1:10,nclust)
# hid <- str_pad(hid,width=2,pad="0")
cid <- rep(1:nclust,each=nhhpc)
# cid <- str_pad(cid,width=3,pad="0")
# fid <- paste(cid,hid,sep="")
fid <- cid*100+hid

set.seed(100)
xvec <- rnorm(nclust*nhhpc,m=2,s=sqrt(1))
eps <- rnorm(nclust*nhhpc,m=0,s=sqrt(0.93))
eta1 <-rnorm(nclust,m=0.32,s=sqrt(0.05))
eta2 <- rnorm(nclust,m=-0.24,s=sqrt(0.05))
pieta1 <- 0.3/0.7
pieta2 <- 0.4/0.7
eta.comp <- rbinom(nclust,1,pieta1)
eta <- rep(eta1*eta.comp+eta2*(1-eta.comp),each=nhhpc)
y <- 3.0+xvec+eta+eps
# pieta1*0.4+pieta2*-0.3
# var(eta)
# var(eta)/(var(eta)+var(eps))

data100 <- data.frame(id=fid,cid=cid,hid=hid,y=y,x=xvec,eta=eta,eps=eps)
data100[1:24,]
# hist(eta,prob=T)
# plot(function(t) dnorm(t,m=0,s=sd(eta)),-1,1,add=T)
## attr(data300,"var.eps") <- 0.93
## attr(data300,"var.eta") <- c(0.05,0.05)
## attr(data300,"mu.eta") <- c(0.32,-0.24)
## attr(data300,"pi.eta") <- c(3.0/7.0,4.0/7.0)

write.dta(data100,file="data100.dta",version=10,convert.factors="string")

# 3. 900 clusters

nclust <- 900
nhhpc <- 10 # max 99

hid <- rep(1:10,nclust)
# hid <- str_pad(hid,width=2,pad="0")
cid <- rep(1:nclust,each=nhhpc)
# cid <- str_pad(cid,width=3,pad="0")
# fid <- paste(cid,hid,sep="")
fid <- cid*100+hid

set.seed(900)
xvec <- rnorm(nclust*nhhpc,m=2,s=sqrt(1))
eps <- rnorm(nclust*nhhpc,m=0,s=sqrt(0.93))
eta1 <-rnorm(nclust,m=0.32,s=sqrt(0.05))
eta2 <- rnorm(nclust,m=-0.24,s=sqrt(0.05))
pieta1 <- 0.3/0.7
pieta2 <- 0.4/0.7
eta.comp <- rbinom(nclust,1,pieta1)
eta <- rep(eta1*eta.comp+eta2*(1-eta.comp),each=nhhpc)
y <- 3.0+xvec+eta+eps
# pieta1*0.4+pieta2*-0.3
# var(eta)
# var(eta)/(var(eta)+var(eps))

data900 <- data.frame(id=fid,cid=cid,hid=hid,y=y,x=xvec,eta=eta,eps=eps)
data900[1:24,]
# hist(eta,prob=T)
# plot(function(t) dnorm(t,m=0,s=sd(eta)),-1,1,add=T)
## attr(data300,"var.eps") <- 0.93
## attr(data300,"var.eta") <- c(0.05,0.05)
## attr(data300,"mu.eta") <- c(0.32,-0.24)
## attr(data300,"pi.eta") <- c(3.0/7.0,4.0/7.0)

write.dta(data900,file="data900.dta",version=10,convert.factors="string")


# 4. 1000 clusters

nclust <- 1000
nhhpc <- 30 # max 99

hid <- rep(1:nhhpc,nclust)
# hid <- str_pad(hid,width=2,pad="0")
cid <- rep(1:nclust,each=nhhpc)
# cid <- str_pad(cid,width=3,pad="0")
# fid <- paste(cid,hid,sep="")
fid <- cid*100+hid

set.seed(1000)
xvec <- rnorm(nclust*nhhpc,m=2,s=sqrt(1))
# eps <- rnorm(nclust*nhhpc,m=0,s=sqrt(0.93))
eta1 <-rnorm(nclust,m=-0.22,s=sqrt(0.04))
eta2 <- rnorm(nclust,m=0.44,s=sqrt(0.04))
eta0 <- rnorm(nclust,m=0,s=sqrt(0.04))
pieta1 <- 0.4
pieta2 <- 0.2
pieta0 <- 0.4
eta.comp <- rmultinom(nclust,1,c(pieta1,pieta0,pieta2))
eta <- rep(eta.comp[1,]*eta1+eta.comp[2,]*eta0+eta.comp[3,]*eta2,each=nhhpc)

eps1 <-rnorm(nclust*nhhpc,m=-1.155,s=sqrt(0.45))
eps2 <- rnorm(nclust*nhhpc,m=0.77,s=sqrt(0.45))
eps0 <- rnorm(nclust*nhhpc,m=0,s=sqrt(0.45))
peps1 <- 0.2
peps2 <- 0.3
peps0 <- 0.5
peps.comp <- rmultinom(nclust*nhhpc,1,c(peps1,peps0,peps2))
eps <- peps.comp[1,]*eps1+peps.comp[2,]*eps0+peps.comp[3,]*eps2

y <- 3.0+xvec+eta+eps
# pieta1*0.4+pieta2*-0.3
# var(eta)
# var(eta)/(var(eta)+var(eps))

data1000 <- data.frame(id=fid,cid=cid,hid=hid,y=y,x=xvec,eta=eta,eps=eps)
data1000[1:40,]

write.dta(data1000,file="data1000.dta",version=10,convert.factors="string")
