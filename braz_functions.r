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
        pky[k,i] <- integrate(fun,-1.9,1.9)$value
        fun1 <- function(t) {
          r <- (y.sample[i]-t) * taufun(t,y.sample[i])[k]
          if(is.nan(r)) return(0) else return(r)
        }
        fun <- Vectorize(function(t) fun1(t)*fun2(t))
        mky[k,i] <- integrate(fun,-1.9,1.9)$value
        fun1 <- function(t) {
          r <- (y.sample[i]-t-mu.old[k])^2* taufun(t,y.sample[i])[k]
          if(is.nan(r)) return(0) else return(r)
        }
        fun <- Vectorize(function(t) fun1(t)*fun2(t))
        s2ky[k,i] <- integrate(fun,-1.9,1.9)$value    
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


run.prog.epsdist <- function(niter=21,from=-5,to=5){
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
          pkyy[k,i2,i1] <- integrate(fun,from,to)$value
          fun1 <- function(t) {
            r <- (y.sample[i2]-t) * taufun.eps(t,y.sample[i2])[k]
            if(is.nan(r)) return(0) else return(r)
          }
          fun <- Vectorize(function(t) fun1(t)*fun2(t))
          mkyy[k,i2,i1] <- integrate(fun,from,to)$value
          fun1 <- function(t) {
            r <- (y.sample[i2]-t-nu.old[k])^2* taufun.eps(t,y.sample[i2])[k]
            if(is.nan(r)) return(0) else return(r)
          }
          fun <- Vectorize(function(t) fun1(t)*fun2(t))
          s2kyy[k,i2,i1] <- integrate(fun,from,to)$value    
        }
      }
    }
## plot(y.sample,pky[5,])
  
    s2ky.p <- mky.p <- pky.p <- matrix(0,nrow=length(tau.old),ncol=length(ybar))
#    for(a in unique(ydata$area)){
    for(a in (1:length(ybar))){
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

est.pov.ind <- function(hh,z=pov.line,ppi=ppi.new,mu=mu.new,sig=sig.new,
                    lam=lam.new,nu=nu.new,tau=tau.new,beta=beta.gls){
  ## x does not include the intercept
  yhat.loc <- yhat[hh]
  ppilam <- as.numeric(ppi %o% lam)
  munu <- as.numeric(outer(mu,nu,FUN="+"))
  sigtau <- as.numeric(sqrt(outer(sig^2,tau^2,"+")))
  return(as.numeric(ppilam %*% pnorm(z-yhat.loc,munu,sigtau)))
}

est.pov.ind.normal <- function(hh,z=pov.line,sig2.eps=v.eps,sig2.u=v.u,beta=beta.gls){
  ## x does not include the intercept
  yhat.loc <- yhat[hh]
  return(pnorm(z-yhat.loc,0,sqrt(sig2.eps+sig2.u)))
}

c.est.pov.ind <- function(hh,ar,z=pov.line,ppi=ppi.new,mu=mu.new,sig=sig.new,
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
  return(est.pov.ind(hh,pov.line,omeg,cmu.u,csig.u,lam,nu,tau,beta))
##  return(list(omeg=omeg,cmu=cmu.u,csig=csig.u,csig2=csig2.u,exp.pov=est.pov.ind(x,pov.line,omeg,cmu.u,csig.u,lam,nu,tau,beta)))
}

c.est.pov.ind.normal <- function(hh,ar,beta=beta.gls){
  ## x does not include the constant term
  ## ar is area indicator
  return(c.est.pov.ind(hh,ar,pov.line,1,0,sqrt(v.u),1,0,sqrt(v.eps),beta))
}


msdiff <- function(v1,v2){ ## mean square difference
  mean((v1-v2)^2)
}
