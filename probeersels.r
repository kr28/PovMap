# mix model met twee populaties

# superpopulation paramters

mu1 <- 3
mu2 <- -1
sig1 <- 0.5
sig2 <- 0.7
npop <- 100000 # 100000
sh2 <- 0.02 # share group 2
area <- (rep(1:100,each=1000))
sigu <- sqrt(0.01*sig1^2)


# generate data

# population
group2 <- rbinom(npop,1,sh2)
group2 <- 0

# data
data.1 <- group2 * rnorm(npop,mu2,sig2)+(1-group2)*rnorm(npop,mu1,sig1)
area.eff <- rep(rnorm(100,m=0,s=sigu),each=1000)
data <- data.1+area.eff
qq <- qqnorm(data)
#hist(data,n=20)


# sample
# samp <- sample(data,300,repl=T) # case without areas
samp.ar <- as.numeric(sapply(1:100,function(a) sample(data[area==a],15)))

# strategy 1:

                                        # estimate model by pooling data, leading to muhat and sighat
reg1 <- lm(samp~1)
muhat <- reg1$coef
sighat <- sqrt(vcov(reg1)[1,1])

                                        # simulate data by picking muhat and drawing from residuals
mu.sim <-rnorm(1,muhat,sighat)
eps.sim <- sample(reg1$resid,npop,repl=T)
data.sim.1 <- mu.sim+(eps.sim) # consider jitter

                                        # do qq plot comparing with true data
qq <- qqplot(data,data.sim.1)

# with groups

samp.gr <- as.factor(rep(1:100,each=15))

sigu <- sqrt(0.01*sig1^2)
area.eff <- rep(rnorm(100,m=0,s=sigu),each=1000)
data <- data.1+area.eff
samp.ar <- as.numeric(sapply(1:100,function(a) sample(data[area==a],15)))
reg.lme <- lme(samp.ar ~ 1, random= ~1 |  samp.gr)

data.sim.1 <- mu.sim+(eps.sim) # consider jitter

                                        # do qq plot comparing with true data
qq <- qqplot(data,data.sim.1)

# with groups

samp.gr <- as.factor(rep(1:100,each=15))


estimates <- numeric(100)

sigu <- sqrt(0.001*sig1^2)

for(i in 1:100) { 
area.eff <- rep(rnorm(100,m=0,s=sigu),each=1000)
data <- data.1+area.eff
samp.ar <- as.numeric(sapply(1:100,function(a) sample(data[area==a],15)))
reg.lme <- lme(samp.ar ~ 1, random= ~1 |  samp.gr)
rrr <- reg.lme$apVar
if(length(rrr)==1) estimates[i] <- 0 else
estimates[i] <- as.numeric(exp(attr(reg.lme$apVar,"Pars")[1]))
}

table(estimates==0)
sqrt(mean(estimates^2))
sigu

# RMS estimated sigu is 0.062, with actual sigu 0.05
# RMS estimated sigu is 0.0365, with actual sigu 0.0354
# RMS estimated sigu is 0.050, with actual sigu 0.0158




# strategy 2: robust center: median

muhat.2 <- median(samp)
eps.sim.2 <- samp - muhat.2
data.sim.2 <- muhat.2+sample(eps.sim.2,npop,repl=T)
qq <- qqplot(data,data.sim.2)

# strategy 3: leave out outliers

reg2 <- lm(samp~1)
bnds <- quantile(reg2$resid,prob=c(0.025,0.975))
keep <- reg2$resid>bnds[1] & reg2$resid<bnds[2]
reg2.a <- lm(samp[keep]~1)

muhat <- reg2.a$coef
sighat <- sqrt(vcov(reg2.a)[1,1])

                                        # simulate data by picking muhat and drawing from residuals
mu.sim <-rnorm(1,muhat,sighat)
eps.sim <- sample(samp-mu.sim,npop,repl=T)
data.sim.3 <- mu.sim+(eps.sim) # consider jitter

qq <- qqplot(data,data.sim.1)
qq <- qqplot(data,data.sim.2)
qq <- qqplot(data,data.sim.3)

qqnorm(reg2.a$res)
qqnorm(reg1$res)

sd(data.sim.1)
sd(data.sim.2)
sd(data.sim.3)

