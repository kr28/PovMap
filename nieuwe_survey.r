#### make new survey #### 
library(foreign)

cen.data <- read.csv("../Shared with Chris/sim_census_ah.csv",header=T)
cen.data.u <- read.csv("../Shared with Chris/sim_census_a.csv",header=T)
u.expand <- rep(cen.data.u$u.a,each=3000)
eps.expand <- cen.data$y.ah-cen.data$x.ah-u.expand
x.new <- sort(cen.data$x.ah)
y.new <- x.new+u.expand-eps.expand ## notice minus in front of eps
cen.data.new=data.frame(X=cen.data$X,area.id=cen.data$area.id,y.ah=y.new,x.ah=x.new,pov.ah=ifelse(y.new < -0.5,1,0))

set.seed(112)
my.samp <- numeric(0)
for(area in 0:499){
  hh.sel <- sample(1:3000,15)
  my.samp <- c(my.samp,area*3000+hh.sel)
}

sur.data.new <- data.frame(X=1:7500,area.id=rep(1:500,each=15),y.ah=y.new[my.samp],x.ah=x.new[my.samp])

write.table(cen.data.new,file="/Users/Chris/Documents/PovMap/sim2_census_ah.csv",row.names=F,sep=",")

write.table(sur.data.new,file="/Users/Chris/Documents/PovMap/sim2_survey_ah.csv",row.names=F,sep=",")

## x not sorted, rest same:

x.new <- (cen.data$x.ah)
y.new <- x.new+u.expand-eps.expand ## notice minus in front of eps
cen.data.new=data.frame(X=cen.data$X,area.id=cen.data$area.id,y.ah=y.new,x.ah=x.new,pov.ah=ifelse(y.new < -0.5,1,0))

set.seed(112)
my.samp <- numeric(0)
for(area in 0:499){
  hh.sel <- sample(1:3000,15)
  my.samp <- c(my.samp,area*3000+hh.sel)
}

sur.data.new <- data.frame(X=1:7500,area.id=rep(1:500,each=15),y.ah=y.new[my.samp],x.ah=x.new[my.samp])

write.table(cen.data.new,file="/Users/Chris/Documents/PovMap/sim3_census_ah.csv",row.names=F,sep=",")

write.table(sur.data.new,file="/Users/Chris/Documents/PovMap/sim3_survey_ah.csv",row.names=F,sep=",")
