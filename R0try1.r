###############################################
##
###
#### And Now for R0
###
##
################################################
###using data from CopyOfExplorationR02.r

for (i in 1:206){
  
  data$ro[i]<-(data$k[i] * ( (1/(1-(data$prevalence[i]/data$N[i]))) ^ (1/data$k[i]) ) ) - data$k[i]
}
plot(data$ro,data$T20,ylim=c(0,1),xlim=c(0,2),xaxt="n",
     ylab="Proportion of parasites in top most infected hosts",
     xlab=expression(paste(R[0])))
axis(1,at=seq(0,2,0.2),labels=seq(0,2,0.2),par(las=1))
log.binom<-function(p.vec){
  
  a<-p.vec[1]
  b<-p.vec[2]
  
  pred1a<- ((exp(a + b * data$ro)) / (1 + exp(a + b * data$ro)) ) 
  prev1<-data$T20
  
  loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
  -sum(loglik1a,  na.rm=T)
}
n.param<-2
logmod<-optim(c(0,0),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,10))
logmod
nc<-seq(0,2,0.01)
pred2<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) )
lines(nc,pred2,lwd=2,lty=2,col="black")

log.binom<-function(p.vec){
  
  a<-p.vec[1]
  b<-p.vec[2]
  
  pred1a<- ((exp(a + b * data$ro)) / (1 + exp(a + b * data$ro)) ) 
  prev1<-data$T60
  
  loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
  -sum(loglik1a,  na.rm=T)
}
n.param<-2
logmod<-optim(c(0,0),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,10))
logmod
nc<-seq(0,2,0.01)
pred2<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) )
lines(nc,pred2,lwd=1,lty=2,col="grey")

##What happen to R0 when the popualtion is treated?
k01<-numeric(206)
for (i in 1:ncol(distrib01)){
  sampmean<-mean(distrib01[,i])
  var<-var(distrib01[,i])
  k01[i]<-(sampmean^2-(var/nrow(distrib01)))/(var-sampmean)
}

for (i in 1:ncol(distrib10)){
  sampmean<-mean(distrib10[,i])
  var<-var(distrib10[,i])
  data$k10[i]<-(sampmean^2-(var/nrow(distrib10)))/(var-sampmean)
}
data$ro10[155]<-data$k10[155]<-NA

data$k01<-k01;data$k05<-k05;data$k10<-k10;data$k15<-k15;data$k20<-k20;
data$k25<-k25;data$k30<-k30;data$k35<-k35;data$k40<-k40;data$k45<-k45
data$k30[is.infinite(data$k30)] <- NA
data$k35[is.infinite(data$k35)] <- NA
data$k40[is.infinite(data$k40)] <- NA
data$k45[is.infinite(data$k45)] <- NA
data$k50[is.infinite(data$k50)] <- NA
data$k55[is.infinite(data$k55)] <- NA
data$k60[is.infinite(data$k60)] <- NA
data$k65[is.infinite(data$k65)] <- NA
data$k70[is.infinite(data$k70)] <- NA
data$k75[is.infinite(data$k75)] <- NA
data$k80[is.infinite(data$k80)] <- NA
data$k85[is.infinite(data$k85)] <- NA

for (i in 1:206){
  data$ro[i]<-(data$k[i] * ( (1/(1-(data$prevalence[i]/data$N[i]))) ^ (1/data$k[i]) ) ) - data$k[i]
  
  data$ro01[i]<-(data$k01[i] * ( (1/(1-(sum(prev01[i])/21))) ^ (1/data$k01[i]) ) )  - data$k01[i]
  data$ro05[i]<-(data$k05[i] * ( (1/(1-(sum(prev05[i])/21))) ^ (1/data$k05[i]) ) ) - data$k05[i]
  data$ro10[i]<-(data$k10[i] * ( (1/(1-(sum(prev10[i])/21))) ^ (1/data$k10[i]) ) ) - data$k10[i]
  data$ro15[i]<-(data$k15[i] * ( (1/(1-(sum(prev15[i])/21))) ^ (1/data$k15[i]) ) ) - data$k15[i]
  data$ro20[i]<-(data$k20[i] * ( (1/(1-(sum(prev20[i])/21))) ^ (1/data$k20[i]) ) ) - data$k20[i]  
  data$ro25[i]<-(data$k25[i] * ( (1/(1-(sum(prev25[i])/21))) ^ (1/data$k25[i]) ) ) - data$k25[i]
  data$ro30[i]<-(data$k30[i] * ( (1/(1-(sum(prev30[i])/21))) ^ (1/data$k30[i]) ) ) - data$k30[i]
  data$ro35[i]<-(data$k35[i] * ( (1/(1-(sum(prev35[i])/21))) ^ (1/data$k35[i]) ) ) - data$k35[i]
  data$ro40[i]<-(data$k40[i] * ( (1/(1-(sum(prev40[i])/21))) ^ (1/data$k40[i]) ) ) - data$k40[i]
  data$ro45[i]<-(data$k45[i] * ( (1/(1-(sum(prev45[i])/21))) ^ (1/data$k45[i]) ) ) - data$k45[i]
  data$ro50[i]<-(data$k50[i] * ( (1/(1-(sum(prev50[i])/21))) ^ (1/data$k50[i]) ) ) - data$k50[i]
  data$ro55[i]<-(data$k55[i] * ( (1/(1-(sum(prev55[i])/21))) ^ (1/data$k55[i]) ) ) - data$k55[i]
  data$ro60[i]<-(data$k60[i] * ( (1/(1-(sum(prev60[i])/21))) ^ (1/data$k60[i]) ) ) - data$k60[i]  
  data$ro65[i]<-(data$k65[i] * ( (1/(1-(sum(prev65[i])/21))) ^ (1/data$k65[i]) ) ) - data$k65[i]
  data$ro70[i]<-(data$k70[i] * ( (1/(1-(sum(prev70[i])/21))) ^ (1/data$k70[i]) ) ) - data$k70[i]
  data$ro75[i]<-(data$k75[i] * ( (1/(1-(sum(prev75[i])/21))) ^ (1/data$k75[i]) ) ) - data$k75[i]
  data$ro80[i]<-(data$k80[i] * ( (1/(1-(sum(prev80[i])/21))) ^ (1/data$k80[i]) ) ) - data$k80[i]
}
redr0<-c(mean(data$ro-data$ro01,na.rm=TRUE),mean(data$ro-data$ro05,na.rm=TRUE),mean(data$ro-data$ro10,na.rm=TRUE),
         mean(data$ro-data$ro15,na.rm=TRUE),mean(data$ro-data$ro20,na.rm=TRUE),mean(data$ro-data$ro25,na.rm=TRUE),
         mean(data$ro-data$ro30,na.rm=TRUE),mean(data$ro-data$ro35,na.rm=TRUE),mean(data$ro-data$ro40,na.rm=TRUE),
         mean(data$ro-data$ro45,na.rm=TRUE),mean(data$ro-data$ro50,na.rm=TRUE),mean(data$ro-data$ro55,na.rm=TRUE),
         mean(data$ro-data$ro60,na.rm=TRUE),mean(data$ro-data$ro65,na.rm=TRUE),mean(data$ro-data$ro70,na.rm=TRUE),
         mean(data$ro-data$ro75,na.rm=TRUE),rep(1,4))

redr0sd<-c(sd(data$ro-data$ro01,na.rm=TRUE),sd(data$ro-data$ro05,na.rm=TRUE),sd(data$ro-data$ro10,na.rm=TRUE),
           sd(data$ro-data$ro15,na.rm=TRUE),sd(data$ro-data$ro20,na.rm=TRUE),sd(data$ro-data$ro25,na.rm=TRUE),
           sd(data$ro-data$ro30,na.rm=TRUE),sd(data$ro-data$ro35,na.rm=TRUE),sd(data$ro-data$ro40,na.rm=TRUE),
           sd(data$ro-data$ro45,na.rm=TRUE),sd(data$ro-data$ro50,na.rm=TRUE),sd(data$ro-data$ro55,na.rm=TRUE),
           sd(data$ro-data$ro60,na.rm=TRUE),sd(data$ro-data$ro65,na.rm=TRUE),sd(data$ro-data$ro70,na.rm=TRUE),
           sd(data$ro-data$ro75,na.rm=TRUE),sd(data$ro-data$ro80,na.rm=TRUE),0.14,0.12,0.1)

par(mfrow=c(1,1))
plot(proportion,redr0,ylim=c(0,1),xlim=c(0,1),pch="",
     ylab=expression(paste("Reduction in  ",  R [0])),xaxt="n",
     xlab="Proportion of most infected hosts treated")
axis(1,at=seq(0,1,0.2),par(las=1),labels=seq(0,1,0.2))
uppred<-redr0+redr0sd;uppred<-ifelse(uppred>1,1,uppred)
#lines(uppred~proportion)
lowred<-redr0-redr0sd;lowred<-ifelse(lowred<0,0,lowred)
#lines(lowred~proportion)
polygon(c(proportion, rev(proportion)),c(lowred,rev(uppred)),border=NA, col=transp("darkseagreen1",alpha=0.3))
lines(redr0~proportion,lty=2)


## Least infected...
for (i in 1:ncol(distleast01)){
  sampmean<-mean(distleast01[,i])
  var<-var(distleast01[,i])
  data$k01least[i]<-(sampmean^2-(var/nrow(distleast01)))/(var-sampmean)
}
for (i in 1:ncol(distleast05)){
  sampmean<-mean(distleast05[,i])
  var<-var(distleast05[,i])
  data$k05least[i]<-(sampmean^2-(var/nrow(distleast05)))/(var-sampmean)
}
for (i in 1:ncol(distleast10)){
  sampmean<-mean(distleast10[,i])
  var<-var(distleast10[,i])
  data$k10least[i]<-(sampmean^2-(var/nrow(distleast10)))/(var-sampmean)
}
for (i in 1:ncol(distleast15)){
  sampmean<-mean(distleast15[,i])
  var<-var(distleast15[,i])
  data$k15least[i]<-(sampmean^2-(var/nrow(distleast15)))/(var-sampmean)
}
for (i in 1:ncol(distleast20)){
  sampmean<-mean(distleast20[,i])
  var<-var(distleast20[,i])
  data$k20least[i]<-(sampmean^2-(var/nrow(distleast20)))/(var-sampmean)
}
for (i in 1:ncol(distleast25)){
  sampmean<-mean(distleast25[,i])
  var<-var(distleast25[,i])
  data$k25least[i]<-(sampmean^2-(var/nrow(distleast25)))/(var-sampmean)
}
for (i in 1:ncol(distleast30)){
  sampmean<-mean(distleast30[,i])
  var<-var(distleast30[,i])
  data$k30least[i]<-(sampmean^2-(var/nrow(distleast30)))/(var-sampmean)
}
for (i in 1:ncol(distleast35)){
  sampmean<-mean(distleast35[,i])
  var<-var(distleast35[,i])
  data$k35least[i]<-(sampmean^2-(var/nrow(distleast35)))/(var-sampmean)
}
for (i in 1:ncol(distleast40)){
  sampmean<-mean(distleast40[,i])
  var<-var(distleast40[,i])
  data$k40least[i]<-(sampmean^2-(var/nrow(distleast40)))/(var-sampmean)
}
for (i in 1:ncol(distleast50)){
  sampmean<-mean(distleast50[,i])
  var<-var(distleast50[,i])
  data$k50least[i]<-(sampmean^2-(var/nrow(distleast50)))/(var-sampmean)
}
for (i in 1:ncol(distleast55)){
  sampmean<-mean(distleast55[,i])
  var<-var(distleast55[,i])
  data$k55least[i]<-(sampmean^2-(var/nrow(distleast55)))/(var-sampmean)
}
for (i in 1:ncol(distleast45)){
  sampmean<-mean(distleast45[,i])
  var<-var(distleast45[,i])
  data$k45least[i]<-(sampmean^2-(var/nrow(distleast45)))/(var-sampmean)
}
for (i in 1:ncol(distleast60)){
  sampmean<-mean(distleast60[,i])
  var<-var(distleast60[,i])
  data$k60least[i]<-(sampmean^2-(var/nrow(distleast60)))/(var-sampmean)
}
for (i in 1:ncol(distleast65)){
  sampmean<-mean(distleast65[,i])
  var<-var(distleast65[,i])
  data$k65least[i]<-(sampmean^2-(var/nrow(distleast65)))/(var-sampmean)
}
for (i in 1:ncol(distleast70)){
  sampmean<-mean(distleast70[,i])
  var<-var(distleast70[,i])
  data$k70least[i]<-(sampmean^2-(var/nrow(distleast70)))/(var-sampmean)
}
for (i in 1:ncol(distleast75)){
  sampmean<-mean(distleast75[,i])
  var<-var(distleast75[,i])
  data$k75least[i]<-(sampmean^2-(var/nrow(distleast75)))/(var-sampmean)
}
for (i in 1:ncol(distleast80)){
  sampmean<-mean(distleast80[,i])
  var<-var(distleast80[,i])
  data$k80least[i]<-(sampmean^2-(var/nrow(distleast80)))/(var-sampmean)
}
for (i in 1:ncol(distleast85)){
  sampmean<-mean(distleast85[,i])
  var<-var(distleast85[,i])
  data$k85least[i]<-(sampmean^2-(var/nrow(distleast85)))/(var-sampmean)
}
for (i in 1:ncol(distleast90)){
  sampmean<-mean(distleast90[,i])
  var<-var(distleast90[,i])
  data$k90least[i]<-(sampmean^2-(var/nrow(distleast90)))/(var-sampmean)
}
for (i in 1:ncol(distleast95)){
  sampmean<-mean(distleast95[,i])
  var<-var(distleast95[,i])
  data$k95least[i]<-(sampmean^2-(var/nrow(distleast95)))/(var-sampmean)
}
for (i in 1:ncol(distleast99)){
  sampmean<-mean(distleast99[,i])
  var<-var(distleast99[,i])
  data$k99least[i]<-(sampmean^2-(var/nrow(distleast99)))/(var-sampmean)
}
data$k01least[is.infinite(data$k01least)] <- NA
data$k05least[is.infinite(data$k05least)] <- NA
data$k10least[is.infinite(data$k10least)] <- NA
data$k15least[is.infinite(data$k15least)] <- NA
data$k20least[is.infinite(data$k20least)] <- NA
data$k25least[is.infinite(data$k25least)] <- NA
data$k30least[is.infinite(data$k30least)] <- NA
data$k35least[is.infinite(data$k35least)] <- NA
data$k40least[is.infinite(data$k40least)] <- NA
data$k45least[is.infinite(data$k45least)] <- NA
data$k50least[is.infinite(data$k50least)] <- NA
data$k55least[is.infinite(data$k55least)] <- NA
data$k60least[is.infinite(data$k60least)] <- NA
data$k65least[is.infinite(data$k65least)] <- NA
data$k70least[is.infinite(data$k70least)] <- NA
data$k75least[is.infinite(data$k75least)] <- NA
data$k80least[is.infinite(data$k80least)] <- NA
data$k85least[is.infinite(data$k85least)] <- NA
data$k90least[is.infinite(data$k75least)] <- NA
data$k95least[is.infinite(data$k80least)] <- NA
data$k99least[is.infinite(data$k85least)] <- NA


for (i in 1:206){
  data$lro[i]<-(data$k[i] * ( (1/(1-(data$prevalence[i]/data$N[i]))) ^ (1/data$k[i]) ) ) - data$k[i]
  
  data$lro01[i]<-(data$k01least[i] * ( (1/(1-(sum(leastgood01[i])/21))) ^ (1/data$k01least[i]) ) )  - data$k01least[i]
  data$lro05[i]<-(data$k05least[i] * ( (1/(1-(sum(leastgood05[i])/21))) ^ (1/data$k05least[i]) ) ) - data$k05least[i]
  data$lro10[i]<-(data$k10least[i] * ( (1/(1-(sum(leastgood10[i])/21))) ^ (1/data$k10least[i]) ) ) - data$k10least[i]
  data$lro15[i]<-(data$k15least[i] * ( (1/(1-(sum(leastgood15[i])/21))) ^ (1/data$k15least[i]) ) ) - data$k15least[i]
  data$lro20[i]<-(data$k20least[i] * ( (1/(1-(sum(leastgood20[i])/21))) ^ (1/data$k20least[i]) ) ) - data$k20least[i]  
  data$lro25[i]<-(data$k25least[i] * ( (1/(1-(sum(leastgood25[i])/21))) ^ (1/data$k25least[i]) ) ) - data$k25least[i]
  data$lro30[i]<-(data$k30least[i] * ( (1/(1-(sum(leastgood30[i])/21))) ^ (1/data$k30least[i]) ) ) - data$k30least[i]
  data$lro35[i]<-(data$k35least[i] * ( (1/(1-(sum(leastgood35[i])/21))) ^ (1/data$k35least[i]) ) ) - data$k35least[i]
  data$lro40[i]<-(data$k40least[i] * ( (1/(1-(sum(leastgood40[i])/21))) ^ (1/data$k40least[i]) ) ) - data$k40least[i]
  data$lro45[i]<-(data$k45least[i] * ( (1/(1-(sum(leastgood45[i])/21))) ^ (1/data$k45least[i]) ) ) - data$k45least[i]
  data$lro50[i]<-(data$k50least[i] * ( (1/(1-(sum(leastgood50[i])/21))) ^ (1/data$k50least[i]) ) ) - data$k50least[i]
  data$lro55[i]<-(data$k55least[i] * ( (1/(1-(sum(leastgood55[i])/21))) ^ (1/data$k55least[i]) ) ) - data$k55least[i]
  data$lro60[i]<-(data$k60least[i] * ( (1/(1-(sum(leastgood60[i])/21))) ^ (1/data$k60least[i]) ) ) - data$k60least[i]  
  data$lro65[i]<-(data$k65least[i] * ( (1/(1-(sum(leastgood65[i])/21))) ^ (1/data$k65least[i]) ) ) - data$k65least[i]
  data$lro70[i]<-(data$k70least[i] * ( (1/(1-(sum(leastgood70[i])/21))) ^ (1/data$k70least[i]) ) ) - data$k70least[i]
  data$lro75[i]<-(data$k75least[i] * ( (1/(1-(sum(leastgood75[i])/21))) ^ (1/data$k75least[i]) ) ) - data$k75least[i]
  data$lro80[i]<-(data$k80least[i] * ( (1/(1-(sum(leastgood80[i])/21))) ^ (1/data$k80least[i]) ) ) - data$k80least[i]
}
redr02<-c(mean(data$lro-data$lro01,na.rm=TRUE),mean(data$lro-data$lro05,na.rm=TRUE),mean(data$lro-data$lro10,na.rm=TRUE),
         mean(data$lro-data$lro15,na.rm=TRUE),mean(data$lro-data$lro20,na.rm=TRUE),mean(data$lro-data$lro25,na.rm=TRUE),
         mean(data$lro-data$lro30,na.rm=TRUE),mean(data$lro-data$lro35,na.rm=TRUE),mean(data$lro-data$lro40,na.rm=TRUE),
         mean(data$lro-data$lro45,na.rm=TRUE),mean(data$lro-data$lro50,na.rm=TRUE),mean(data$lro-data$lro55,na.rm=TRUE),
         mean(data$lro-data$lro60,na.rm=TRUE),mean(data$lro-data$lro65,na.rm=TRUE),mean(data$lro-data$lro70,na.rm=TRUE),
         mean(data$lro-data$lro75,na.rm=TRUE),mean(data$lro-data$lro80,na.rm=TRUE),mean(data$lro-data$lro85,na.rm=TRUE),
         mean(data$lro-data$lro90,na.rm=TRUE),mean(data$lro-data$lro95,na.rm=TRUE),mean(data$lro-data$lro99,na.rm=TRUE))

redr0sd2<-c(sd(data$lro-data$lro01,na.rm=TRUE),sd(data$lro-data$lro05,na.rm=TRUE),sd(data$lro-data$lro10,na.rm=TRUE),
           sd(data$lro-data$lro15,na.rm=TRUE),sd(data$lro-data$lro20,na.rm=TRUE),sd(data$lro-data$lro25,na.rm=TRUE),
           sd(data$lro-data$lro30,na.rm=TRUE),sd(data$lro-data$lro35,na.rm=TRUE),sd(data$lro-data$lro40,na.rm=TRUE),
           sd(data$lro-data$lro45,na.rm=TRUE),sd(data$lro-data$lro50,na.rm=TRUE),sd(data$lro-data$lro55,na.rm=TRUE),
           sd(data$lro-data$lro60,na.rm=TRUE),sd(data$lro-data$lro65,na.rm=TRUE),sd(data$lro-data$lro70,na.rm=TRUE),
           sd(data$lro-data$lro75,na.rm=TRUE),sd(data$lro-data$lro80,na.rm=TRUE),0.14,0.12,0.1)

uppred2<-redr02+redr0sd2;uppred2<-ifelse(uppred2>1,1,uppred2)
lowred2<-redr02-redr0sd2;lowred2<-ifelse(lowred2<0,0,lowred2)
polygon(c(proportion, rev(proportion)),c(lowred2,rev(uppred2)),border=NA, col=transp("darkseagreen1",alpha=0.1))
lines(redr0~proportion,lty=2)