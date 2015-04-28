
library(rstan)
library(MASS)
library(boot)
library(coda)
library(R2OpenBUGS)
library(ggplot2)
library(reshape2)

data<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\Exploring R0.csv",header=TRUE)
head(data);summary(data)

############################################################
## What impact on prevalence occurs if a given proportion of
## the hosts are treated (ie the most infection is removed)
##
##

data$prev01T<-(data$prevalence-(0.01*data$N))/(data$N-data$N*0.01);data$prev01T<-ifelse(data$prev01T<0,0,data$prev01T)
data$prev05T<-(data$prevalence-(0.05*data$N))/(data$N-data$N*0.05);data$prev05T<-ifelse(data$prev05T<0,0,data$prev05T)

data$prev10T<-(data$prevalence-(0.1*data$N))/(data$N-data$N*0.1);data$prev10T<-ifelse(data$prev10T<0,0,data$prev10T)
data$prev15T<-(data$prevalence-(0.15*data$N))/(data$N-data$N*0.15);data$prev15T<-ifelse(data$prev15T<0,0,data$prev15T)

data$prev20T<-(data$prevalence-(0.2*data$N))/(data$N-data$N*0.2);data$prev20T<-ifelse(data$prev20T<0,0,data$prev20T)
data$prev25T<-(data$prevalence-(0.25*data$N))/(data$N-data$N*0.25);data$prev25T<-ifelse(data$prev25T<0,0,data$prev25T)

data$prev30T<-(data$prevalence-(0.3*data$N))/(data$N-data$N*0.3);data$prev30T<-ifelse(data$prev30T<0,0,data$prev30T)
data$prev35T<-(data$prevalence-(0.35*data$N))/(data$N-data$N*0.35);data$prev35T<-ifelse(data$prev35T<0,0,data$prev35T)

data$prev40T<-(data$prevalence-(0.4*data$N))/(data$N-data$N*0.4);data$prev40T<-ifelse(data$prev40T<0,0,data$prev40T)
data$prev45T<-(data$prevalence-(0.45*data$N))/(data$N-data$N*0.45);data$prev45T<-ifelse(data$prev45T<0,0,data$prev45T)

data$prev50T<-(data$prevalence-(0.5*data$N))/(data$N-data$N*0.5);data$prev50T<-ifelse(data$prev50T<0,0,data$prev50T)
data$prev55T<-(data$prevalence-(0.55*data$N))/(data$N-data$N*0.55);data$prev55T<-ifelse(data$prev55T<0,0,data$prev55T)

data$prev60T<-(data$prevalence-(0.6*data$N))/(data$N-data$N*0.6);data$prev60T<-ifelse(data$prev60T<0,0,data$prev60T)
data$prev65T<-(data$prevalence-(0.65*data$N))/(data$N-data$N*0.65);data$prev65T<-ifelse(data$prev65T<0,0,data$prev65T)

data$prev70T<-(data$prevalence-(0.70*data$N))/(data$N-data$N*0.70);data$prev70T<-ifelse(data$prev70T<0,0,data$prev70T)
data$prev75T<-(data$prevalence-(0.75*data$N))/(data$N-data$N*0.75);data$prev75T<-ifelse(data$prev75T<0,0,data$prev75T)

data$prev80T<-(data$prevalence-(0.80*data$N))/(data$N-data$N*0.80);data$prev80T<-ifelse(data$prev80T<0,0,data$prev80T)
data$prev85T<-(data$prevalence-(0.85*data$N))/(data$N-data$N*0.85);data$prev85T<-ifelse(data$prev85T<0,0,data$prev85T)

data$prev90T<-(data$prevalence-(0.90*data$N))/(data$N-data$N*0.90);data$prev90T<-ifelse(data$prev90T<0,0,data$prev90T)
data$prev95T<-(data$prevalence-(0.95*data$N))/(data$N-data$N*0.95);data$prev95T<-ifelse(data$prev95T<0,0,data$prev95T)

data$prev99T<-(data$prevalence-(0.99*data$N))/(data$N-data$N*0.99);data$prev99T<-ifelse(data$prev99T<0,0,data$prev99T)

dim(data)
prevalence<-c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99)
reduction<-as.numeric(c(data[1,35:55]))
plot(reduction~prevalence,ylim=c(0,1))
redu<-matrix(nrow=length(prevalence),ncol=length(data$N),data=NA)
for (i in 1:length(data$N)){
redu[,i]<-as.numeric(c(data[i,35:55]))
lines(redu[,i]~prevalence)
}


###############################################
## Figure 1

par(mfrow=c(1,2))
par(mar=c(5,5,2,2))
parasitesTreated<-as.numeric(data[1,14:34])
plot(parasitesTreated~prevalence,pch="",
     ylim=c(0,1),xlim=c(0,1),cex.lab=1.1,
     ylab="Proportion of parasites",xlab="Proportion of hosts treated")
treated<-matrix(nrow=length(prevalence),ncol=length(data$N),data=NA)
for (i in 1:length(data$N)){
  treated[,i]<-as.numeric(data[i,14:34])
  lines(c(treated[,i])~prevalence,col="gray",)
}
meantreated<-as.numeric(21)
for (i in 1:21){ 
  meantreated[i]<-mean(treated[i,])
}
lines(c(0,meantreated)~c(0,prevalence),lty=2,lwd=2,col="red")
    
##Create a distribution of the count data estimated for each population
dim(data);head(data)
distribs<-matrix(nrow=21,ncol=206)
colnames(distribs)<-data$Label[1:206]
for (j in 1:206){
for (i in 14:34){
  distribs[i-13,j]<-data$N[j]-round(data$N[j]*data[j,i])
}
}
distribsALL<-distribs

##and the distribution were x % treated
disttemp<-distribs;disttemp[1,]<-0;distrib01<-disttemp
disttemp[1:2,]<-0;distrib05<-disttemp
disttemp[1:3,]<-0;distrib10<-disttemp
disttemp[1:4,]<-0;distrib15<-disttemp
disttemp[1:5,]<-0;distrib20<-disttemp
disttemp[1:6,]<-0;distrib25<-disttemp
disttemp[1:7,]<-0;distrib30<-disttemp
disttemp[1:8,]<-0;distrib35<-disttemp
disttemp[1:9,]<-0;distrib40<-disttemp
disttemp[1:10,]<-0;distrib45<-disttemp
disttemp[1:11,]<-0;distrib50<-disttemp
disttemp[1:12,]<-0;distrib55<-disttemp
disttemp[1:13,]<-0;distrib60<-disttemp
disttemp[1:14,]<-0;distrib65<-disttemp
disttemp[1:15,]<-0;distrib70<-disttemp
disttemp[1:16,]<-0;distrib75<-disttemp
disttemp[1:17,]<-0;distrib80<-disttemp
disttemp[1:18,]<-0;distrib85<-disttemp
disttemp[1:19,]<-0;distrib90<-disttemp
disttemp[1:20,]<-0;distrib95<-disttemp
disttemp[1:21,]<-0;distrib99<-disttemp

prevdist<-ifelse(distribs==0,0,1)
prevs<-numeric(ncol(distribs))
for (i in 1:ncol(distribs)) {
  prevs[i]<-sum(prevdist[,i])
}

##and the distribution were x % treated
disttemp<-prevdist;disttemp[1,]<-0;prev01<-disttemp ##1%
disttemp[1:2,]<-0;prev05<-disttemp
disttemp[1:3,]<-0;prev10<-disttemp
disttemp[1:4,]<-0;prev15<-disttemp
disttemp[1:5,]<-0;prev20<-disttemp ##20% target treated
disttemp[1:6,]<-0;prev25<-disttemp
disttemp[1:7,]<-0;prev30<-disttemp
disttemp[1:8,]<-0;prev35<-disttemp
disttemp[1:9,]<-0;prev40<-disttemp
disttemp[1:10,]<-0;prev45<-disttemp
disttemp[1:11,]<-0;prev50<-disttemp
disttemp[1:12,]<-0;prev55<-disttemp
disttemp[1:13,]<-0;prev60<-disttemp
disttemp[1:14,]<-0;prev65<-disttemp
disttemp[1:15,]<-0;prev70<-disttemp
disttemp[1:16,]<-0;prev75<-disttemp
disttemp[1:17,]<-0;prev80<-disttemp
disttemp[1:18,]<-0;prev85<-disttemp
disttemp[1:19,]<-0;prev90<-disttemp
disttemp[1:20,]<-0;prev95<-disttemp
disttemp[1:21,]<-0;prev99<-disttemp

###############################################################
## Target treating the most infected xx% as above
disttemp<-distribs;disttemp[1:5,]<-0;distrib20<-disttemp
disttemp<-prevdist;disttemp[1:5,]<-0;prev20<-disttemp ##20% target treated
###############################################################
## Target treating the least infected xx% 
disttemp<-distribs;disttemp[18:21,]<-0;distleast80<-disttemp
disttemp<-prevdist;disttemp[18:21,]<-0;leastgood20<-disttemp
###############################################################
## Randomly treating xx% infected
datatest<-randtest<-expand.grid(c(1:21))
disttemp<-distribs;for (i in 1:206){
  datatest[,i]<-sample(disttemp[,i])
  datatest[1:5,i]<-0;randtest[,i]<-datatest[,i]
  }
randprev<-expand.grid(c(1:21))
disttemp<-prevdist;for (i in 1:206){
  datatest[,i]<-sample(disttemp[,i])
  datatest[1:5,i]<-0;randprev[,i]<-datatest[,i]
}
colnames(randtest)<-colnames(randprev)<-data$Label
###Take for example the first case


###################################################################
##################################################################### AUTOMATED
#######################################################################
dataoutmean_i<-expand.grid(c(1:4))
dataout95upper_i<-expand.grid(c(1:4))
dataout95lower_i<-expand.grid(c(1:4))
for (i in 1:206){
  data2<-list(N_st=4,
              N_counts=21,
              para_count = structure(.Data = c(distribsALL[,i],distrib20[,i],distleast80[,i],randtest[,i]),
                                     .Dim=c(21,4)),###[N_counts,N_st]
              prev = structure(.Data =c(prevdist[,i],prev20[,i],leastgood20[,i],randtest[,i]),
                               .Dim=c(21,4)))

  fit1 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080\\modelA2.stan", data=data2,
               iter=1000, chains=2)
  
  #print(fit1)
  #data$Label
  #proportions<-c(0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
params = extract(fit1);names(params)
  dataoutmean_i[,i]<-c(mean(params$theta[,1]),mean(params$theta[,2]),mean(params$theta[,3]),mean(params$theta[,4]))
  dataout95upper_i[,i]<-as.numeric(c(quantile(params$theta[,1],0.975),quantile(params$theta[,2],0.975),quantile(params$theta[,3],0.975),quantile(params$theta[,4],0.975)))
  dataout95lower_i[,i]<-as.numeric(c(quantile(params$theta[,1],0.025),quantile(params$theta[,2],0.025),quantile(params$theta[,3],0.025),quantile(params$theta[,4],0.025)))
  
}

#write.csv(dataoutmean_i,"C:\\Users\\Ellie\\Documents\\2080\\model_outputMEAN_i2.csv")
#write.csv(dataout95upper_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95upper_i2.csv")
#write.csv(dataout95lower_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95lower_i2.csv")
treatments<-rep(c("0","1",seq(from=5,to=95,by=5)),each=21)


dataoutmean<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\model_outputMEAN_i22.csv")
dataout95upper<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\model_output95upper_i22.csv")
dataout95lower<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\model_output95lower_i22.csv")
colnames(dataoutmean)<-colnames(dataout95lower)<-colnames(dataout95upper)<-c("variable",data$Label[1:206])

datacheck2<-melt(dataoutmean,id=1)##
colnames(datacheck2)<-c("prop","study","theta")
datacheck2$proportion20<-rep(data$T20,each="4")
ggplot(datacheck2) + geom_point(aes(x=proportion20,y=theta, col=prop))


alldata<-as.numeric(dataoutmean[1,2:207])
rem20<-as.numeric(dataoutmean[2,2:207])
least20<-as.numeric(dataoutmean[3,2:207])
rand20<-as.numeric(dataoutmean[4,2:207])

alldataU<-as.numeric(dataout95upper[1,2:207])
rem20U<-as.numeric(dataout95upper[2,2:207])
least20U<-as.numeric(dataout95upper[3,2:207])
rand20U<-as.numeric(dataout95upper[4,2:207])

alldataL<-as.numeric(dataout95lower[1,2:207])
rem20L<-as.numeric(dataout95lower[2,2:207])
least20L<-as.numeric(dataout95lower[3,2:207])
rand20L<-as.numeric(dataout95lower[4,2:207])

par(mfrow=c(1,1))
par(mar=c(5,5,2,2))
xv<-data$T20[1:206]
#xv<-data$prevalence[1:206]/data$N[1:206]
plot(alldata~xv,ylab="Probability of infection",ylim=c(0,1),
     xlab="Proportion of parasites in 20% of hosts",xlim=c(0,1),pch="",cex.lab=1.2)

for (i in 1:206){
  segments(xv[i], alldataU[i], x1 = xv[i], y1 = alldataL[i],
           col  = terrain.colors(2,alpha = 0.1), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(xv[i], rem20U[i], x1 = xv[i], y1 = rem20L[i],
           col  = terrain.colors(10,alpha = 0.5), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(xv[i], least20U[i], x1 = xv[i], y1 = least20L[i],
           col  = terrain.colors(25,alpha = 0.8), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(xv[i], rand20U[i], x1 = xv[i], y1 = rand20L[i],
           col  = terrain.colors(25,alpha = 0.8), lty = 1, lwd = 5)
}

points(alldata~xv,col="lightblue",pch=20);points(least20~xv,col="black",pch=20)
points(rem20~xv,col="blue",pch=20);points(rand20~xv,col="blue",pch=20)

legend(0,1,legend=c("No treatment","Treat 20%","Treat 50%"),
       col=c("lightblue","blue","black"),
       pch=20,bty="n")


###################################################################
##
###
####  Now think about the reduction in transmission potential that is acheived with x% treatment
###
##

reductionacheived<-expand.grid(seq(1,3))
for(i in 2:207){
  for (j in 2:4){
    reductionacheived[j-1,i-1]<-dataoutmean[1,i]-dataoutmean[j,i]
  }
}
colnames(reductionacheived)<-data$Label[1:206]

reductionach95U<-expand.grid(seq(1,3))
for(i in 2:207){
  for (j in 2:4){
    reductionach95U[j-1,i-1]<-dataout95upper[1,i]-dataout95upper[j,i]
  }
}
colnames(reductionach95U)<-data$Label[1:206]

reductionach95L<-expand.grid(seq(1,3))
for(i in 2:207){
  for (j in 2:4){
    reductionach95L[j-1,i-1]<-dataout95lower[1,i]-dataout95lower[j,i]
  }
}
colnames(reductionach95L)<-data$Label[1:206]

##percentage reduction acheived
perredacheived<-expand.grid(seq(1,3))
for(i in 2:207){
  for (j in 2:4){
    perredacheived[j-1,i-1]<-(dataoutmean[1,i]-dataoutmean[j,i])/dataoutmean[1,i]
  }
}

perredach95U<-expand.grid(seq(1,3))
for(i in 2:207){
  for (j in 2:4){
    perredach95U[j-1,i-1]<-(dataout95upper[1,i]-dataout95upper[j,i])/dataout95upper[1,i]
  }
}

perredach95L<-expand.grid(seq(1,3))
for(i in 2:207){
  for (j in 2:4){
    perredach95L[j-1,i-1]<-(dataout95lower[1,i]-dataout95lower[j,i])/dataout95lower[1,i]
  }
}
colnames(perredacheived)<-colnames(perredach95U)<-colnames(perredach95L)<-data$Label[1:206]

rem20red<-as.numeric(perredacheived[1,1:206])
remleast20red<-as.numeric(perredacheived[2,1:206])
remrand20red<-as.numeric(perredacheived[3,1:206])

rem20redU<-as.numeric(perredach95U[1,1:206])
remleast20redU<-as.numeric(perredach95U[2,1:206])
remrand20redU<-as.numeric(perredach95U[3,1:206])

rem20redL<-as.numeric(perredach95L[1,1:206])
remleast20redL<-as.numeric(perredach95L[2,1:206])
remrand20redL<-as.numeric(perredach95L[3,1:206])

D1<-data.frame(rem20red,rem20redU,rem20redL,remleast20red,remleast20redU,remleast20redL,
               remrand20red,remrand20redU,remrand20redL)
##D2<-ordered(D1,D1$rem20red) ## order the data by the ones in top 20

par(mfrow=c(1,1))
par(mar=c(5,5,2,2))
xv<-data$Label[1:206]
#xv<-data$prevalence[1:206]/data$N[1:206]
plot(rem20red~xv,ylab=expression(paste("Proportionate reduction in  ",theta)),
     ylim=c(0,1),
     xlab="Studies",
     pch="",cex.lab=1.2)

for (i in 1:206){
  segments(xv[i], rem20redU[i], x1 = xv[i], y1 = rem20redL[i],
           col  = terrain.colors(10,alpha = 0.7), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(xv[i], remleast20redU[i], x1 = xv[i], y1 = remleast20redL[i],
           col  = terrain.colors(10,alpha = 0.3), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(xv[i], remrand20redU[i], x1 = xv[i], y1 = remrand20redL[i],
           col  = terrain.colors(10,alpha = 0.1), lty = 1, lwd = 5)
}

points(rem20red~xv,col="lightblue",pch=20);points(remleast20red~xv,col="black",pch=20);points(remrand20red~xv,col="blue",pch=20)

legend(250,0.4,legend=c("Targeted most infected","Least infected","Randomly infected"),
       col=c("lightblue","blue","black"),
       pch=20,bty="n")

sum(ifelse(rem20red>0.95,1,0))/206 ##A 95% or better reduction in transmission probability can be acheived in 32% of the studies
sum(ifelse(rem50red>0.95,1,0))/206 

###################################################################################
## Figure 2a
par(mfrow=c(1,1))
plot(log(data$k+1),data$T10,pch=20,col="grey15",
     ylab="Proportion of parasites in x% of hosts",
     xlab="Log k",cex.lab=1.1,log="x",xlim=c(0.01,5))
points(log(data$k+1),data$T20,col="grey30",pch=20,cex=2)
points(log(data$k+1),data$T30,col="grey50",pch=17)
points(log(data$k+1),data$T50,col="grey85",pch=15)

legend(0.01,0.6,
       legend=c(expression("t"[10]),expression("t"[20]),
       expression("t"[30]),expression("t"[50])),
       col=c("grey15","grey30","grey50","grey85"),
       pch=c(20,20,17,15),cex=1.4)

abline(v=0.055,lty=2);text(0.045,0.2,"k: 1.057")

abline(v=0.03,lty=2);text(0.024,0.2,"k: 1.030")

abline(v=0.1,lty=2);text(0.124,0.2,"k: 1.105")

abline(v=0.2,lty=2);text(0.245,0.2,"k: 1.221")


names(data)
dat2<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\DataFinalApril2015.csv",header=TRUE)
names(dat2)
data3<-merge(data,dat2,by.x="Label",by.y="Label")
dim(data3)
names(data3)
## Figure 2b
par(mfrow=c(1,1))
plot(1-data3$Gini.Co.efficient,data$T10,pch=20,col="grey15",
     ylab="Proportion of parasites in x% of hosts",
     xlab="Gini Co-efficient",cex.lab=1.1,xlim=c(0,1))
points(1-data3$Gini.Co.efficient,data$T20,col="grey30",pch=20,cex=2)
points(1-data3$Gini.Co.efficient,data$T30,pch=20,col="grey50")
points(1-data3$Gini.Co.efficient,data$T50,pch=20,col="grey80")

abline(v=0.06,lty=2);text(0.04,0.2,"0.06")

abline(v=0.1,lty=2);text(0.08,0.25,"0.10")

abline(v=0.14,lty=2);text(0.12,0.3,"0.14")

abline(v=0.22,lty=2);text(0.2,0.35,"0.22")

legend(0.84,1,
       legend=c(expression("t"[10]),expression("t"[20]),
                expression("t"[30]),expression("t"[50])),
       col=c("grey15","grey30","grey50","grey85"),
       pch=c(20,20,17,15),cex=1.4)



ggplot(plotdat) + geom_point(aes(x=prop,y=Theta, col=treatments), alpha = 0.5)
ggplot(plotdat) + geom_jitter(aes(x=prop,y=Theta, col=treatments))
ggplot(plotdat) + geom_violin(aes(x=prop,y=Theta, col=treatments))