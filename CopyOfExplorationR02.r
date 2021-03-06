
library(rstan)
library(MASS)
library(boot)
library(coda)
library(R2OpenBUGS)
library(ggplot2)
library(reshape2)
library(adegenet)

data<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\Exploring R0.csv",header=TRUE)
head(data);summary(data);dim(data)

############################################################
## What impact on prevalence occurs if a given proportion of
## the hosts are treated (ie the most infection is removed)
##
##
                    
###############################################
## Figure 1
prevalence<-c(0.01,seq(0.05,0.95,0.05),0.99)

par(mfrow=c(1,2))
par(mar=c(5,5,2,2))
parasitesTreated<-as.numeric(data[1,14:34])
plot(parasitesTreated~prevalence,pch="",frame.plot=F,
     ylim=c(0,1),xlim=c(0,1),cex.lab=1.1,xaxt="n",
     ylab="Proportion of parasites",xlab="Proportion of hosts treated")
par(las=1)
axis(1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,0.2,0.4,0.6,0.8,1.0))
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
disttemp<-distribs;disttemp[21,]<-0;distleast01<-disttemp
disttemp<-distribs;disttemp[20:21,]<-0;distleast05<-disttemp
disttemp<-distribs;disttemp[19:21,]<-0;distleast10<-disttemp
disttemp<-distribs;disttemp[18:21,]<-0;distleast15<-disttemp
disttemp<-distribs;disttemp[17:21,]<-0;distleast20<-disttemp
disttemp<-distribs;disttemp[16:21,]<-0;distleast25<-disttemp
disttemp<-distribs;disttemp[15:21,]<-0;distleast30<-disttemp
disttemp<-distribs;disttemp[14:21,]<-0;distleast35<-disttemp
disttemp<-distribs;disttemp[13:21,]<-0;distleast40<-disttemp
disttemp<-distribs;disttemp[12:21,]<-0;distleast45<-disttemp
disttemp<-distribs;disttemp[11:21,]<-0;distleast50<-disttemp
disttemp<-distribs;disttemp[10:21,]<-0;distleast55<-disttemp
disttemp<-distribs;disttemp[9:21,]<-0;distleast60<-disttemp
disttemp<-distribs;disttemp[8:21,]<-0;distleast65<-disttemp
disttemp<-distribs;disttemp[7:21,]<-0;distleast70<-disttemp
disttemp<-distribs;disttemp[6:21,]<-0;distleast75<-disttemp
disttemp<-distribs;disttemp[5:21,]<-0;distleast80<-disttemp
disttemp<-distribs;disttemp[4:21,]<-0;distleast85<-disttemp
disttemp<-distribs;disttemp[3:21,]<-0;distleast90<-disttemp
disttemp<-distribs;disttemp[2:21,]<-0;distleast95<-disttemp
disttemp<-distribs;disttemp[1:21,]<-0;distleast99<-disttemp

disttemp<-prevdist;disttemp[21,]<-0;leastgood01<-disttemp
disttemp<-prevdist;disttemp[20:21,]<-0;leastgood05<-disttemp
disttemp<-prevdist;disttemp[19:21,]<-0;leastgood10<-disttemp
disttemp<-prevdist;disttemp[18:21,]<-0;leastgood15<-disttemp
disttemp<-prevdist;disttemp[17:21,]<-0;leastgood20<-disttemp
disttemp<-prevdist;disttemp[16:21,]<-0;leastgood25<-disttemp
disttemp<-prevdist;disttemp[15:21,]<-0;leastgood30<-disttemp
disttemp<-prevdist;disttemp[14:21,]<-0;leastgood35<-disttemp
disttemp<-prevdist;disttemp[13:21,]<-0;leastgood40<-disttemp
disttemp<-prevdist;disttemp[12:21,]<-0;leastgood45<-disttemp
disttemp<-prevdist;disttemp[11:21,]<-0;leastgood50<-disttemp
disttemp<-prevdist;disttemp[10:21,]<-0;leastgood55<-disttemp
disttemp<-prevdist;disttemp[9:21,]<-0;leastgood60<-disttemp
disttemp<-prevdist;disttemp[8:21,]<-0;leastgood65<-disttemp
disttemp<-prevdist;disttemp[7:21,]<-0;leastgood70<-disttemp
disttemp<-prevdist;disttemp[6:21,]<-0;leastgood75<-disttemp
disttemp<-prevdist;disttemp[5:21,]<-0;leastgood80<-disttemp
disttemp<-prevdist;disttemp[4:21,]<-0;leastgood85<-disttemp
disttemp<-prevdist;disttemp[3:21,]<-0;leastgood90<-disttemp
disttemp<-prevdist;disttemp[2:21,]<-0;leastgood95<-disttemp
disttemp<-prevdist;disttemp[1:21,]<-0;leastgood99<-disttemp

###############################################################
## Randomly treating xx% infected
datatest<-dataprev<-expand.grid(c(1:21))
randtest01<-randtest05<-randtest10<-randtest15<-randtest20<-randtest25<-expand.grid(c(1:21))
randtest30<-randtest35<-randtest40<-randtest45<-randtest50<-expand.grid(c(1:21))
randtest55<-randtest60<-randtest65<-randtest70<-randtest75<-expand.grid(c(1:21))
randtest80<-randtest85<-randtest90<-randtest95<-randtest99<-expand.grid(c(1:21))
disttemp<-distribs;for (i in 1:206){
  datatest[,i]<-sample(disttemp[,i])
  datatest[1,i]<-0;randtest01[,i]<-datatest[,i]
  datatest[1:2,i]<-0;randtest05[,i]<-datatest[,i]
  datatest[1:3,i]<-0;randtest10[,i]<-datatest[,i]
  datatest[1:4,i]<-0;randtest15[,i]<-datatest[,i]
  datatest[1:5,i]<-0;randtest20[,i]<-datatest[,i]
  datatest[1:6,i]<-0;randtest25[,i]<-datatest[,i]
  datatest[1:7,i]<-0;randtest30[,i]<-datatest[,i]
  datatest[1:8,i]<-0;randtest35[,i]<-datatest[,i]
  datatest[1:9,i]<-0;randtest40[,i]<-datatest[,i]
  datatest[1:10,i]<-0;randtest45[,i]<-datatest[,i]
  datatest[1:11,i]<-0;randtest50[,i]<-datatest[,i]
  datatest[1:12,i]<-0;randtest55[,i]<-datatest[,i]
  datatest[1:13,i]<-0;randtest60[,i]<-datatest[,i]
  datatest[1:14,i]<-0;randtest65[,i]<-datatest[,i]
  datatest[1:15,i]<-0;randtest70[,i]<-datatest[,i]
  datatest[1:16,i]<-0;randtest75[,i]<-datatest[,i]
  datatest[1:17,i]<-0;randtest80[,i]<-datatest[,i]
  datatest[1:18,i]<-0;randtest85[,i]<-datatest[,i]
  datatest[1:19,i]<-0;randtest90[,i]<-datatest[,i]
  datatest[1:20,i]<-0;randtest95[,i]<-datatest[,i]
  datatest[1:21,i]<-0;randtest99[,i]<-datatest[,i]
  }
randprev01<-randprev05<-randprev10<-randprev15<-randprev20<-randprev25<-expand.grid(c(1:21))
randprev30<-randprev35<-randprev40<-randprev45<-randprev50<-expand.grid(c(1:21))
randprev55<-randprev60<-randprev65<-randprev70<-randprev75<-expand.grid(c(1:21))
randprev80<-randprev85<-randprev90<-randprev95<-randprev99<-expand.grid(c(1:21))
disttemp<-prevdist;for (i in 1:206){
  dataprev[,i]<-sample(disttemp[,i])
  dataprev[1,i]<-0;randprev01[,i]<-dataprev[,i]
  dataprev[1:2,i]<-0;randprev05[,i]<-dataprev[,i]
  dataprev[1:3,i]<-0;randprev10[,i]<-dataprev[,i]
  dataprev[1:4,i]<-0;randprev15[,i]<-dataprev[,i]
  dataprev[1:5,i]<-0;randprev20[,i]<-dataprev[,i]
  dataprev[1:6,i]<-0;randprev25[,i]<-dataprev[,i]
  dataprev[1:7,i]<-0;randprev30[,i]<-dataprev[,i]
  dataprev[1:8,i]<-0;randprev35[,i]<-dataprev[,i]
  dataprev[1:9,i]<-0;randprev40[,i]<-dataprev[,i]
  dataprev[1:10,i]<-0;randprev45[,i]<-dataprev[,i]
  dataprev[1:11,i]<-0;randprev50[,i]<-dataprev[,i]
  dataprev[1:12,i]<-0;randprev55[,i]<-dataprev[,i]
  dataprev[1:13,i]<-0;randprev60[,i]<-dataprev[,i]
  dataprev[1:14,i]<-0;randprev65[,i]<-dataprev[,i]
  dataprev[1:15,i]<-0;randprev70[,i]<-dataprev[,i]
  dataprev[1:16,i]<-0;randprev75[,i]<-dataprev[,i]
  dataprev[1:17,i]<-0;randprev80[,i]<-dataprev[,i]
  dataprev[1:18,i]<-0;randprev85[,i]<-dataprev[,i]
  dataprev[1:19,i]<-0;randprev90[,i]<-dataprev[,i]
  dataprev[1:20,i]<-0;randprev95[,i]<-dataprev[,i]
  dataprev[1:21,i]<-0;randprev99[,i]<-dataprev[,i]
}
#colnames(randtest)<-colnames(randprev)<-data$Label



###################################################################
##################################################################### AUTOMATED
#######################################################################
dataoutmean_i<-data3outmean_i<-dataoutvar_i<-data4outmean_i<-dataoutupper_i<-dataoutlower_i<-expand.grid(c(1:21))
dataout95upper_i<-data3out95upper_i<-data4out95upper_i<-expand.grid(c(1:21))
dataout95lower_i<-data3out95lower_i<-data4out95lower_i<-expand.grid(c(1:21))
for (i in 1:14){
  data2<-list(N_st=21,
              N_counts=21,
              para_count = structure(.Data = c(distribsALL[,i],distrib01[,i],distrib05[,i],distrib10[,i],distrib15[,i],
                                               distrib20[,i],distrib25[,i],distrib30[,i],distrib35[,i],
                                               distrib40[,i],distrib45[,i],distrib50[,i],distrib55[,i],
                                               distrib60[,i],distrib65[,i],distrib70[,i],distrib75[,i],
                                               distrib80[,i],distrib85[,i],distrib90[,i],distrib95[,i]),
                                     .Dim=c(21,21)),###[N_counts,N_st]
              prev = structure(.Data =c(prevdist[,i],prev01[,i],prev05[,i],prev10[,i],prev15[,i],
                                        prev20[,i],prev25[,i],prev30[,i],prev35[,i],
                                        prev40[,i],prev45[,i],prev50[,i],prev55[,i],
                                        prev60[,i],prev65[,i],prev70[,i],prev75[,i],
                                        prev80[,i],prev85[,i],prev90[,i],prev95[,i]),
                               .Dim=c(21,21)))
  
  data3<-list(N_st=21,
              N_counts=21,
              para_count = structure(.Data = c(distribsALL[,i],
    distleast01[,i],distleast05[,i],distleast10[,i],distleast15[,i],
  distleast20[,i],distleast25[,i],distleast30[,i],distleast35[,i],
  distleast40[,i],distleast45[,i],distleast50[,i],distleast55[,i],
  distleast60[,i],distleast65[,i],distleast70[,i],distleast75[,i],
  distleast80[,i],distleast85[,i],distleast90[,i],distleast95[,i]),
  .Dim=c(21,21)),###[N_counts,N_st]
  prev = structure(.Data =c(prevdist[,i],
    leastgood01[,i],leastgood05[,i],leastgood10[,i],leastgood15[,i],
    leastgood20[,i],leastgood25[,i],leastgood30[,i],leastgood35[,i],
    leastgood40[,i],leastgood45[,i],leastgood50[,i],leastgood55[,i],
    leastgood60[,i],leastgood65[,i],leastgood70[,i],leastgood75[,i],
    leastgood80[,i],leastgood85[,i],leastgood90[,i],leastgood95[,i]),
    .Dim=c(21,21)))
  
  
  
data4<-list(N_st=21,
            N_counts=21,
            para_count = structure(.Data = c(distribsALL[,i],
  randtest01[,i],randtest05[,i],randtest10[,i],randtest15[,i],
  randtest20[,i],randtest25[,i],randtest30[,i],randtest35[,i],
  randtest40[,i],randtest45[,i],randtest50[,i],randtest55[,i],
  randtest60[,i],randtest65[,i],randtest70[,i],randtest75[,i],
  randtest80[,i],randtest85[,i],randtest90[,i],randtest95[,i]),
  .Dim=c(21,21)),###[N_counts,N_st]
  prev = structure(.Data =c(prevdist[,i],                            
                            randprev01[,i],randprev05[,i],randprev10[,i],randprev15[,i],
                            randprev20[,i],randprev25[,i],randprev30[,i],randprev35[,i],
                            randprev40[,i],randprev45[,i],randprev50[,i],randprev55[,i],
                            randprev60[,i],randprev65[,i],randprev70[,i],randprev75[,i],
                            randprev80[,i],randprev85[,i],randprev90[,i],randprev95[,i]),
                   .Dim=c(21,21)))
  fit1 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\modelA2c.stan", data=data2,
               iter=1000, chains=2)
fit2 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080\\modelA2.stan", data=data3,
             iter=1000, chains=2)
fit3 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080\\modelA2.stan", data=data4,
             iter=1000, chains=2)

  #print(fit1)
  #data$Label
  #proportions<-c(0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
params = extract(fit1);names(params)

##for (j in 1:21){
##  dataoutmean_i[j,i]<-c(mean(params$theta[,j])
##}
  dataoutmean_i[,i]<-c(mean(params$theta[501:1000,1]),mean(params$theta[501:1000,2]),mean(params$theta[501:1000,3]),mean(params$theta[501:1000,4]),
                       mean(params$theta[501:1000,5]),mean(params$theta[501:1000,6]),mean(params$theta[501:1000,7]),mean(params$theta[501:1000,8]),
                       mean(params$theta[501:1000,9]),mean(params$theta[501:1000,10]),mean(params$theta[501:1000,11]),mean(params$theta[501:1000,12]),
                       mean(params$theta[501:1000,13]),mean(params$theta[501:1000,14]),mean(params$theta[501:1000,15]),mean(params$theta[501:1000,16]),
                       mean(params$theta[501:1000,17]),mean(params$theta[501:1000,18]),mean(params$theta[501:1000,19]),mean(params$theta[501:1000,20]),mean(params$theta[501:1000,21]))
  dataoutvar_i[,i]<-c(var(params$theta[501:1000,1]),var(params$theta[501:1000,2]),var(params$theta[501:1000,3]),var(params$theta[501:1000,4]),
                     var(params$theta[501:1000,5]),var(params$theta[501:1000,6]),var(params$theta[501:1000,7]),var(params$theta[501:1000,8]),
                     var(params$theta[501:1000,9]),var(params$theta[501:1000,10]),var(params$theta[501:1000,11]),var(params$theta[501:1000,12]),
                     var(params$theta[501:1000,13]),var(params$theta[501:1000,14]),var(params$theta[501:1000,15]),var(params$theta[501:1000,16]),
                     var(params$theta[501:1000,17]),var(params$theta[501:1000,18]),var(params$theta[501:1000,19]),var(params$theta[501:1000,20]),var(params$theta[501:1000,21]))#
  dataout95upper_i[,i]<-as.numeric(c(quantile(params$theta[,1],0.975),quantile(params$theta[,2],0.975),quantile(params$theta[,3],0.975),quantile(params$theta[,4],0.975),
                      quantile(params$theta[,5],0.975),quantile(params$theta[,6],0.975),quantile(params$theta[,7],0.975),quantile(params$theta[,8],0.975),
                      quantile(params$theta[,9],0.975),quantile(params$theta[,10],0.975),quantile(params$theta[,11],0.975),quantile(params$theta[,12],0.975),
                      quantile(params$theta[,13],0.975),quantile(params$theta[,14],0.975),quantile(params$theta[,15],0.975),quantile(params$theta[,16],0.975),
                      quantile(params$theta[,17],0.975),quantile(params$theta[,18],0.975),quantile(params$theta[,19],0.975),quantile(params$theta[,20],0.975),quantile(params$theta[,21],0.975)))
  dataout95lower_i[,i]<-as.numeric(c(quantile(params$theta[,1],0.025),quantile(params$theta[,2],0.025),quantile(params$theta[,3],0.025),quantile(params$theta[,4],0.025),
                      quantile(params$theta[,5],0.025),quantile(params$theta[,6],0.025),quantile(params$theta[,7],0.025),quantile(params$theta[,8],0.025),
                      quantile(params$theta[,9],0.025),quantile(params$theta[,10],0.025),quantile(params$theta[,11],0.025),quantile(params$theta[,12],0.025),
                      quantile(params$theta[,13],0.025),quantile(params$theta[,14],0.025),quantile(params$theta[,15],0.025),quantile(params$theta[,16],0.025),
                      quantile(params$theta[,17],0.025),quantile(params$theta[,18],0.025),quantile(params$theta[,19],0.025),quantile(params$theta[,20],0.025),quantile(params$theta[,21],0.025)))

params2 = extract(fit2);names(params2)
data3outmean_i[,i]<-c(mean(params2$theta[,1]),mean(params2$theta[,2]),mean(params2$theta[,3]),mean(params2$theta[,4]),
                     mean(params2$theta[,5]),mean(params2$theta[,6]),mean(params2$theta[,7]),mean(params2$theta[,8]),
                     mean(params2$theta[,9]),mean(params2$theta[,10]),mean(params2$theta[,11]),mean(params2$theta[,12]),
                     mean(params2$theta[,13]),mean(params2$theta[,14]),mean(params2$theta[,15]),mean(params2$theta[,16]),
                     mean(params2$theta[,17]),mean(params2$theta[,18]),mean(params2$theta[,19]),mean(params2$theta[,20]),mean(params2$theta[,21]))
data3out95upper_i[,i]<-as.numeric(c(quantile(params2$theta[,1],0.975),quantile(params2$theta[,2],0.975),quantile(params2$theta[,3],0.975),quantile(params2$theta[,4],0.975),
                                   quantile(params2$theta[,5],0.975),quantile(params2$theta[,6],0.975),quantile(params2$theta[,7],0.975),quantile(params2$theta[,8],0.975),
                                   quantile(params2$theta[,9],0.975),quantile(params2$theta[,10],0.975),quantile(params2$theta[,11],0.975),quantile(params2$theta[,12],0.975),
                                   quantile(params2$theta[,13],0.975),quantile(params2$theta[,14],0.975),quantile(params2$theta[,15],0.975),quantile(params2$theta[,16],0.975),
                                   quantile(params2$theta[,17],0.975),quantile(params2$theta[,18],0.975),quantile(params2$theta[,19],0.975),quantile(params2$theta[,20],0.975),quantile(params2$theta[,21],0.975)))
data3out95lower_i[,i]<-as.numeric(c(quantile(params2$theta[,1],0.025),quantile(params2$theta[,2],0.025),quantile(params2$theta[,3],0.025),quantile(params2$theta[,4],0.025),
                                   quantile(params2$theta[,5],0.025),quantile(params2$theta[,6],0.025),quantile(params2$theta[,7],0.025),quantile(params2$theta[,8],0.025),
                                   quantile(params2$theta[,9],0.025),quantile(params2$theta[,10],0.025),quantile(params2$theta[,11],0.025),quantile(params2$theta[,12],0.025),
                                   quantile(params2$theta[,13],0.025),quantile(params2$theta[,14],0.025),quantile(params2$theta[,15],0.025),quantile(params2$theta[,16],0.025),
                                   quantile(params2$theta[,17],0.025),quantile(params2$theta[,18],0.025),quantile(params2$theta[,19],0.025),quantile(params2$theta[,20],0.025),quantile(params2$theta[,21],0.025)))

params3 = extract(fit3);names(params3)
data4outmean_i[,i]<-c(mean(params3$theta[,1]),mean(params3$theta[,2]),mean(params3$theta[,3]),mean(params3$theta[,4]),
                     mean(params3$theta[,5]),mean(params3$theta[,6]),mean(params3$theta[,7]),mean(params3$theta[,8]),
                     mean(params3$theta[,9]),mean(params3$theta[,10]),mean(params3$theta[,11]),mean(params3$theta[,12]),
                     mean(params3$theta[,13]),mean(params3$theta[,14]),mean(params3$theta[,15]),mean(params3$theta[,16]),
                     mean(params3$theta[,17]),mean(params3$theta[,18]),mean(params3$theta[,19]),mean(params3$theta[,20]),mean(params3$theta[,21]))
data4out95upper_i[,i]<-as.numeric(c(quantile(params3$theta[,1],0.975),quantile(params3$theta[,2],0.975),quantile(params3$theta[,3],0.975),quantile(params3$theta[,4],0.975),
                                   quantile(params3$theta[,5],0.975),quantile(params3$theta[,6],0.975),quantile(params3$theta[,7],0.975),quantile(params3$theta[,8],0.975),
                                   quantile(params3$theta[,9],0.975),quantile(params3$theta[,10],0.975),quantile(params3$theta[,11],0.975),quantile(params3$theta[,12],0.975),
                                   quantile(params3$theta[,13],0.975),quantile(params3$theta[,14],0.975),quantile(params3$theta[,15],0.975),quantile(params3$theta[,16],0.975),
                                   quantile(params3$theta[,17],0.975),quantile(params3$theta[,18],0.975),quantile(params3$theta[,19],0.975),quantile(params3$theta[,20],0.975),quantile(params3$theta[,21],0.975)))
data4out95lower_i[,i]<-as.numeric(c(quantile(params3$theta[,1],0.025),quantile(params3$theta[,2],0.025),quantile(params3$theta[,3],0.025),quantile(params3$theta[,4],0.025),
                                   quantile(params3$theta[,5],0.025),quantile(params3$theta[,6],0.025),quantile(params3$theta[,7],0.025),quantile(params3$theta[,8],0.025),
                                   quantile(params3$theta[,9],0.025),quantile(params3$theta[,10],0.025),quantile(params3$theta[,11],0.025),quantile(params3$theta[,12],0.025),
                                   quantile(params3$theta[,13],0.025),quantile(params3$theta[,14],0.025),quantile(params3$theta[,15],0.025),quantile(params3$theta[,16],0.025),
                                   quantile(params3$theta[,17],0.025),quantile(params3$theta[,18],0.025),quantile(params3$theta[,19],0.025),quantile(params3$theta[,20],0.025),quantile(params3$theta[,21],0.025)))

}
  dataoutupper_i<-  dataoutmean_i+dataoutvar_i
  dataoutlower_i<-  dataoutmean_i-dataoutvar_i

#write.csv(dataoutmean_i,"C:\\Users\\Ellie\\Documents\\2080\\model_outputMEAN_i2.csv")
#write.csv(dataout95upper_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95upper_i2.csv")
#write.csv(dataout95lower_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95lower_i2.csv")

#write.csv(data3outmean_i,"C:\\Users\\Ellie\\Documents\\2080\\least infected data_mean.csv")
#write.csv(data3out95upper_i,"C:\\Users\\Ellie\\Documents\\2080\\least infected data_upper.csv")
#write.csv(data3out95lower_i,"C:\\Users\\Ellie\\Documents\\2080\\least infected data_lower.csv")

#write.csv(data4outmean_i,"C:\\Users\\Ellie\\Documents\\2080\\random data_mean.csv")
#write.csv(data4out95upper_i,"C:\\Users\\Ellie\\Documents\\2080\\random data_upper.csv")
#write.csv(data4out95lower_i,"C:\\Users\\Ellie\\Documents\\2080\\random data_lower.csv")

par(mfrow=c(1,1))
dataoutmean<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\topmost infected data_mean.csv")
dataout95upper<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\topmost infected data_upper.csv")
dataout95lower<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\topmost infected data_lower.csv")

data3outmean<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\least infected data_mean.csv")
data3out95upper<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\least infected data_upper.csv")
data3out95lower<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\least infected data_lower.csv")

data4outmean<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\random data_mean.csv")
data4out95upper<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\random data_upper.csv")
data4out95lower<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\random data_lower.csv")

###check that the first column is DATA else delete ## rename columns
colnames(dataoutmean)<-data$Label[1:206]
colnames(data3outmean)<-colnames(data4outmean)<-c("count",data$Label[1:206])

reductionacheivedt20<-expand.grid(seq(1,20))##reduction acheived most infected hosts
for(i in 1:206){
  for (j in 2:21){
    reductionacheivedt20[j-1,i]<-(dataoutmean[1,i]-dataoutmean[j,i])/dataoutmean[1,i]
  }
};colnames(reductionacheivedt20)<-data$Label[1:206] 
#reductionacheivedt20[21,]<-reductionacheivedw20[21,]<-1

reductionacheivedw20<-expand.grid(seq(1,20))##reduction for least infected hosts
for(i in 2:207){
  for (j in 2:21){
    reductionacheivedw20[j-1,i]<-(data3outmean[1,i]-data3outmean[j,i])/data3outmean[1,i]
  }
};colnames(reductionacheivedw20)<-c("count",data$Label[1:206])

  for (j in 2:207){
    for(i in 1:20){
reductionacheivedw20[i,j]<-ifelse(reductionacheivedw20[i,j]<0,0,reductionacheivedw20[i,j])
}}

reductionacheivedr20<-expand.grid(seq(1,20))##reduction for least infected hosts
for(i in 2:207){
  for (j in 2:21){
    reductionacheivedr20[j-1,i]<-(data4outmean[1,i]-data4outmean[j,i])/data4outmean[1,i]
  }
};colnames(reductionacheivedr20)<-c("count",data$Label[1:206]) 

for (j in 2:207){
  for(i in 1:20){
    reductionacheivedr20[i,j]<-ifelse(reductionacheivedr20[i,j]<0,0,reductionacheivedr20[i,j])
  }}
#################
## Figure 4
par(new=TRUE)
proportion<-c(0.01,seq(0.05,0.95,0.05))
par(mar=c(10,25,8,5))
plot(reductionacheivedt20[,1]~proportion,pch="",ylim=c(0,1),xlim=c(0,1),par(las=1),
     xlab="Proportion of the host population treated",cex.lab=1.1,bty="n",yaxt="n",
     ylab=expression(paste("Effective reduction in  ", theta)))
axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),par(las=2))

#for (i in 1:206){
#lines(reductionacheivedt20[,i]~proportion,lty=3,col="grey65")
#lines(reductionacheivedw20[,i]~proportion,lty=3,col="grey80")
#lines(reductionacheivedr20[,i]~proportion,lty=3,col="grey60")
#}
meansdat<-means3dat<-means4dat<-maxdat<-max3dat<-max4dat<-mindat<-min3dat<-min4dat<-sddat<-sd3dat<-sd4dat<-numeric(20)
for (i in 1:20){
  meansdat[i]<-(sum(reductionacheivedt20[i,])/206)
  means3dat[i]<-(sum(reductionacheivedw20[i,2:207],na.rm=TRUE)/204)
  means4dat[i]<-(sum(reductionacheivedr20[i,2:207],na.rm=TRUE)/204)
  
  maxdat[i]<-max(reductionacheivedt20[i,])
  max3dat[i]<-max(reductionacheivedw20[i,2:207],na.rm=TRUE)
  max4dat[i]<-max(reductionacheivedr20[i,2:207],na.rm=TRUE)
  
  mindat[i]<-min(reductionacheivedt20[i,])
  min3dat[1]<-min(reductionacheivedw20[i,2:207],na.rm=TRUE)
  min4dat[i]<-min(reductionacheivedr20[i,2:207],na.rm=TRUE)
  
  sddat[i]<-sd(reductionacheivedt20[i,],na.rm=TRUE)
  sd3dat[i]<-sd(reductionacheivedw20[i,2:207],na.rm=TRUE)
  sd4dat[i]<-sd(reductionacheivedr20[i,2:207],na.rm=TRUE)
  
}
lines(c(0,meansdat)~c(0,proportion),col="black",lty=1,lwd=2)
lines(c(0,means3dat)~c(0,proportion),col="grey50",lty=1,lwd=2)
lines(c(0,means4dat)~c(0,proportion),col="grey50",lty=1,lwd=2)

sdplus<-meansdat+sddat;sdminus<-meansdat-sddat
sdplus<-ifelse(sdplus<1,sdplus,1);sdminus<-ifelse(sdminus<0,0,sdminus)
sdplus3<-means3dat+sd3dat;sdminus3<-means3dat-sd3dat
sdplus3<-ifelse(sdplus3<1,sdplus3,1);sdminus3<-ifelse(sdminus3<0,0,sdminus3)
sdplus4<-means4dat+sd4dat;sdminus4<-means4dat-sd4dat
sdplus4<-ifelse(sdplus4<1,sdplus4,1);sdminus4<-ifelse(sdminus4<0,0,sdminus4)

polygon(c(proportion, rev(proportion)),c(sdplus,rev(sdminus)),border=NA, col=transp("darkseagreen1",alpha=0.3))
polygon(c(proportion, rev(proportion)),c(sdplus3,rev(sdminus3)),border=NA, col=transp("darkseagreen4",alpha=0.2))
polygon(c(proportion, rev(proportion)),c(sdplus4,rev(sdminus4)),border=NA, col=transp("darkseagreen3",alpha=0.2))

segments(0.2, y0=-0.2, 0.2, y1 = meansdat[5],lty = 4, lwd = 1)
segments(-0.2, y0=meansdat[5], 0.2, y1 = meansdat[5],lty = 4, lwd = 1)
segments(-0.2, y0=means3dat[5], 0.2, y1 = means3dat[5],lty = 4, lwd = 1)
segments(-0.2, y0=means4dat[5], 0.2, y1 = means4dat[5],lty = 4, lwd = 1)



countzeros<-as.numeric(reductionacheivedw20[5,],na.action=na.exclude)
countzeros<-c(countzeros[2],countzeros[4:80],countzeros[82:207])
sum(ifelse(countzeros==0,1,0))
#write.csv(dataoutmean_i,"C:\\Users\\Ellie\\Documents\\2080\\model_outputMEAN_i2.csv")
#write.csv(dataout95upper_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95upper_i2.csv")
#write.csv(dataout95lower_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95lower_i2.csv")
treatments<-rep(c("0","1",seq(from=5,to=95,by=5)),each=21)


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
plot((1-data3$Gini.Co.efficient),data$T10,pch=20,col="grey15",
     ylab="Proportion of parasites in x% of hosts",
     xlab="Gini Co-efficient",cex.lab=1.1,xlim=c(0,1))
points(1-data3$Gini.Co.efficient,data$T20,col="grey30",pch=20,cex=2)
points(1-data3$Gini.Co.efficient,data$T60,pch=20,col="grey50")
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



############################
## Figure 3
###############################

data<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\Exploring R0.csv",header=TRUE)
data[1:10,];summary(data);dim(data)
##Trying for study 1 (Label 36)
pointsdat<-matrix(nrow=max(data$N),ncol=206)
teststore<-matrix(nrow=max(data$N),ncol=206)

proportions<-seq(0,1,length=data[1,4])
for (j in 1:206){
  
  d1<-sort(rnegbin(data[j,4],data[j,6],data[j,5]),decreasing=TRUE)
  dat<-expand.grid(d1)
  test<-R0<-prev1<-kx<-dat1<-numeric(length(d1))
  
  #R0[1]<- (data[1,5] * ( pr1 ^ (-1/data[1,5]))) - data[1,5]
  
  for(i in 1:nrow(dat)){
    dat[,i+1]<-ifelse(dat[,i]==max(dat[,i]),0,dat[,i])
  }
  
  for (i in 1:length(prev1)){
    prev1[i]<-sum(ifelse(dat[,i]==0,0,1))/length(dat[,i])
  }
  
  kx[1]<-data[j,5]
  for (i in 2:length(kx)){
    kx[i]<-(mean(dat[,i])^2-(var(dat[,i])/length(dat[,1])))/(var(dat[,i])-mean(dat[,i]))
    kx<-ifelse(kx<0,NA,kx)
  }
  
  
  
  for (i in 1:length(R0)){
    R0[i]<- (kx[i] * ( (1-prev1[i]) ^ (-1/kx[i]))) - kx[i]##1-prev as it is those who are not contributing to transmission
    R0[is.na(R0)] <- 0
    #R0[is.infinite(R0)]<- 10
  }
  
  for (i in 2:length(R0)){
    test[1]<-0
    test[i] <-  (R0[1]-R0[i])/R0[1]
  }
  teststore[,j]<-c(test,rep(99999,(2611-length(d1))))
  proports<-seq(0,1,length=data[j,4])
  a<-max(data$N-length(data[j,4]))
  pointsdat[,j]<-c(proports,rep(99999,(2611-length(proports))))

  log.binom<-function(p.vec){
    
    a<-p.vec[1]
    b<-p.vec[2]
    
    pred1a<- ((exp(a + b * proports)) / (1 + exp(a + b * proports)) ) 
    prev1<-test
    
    loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
    -sum(loglik1a,  na.rm=T)
  }
  n.param<-2
  logmod<-optim(c(0,5),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,100))
  logmod
  nc<-seq(0,1,0.001)
  pred2<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) ) 
  lines(nc,pred2,lwd=3,lty=1,col=transp("aquamarine4"))
  
}


t01R0<-t05R0<-t10R0<-t15R0<-t20R0<-t25R0<-t30R0<-t35R0<-
  t40R0<-t45R0<-t50R0<-
  t55R0<-t60R0<-t65R0<-t70R0<-t75R0<-t80R0<-t85R0<-
  t90R0<-t95R0<-t99R0<-numeric(206)
for(i in 1:206){
  a<-(teststore[,i][round(pointsdat[,i],2)==0.01]);t01R0[i]<-a[1]
  b<-(teststore[,i][round(pointsdat[,i],2)==0.05]);t05R0[i]<-b[1]
  d<-(teststore[,i][round(pointsdat[,i],2)==0.1]);t10R0[i]<-d[1]
  ee<-(teststore[,i][round(pointsdat[,i],2)==0.15]);t15R0[i]<-ee[1]
  g<-(teststore[,i][round(pointsdat[,i],2)==0.2]);t20R0[i]<-g[1]
  h<-(teststore[,i][round(pointsdat[,i],2)==0.25]);t25R0[i]<-h[1]
  m<-(teststore[,i][round(pointsdat[,i],2)==0.30]);t30R0[i]<-m[1]
  m1<-(teststore[,i][round(pointsdat[,i],2)==0.35]);t35R0[i]<-m1[1]
  m2<-(teststore[,i][round(pointsdat[,i],2)==0.40]);t40R0[i]<-m2[1]
  m3<-(teststore[,i][round(pointsdat[,i],2)==0.45]);t45R0[i]<-m3[1]
  bg<-(teststore[,i][round(pointsdat[,i],2)==0.50]);t50R0[i]<-bg[1]
  bg1<-(teststore[,i][round(pointsdat[,i],2)==0.55]);t55R0[i]<-bg1[1]
  bg2<-(teststore[,i][round(pointsdat[,i],2)==0.60]);t60R0[i]<-bg2[1]
  bg3<-(teststore[,i][round(pointsdat[,i],2)==0.65]);t65R0[i]<-bg3[1]
  bg4<-(teststore[,i][round(pointsdat[,i],2)==0.70]);t70R0[i]<-bg4[1]
  bg5<-(teststore[,i][round(pointsdat[,i],2)==0.75]);t75R0[i]<-bg5[1]
  bm<-(teststore[,i][round(pointsdat[,i],2)==0.80]);t80R0[i]<-bm[1]
  bm1<-(teststore[,i][round(pointsdat[,i],2)==0.85]);t85R0[i]<-bm1[1]
  bm2<-(teststore[,i][round(pointsdat[,i],2)==0.90]);t90R0[i]<-bm2[1]
  bm3<-(teststore[,i][round(pointsdat[,i],2)==0.95]);t95R0[i]<-bm3[1]
  bm4<-(teststore[,i][round(pointsdat[,i],2)==0.99]);t99R0[i]<-bm4[1]
}
t01R0<-ifelse(t01R0<0,0,t01R0)  

test2<-c(0,mean(t01R0,na.rm=TRUE),mean(t05R0,na.rm=TRUE),mean(t10R0,na.rm=TRUE),mean(t15R0,na.rm=TRUE),mean(t20R0,na.rm=TRUE),
         mean(t25R0,na.rm=TRUE),mean(t30R0,na.rm=TRUE),mean(t35R0,na.rm=TRUE),mean(t40R0,na.rm=TRUE),
         mean(t45R0,na.rm=TRUE),mean(t50R0,na.rm=TRUE),mean(t55R0,na.rm=TRUE),mean(t60R0,na.rm=TRUE),
         mean(t65R0,na.rm=TRUE),mean(t70R0,na.rm=TRUE),mean(t75R0,na.rm=TRUE),mean(t80R0,na.rm=TRUE),
         mean(t85R0,na.rm=TRUE),mean(t90R0,na.rm=TRUE),mean(t95R0,na.rm=TRUE),mean(t99R0,na.rm=TRUE))
test2up<-c(0,quantile(t01R0,0.975,na.rm=TRUE),quantile(t05R0,0.975,na.rm=TRUE),quantile(t10R0,0.975,na.rm=TRUE),quantile(t15R0,0.975,na.rm=TRUE),quantile(t20R0,0.975,na.rm=TRUE),
           quantile(t25R0,0.975,na.rm=TRUE),quantile(t30R0,0.975,na.rm=TRUE),quantile(t35R0,0.975,na.rm=TRUE),
           quantile(t40R0,0.975,na.rm=TRUE),quantile(t45R0,0.975,na.rm=TRUE),quantile(t50R0,0.975,na.rm=TRUE),
           quantile(t55R0,0.975,na.rm=TRUE),quantile(t60R0,0.975,na.rm=TRUE),quantile(t65R0,0.975,na.rm=TRUE),
           quantile(t70R0,0.975,na.rm=TRUE),quantile(t75R0,0.975,na.rm=TRUE),quantile(t80R0,0.975,na.rm=TRUE),
           quantile(t85R0,0.975,na.rm=TRUE),quantile(t90R0,0.975,na.rm=TRUE),quantile(t95R0,0.975,na.rm=TRUE),quantile(t99R0,0.975,na.rm=TRUE))
test2low<-c(0,quantile(t01R0,0.025,na.rm=TRUE),quantile(t05R0,0.025,na.rm=TRUE),quantile(t10R0,0.025,na.rm=TRUE),quantile(t15R0,0.025,na.rm=TRUE),quantile(t20R0,0.025,na.rm=TRUE),
            quantile(t25R0,0.025,na.rm=TRUE),quantile(t30R0,0.025,na.rm=TRUE),quantile(t35R0,0.025,na.rm=TRUE),
            quantile(t40R0,0.025,na.rm=TRUE),quantile(t45R0,0.025,na.rm=TRUE),quantile(t50R0,0.025,na.rm=TRUE),
            quantile(t55R0,0.025,na.rm=TRUE),quantile(t60R0,0.025,na.rm=TRUE),quantile(t65R0,0.025,na.rm=TRUE),
            quantile(t70R0,0.025,na.rm=TRUE),quantile(t75R0,0.025,na.rm=TRUE),quantile(t80R0,0.025,na.rm=TRUE),
            quantile(t85R0,0.025,na.rm=TRUE),quantile(t90R0,0.025,na.rm=TRUE),quantile(t95R0,0.025,na.rm=TRUE),quantile(t99R0,0.025,na.rm=TRUE))
proports2<-c(0,0.01,seq(0.05,0.95,0.05),0.99)

log.binom<-function(p.vec){
  
  a<-p.vec[1]
  b<-p.vec[2]
  
  pred1a<- ((exp(a + b * proports2)) / (1 + exp(a + b * proports2)) ) 
  prev1<-test2
  
  loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
  -sum(loglik1a,  na.rm=T)
}
n.param<-2
logmod<-optim(c(0,5),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,100))
logmod
nc<-seq(0,1,0.001)
pred2<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) ) 


log.binom<-function(p.vec){
  
  a<-p.vec[1]
  b<-p.vec[2]
  
  pred1a<- ((exp(a + b * proports2)) / (1 + exp(a + b * proports2)) ) 
  prev1<-test2up
  
  loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
  -sum(loglik1a,  na.rm=T)
}
n.param<-2
logmod<-optim(c(0,5),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,100))
logmod
nc<-seq(0,1,0.001)
pred2u<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) ) 
pred2u
#lines(nc,pred2u,lwd=3,lty=2,col=transp("aquamarine4"))

log.binom<-function(p.vec){
  
  a<-p.vec[1]
  b<-p.vec[2]
  
  pred1a<- ((exp(a + b * proports2)) / (1 + exp(a + b * proports2)) ) 
  prev1<-test2low
  
  loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
  -sum(loglik1a,  na.rm=T)
}
n.param<-2
logmod<-optim(c(0,5),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,100))
logmod
nc<-seq(0,1,0.001)
pred2l<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) ) 

par(mar=c(5,6,2,2))
plot(nc,pred2,pty="n",pch="",yaxt="n",xaxt="n",ylab="",xlab="")

polygon(c(nc, rev(nc)),c(pred2u,rev(pred2l)),border=NA, col=transp("darkseagreen1",alpha=0.3))
lines(nc,pred2,lwd=3,lty=1,col=transp("aquamarine4"))

par(new=TRUE)
par(bty="n")
par(las=2)
par(mar=c(5,5,2,2))
boxplot(rep(0,206),t01R0,t05R0,t10R0,t15R0,t20R0,t25R0,t30R0,t35R0,t40R0,t45R0,t50R0,
        t55R0,t60R0,t65R0,t70R0,t75R0,t80R0,t85R0,t90R0,t95R0,t99R0,
        na.rm=TRUE,
        ylim=c(0,1),col=transp("aquamarine4",alpha=0.2),cex.lab=1.1,at=seq(-1,20,1),
        ylab=expression(paste("Relative reduction in   ",   R[0])),xaxt="n",
        xlab="Proportion of hosts treated from most to least infected")
par(las=1)
proports3<-seq(0.0,1,0.1)
axis(1,at=seq(-1,20,length=11),labels=proports3,cex.lab=1.5)























##############
## Other info....
####
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

xv<-data$Label[1:206]

D1<-data.frame(xv,rem20red,rem20redU,rem20redL,remleast20red,remleast20redU,remleast20redL,
               remrand20red,remrand20redU,remrand20redL)
I<-ordered(D1$rem20red) ## order the data by the ones in top 20

par(mfrow=c(1,1))
par(mar=c(5,5,2,2))

#xv<-data$prevalence[1:206]/data$N[1:206]
plot(D1$rem20red~D1$xv,ylab=expression(paste("Proportionate reduction in  ",theta)),
     ylim=c(0,1),
     xlab="Studies",
     pch="",cex.lab=1.2)

for (i in 1:206){
  segments(D1$xv[i], D1$rem20redU[i], x1 = D1$xv[i], y1 = D1$rem20redL[i],
           col  = terrain.colors(10,alpha = 0.7), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(D1$xv[i], D1$remleast20redU[i], x1 = D1$xv[i], y1 = D1$remleast20redL[i],
           col  = terrain.colors(10,alpha = 0.3), lty = 1, lwd = 5)
}

for (i in 1:206){
  segments(D1$xv[i], D1$remrand20redU[i], x1 = D1$xv[i], y1 = D1$remrand20redL[i],
           col  = terrain.colors(10,alpha = 0.1), lty = 1, lwd = 5)
}

points(D1$rem20red~D1$xv,col="lightblue",pch=20)
points(D1$remleast20red~D1$xv,col="black",pch=20)
points(D1$remrand20red~D1$xv,col="blue",pch=20)

abline(h=mean(D1$rem20red,na.rm=TRUE),col="lightblue")
abline(h=mean(D1$remleast20red,na.rm=TRUE),col="black",lty=2)
abline(h=mean(D1$remrand20red,na.rm=TRUE),col="blue",lty=2)


legend(250,0.4,legend=c("Targeted most infected","Least infected","Randomly infected"),
       col=c("lightblue","blue","black"),
       pch=20,bty="n")

sum(ifelse(rem20red>0.95,1,0))/206 ##A 95% or better reduction in transmission probability can be acheived in 32% of the studies
sum(ifelse(rem50red>0.95,1,0))/206 