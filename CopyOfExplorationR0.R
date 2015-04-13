
library(rstan)
library(MASS)
library(boot)
library(coda)
library(R2OpenBUGS)

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

parasitesTreated<-as.numeric(c(data[1,14:34]))
plot(parasitesTreated~prevalence,
     ylim=c(0,1),xlim=c(0,1),
     ylab="Proportion of parasites removed per hosts treated",xlab="Proportion of hosts treated")
treated<-matrix(nrow=length(prevalence),ncol=length(data$N),data=NA)
for (i in 1:length(data$N)){
  treated[,i]<-as.numeric(c(data[i,14:34]))
  lines(treated[,i]~prevalence)
}

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
disttemp<-prevdist;disttemp[1,]<-0;prev01<-disttemp
disttemp[1:2,]<-0;prev05<-disttemp
disttemp[1:3,]<-0;prev10<-disttemp
disttemp[1:4,]<-0;prev15<-disttemp
disttemp[1:5,]<-0;prev20<-disttemp
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

###Take for example the first case


data2<-list(N_st=21,
            N_counts=21,
            para_count = structure(.Data = c(distribsALL[,11],distrib01[,11],distrib05[,11],distrib10[,11],
                                             distrib15[,11],distrib20[,11],distrib25[,11],distrib30[,11],
                                             distrib35[,11],distrib40[,11],distrib45[,11],distrib50[,11],
                                             distrib55[,11],distrib60[,11],distrib65[,11],distrib70[,11],
                                             distrib75[,11],distrib80[,11],distrib85[,11],distrib90[,11],
                                             distrib95[,11]),
                                  .Dim=c(21,21)),###[N_counts,N_st]
            prev = structure(.Data =c(prevdist[,11],prev01[,11],prev05[,11],prev10[,11],prev15[,11],prev20[,11],
                                      prev25[,11],prev30[,11],prev35[,11],prev40[,11],prev45[,11],prev50[,11],
                                      prev55[,11],prev60[,11],prev65[,11],prev70[,11],prev75[,11],prev80[,11],
                                      prev85[,11],prev90[,11],prev95[,11]),
                                  .Dim=c(21,21)))

fit1 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080\\modelA2.stan", data=data2,
             iter=1000, chains=2)

print(fit1)
data$Label
proportions<-c(0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
params = extract(fit1);names(params)
theta47<-c(mean(params$theta[,1]),mean(params$theta[,2]),mean(params$theta[,3]),mean(params$theta[,4]),
           mean(params$theta[,5]),mean(params$theta[,6]),mean(params$theta[,7]),mean(params$theta[,8]),
           mean(params$theta[,9]),mean(params$theta[,10]),mean(params$theta[,11]),mean(params$theta[,12]),
           mean(params$theta[,13]),mean(params$theta[,14]),mean(params$theta[,15]),mean(params$theta[,16]),
           mean(params$theta[,17]),mean(params$theta[,18]),mean(params$theta[,19]),mean(params$theta[,20]),
           mean(params$theta[,21]))
theta47Upper<-as.numeric(c(quantile(params$theta[,1],0.975),quantile(params$theta[,2],0.975),quantile(params$theta[,3],0.975),quantile(params$theta[,4],0.975),
           quantile(params$theta[,5],0.975),quantile(params$theta[,6],0.975),quantile(params$theta[,7],0.975),quantile(params$theta[,8],0.975),
           quantile(params$theta[,9],0.975),quantile(params$theta[,10],0.975),quantile(params$theta[,11],0.975),quantile(params$theta[,12],0.975),
           quantile(params$theta[,13],0.975),quantile(params$theta[,14],0.975),quantile(params$theta[,15],0.975),quantile(params$theta[,16],0.975),
           quantile(params$theta[,17],0.975),quantile(params$theta[,18],0.975),quantile(params$theta[,19],0.975),quantile(params$theta[,20],0.975),
           quantile(params$theta[,21],0.975)))
theta47lower<-as.numeric(c(quantile(params$theta[,1],0.025),quantile(params$theta[,2],0.025),quantile(params$theta[,3],0.025),quantile(params$theta[,4],0.025),
                           quantile(params$theta[,5],0.025),quantile(params$theta[,6],0.025),quantile(params$theta[,7],0.025),quantile(params$theta[,8],0.025),
                           quantile(params$theta[,9],0.025),quantile(params$theta[,10],0.025),quantile(params$theta[,11],0.025),quantile(params$theta[,12],0.025),
                           quantile(params$theta[,13],0.025),quantile(params$theta[,14],0.025),quantile(params$theta[,15],0.025),quantile(params$theta[,16],0.025),
                           quantile(params$theta[,17],0.025),quantile(params$theta[,18],0.025),quantile(params$theta[,19],0.025),quantile(params$theta[,20],0.025),
                           quantile(params$theta[,21],0.025)))


dataoutmean<-expand.grid(c(1:21))
dataoutmean[,1]<-theta36;colnames(dataoutmean)<-sub("Var1","Study36",colnames(dataoutmean))
dataoutmean[,2]<-theta37;colnames(dataoutmean)<-sub("V2","Study37",colnames(dataoutmean))
dataoutmean[,3]<-theta38;colnames(dataoutmean)<-sub("V3","Study38",colnames(dataoutmean))
dataoutmean[,4]<-theta40;colnames(dataoutmean)<-sub("V4","Study40",colnames(dataoutmean))
dataoutmean[,5]<-theta41;colnames(dataoutmean)<-sub("V5","Study41",colnames(dataoutmean))
dataoutmean[,6]<-theta42;colnames(dataoutmean)<-sub("V6","Study42",colnames(dataoutmean))
dataoutmean[,7]<-theta43;colnames(dataoutmean)<-sub("V7","Study43",colnames(dataoutmean))
dataoutmean[,8]<-theta44;colnames(dataoutmean)<-sub("V8","Study44",colnames(dataoutmean))
dataoutmean[,9]<-theta45;colnames(dataoutmean)<-sub("V9","Study45",colnames(dataoutmean))
dataoutmean[,10]<-theta46;colnames(dataoutmean)<-sub("V10","Study46",colnames(dataoutmean))
dataoutmean[,11]<-theta47;colnames(dataoutmean)<-sub("V11","Study47",colnames(dataoutmean))


dataout95upper<-expand.grid(c(1:21))
dataout95upper[,1]<-theta36Upper;colnames(dataout95upper)<-sub("Var1","Study36",colnames(dataout95upper))
dataout95upper[,2]<-theta37Upper;colnames(dataout95upper)<-sub("V2","Study37",colnames(dataout95upper))
dataout95upper[,3]<-theta38Upper;colnames(dataout95upper)<-sub("V3","Study38",colnames(dataout95upper))
dataout95upper[,4]<-theta40Upper;colnames(dataout95upper)<-sub("V4","Study40",colnames(dataout95upper))
dataout95upper[,5]<-theta41Upper;colnames(dataout95upper)<-sub("V5","Study41",colnames(dataout95upper))
dataout95upper[,6]<-theta42Upper;colnames(dataout95upper)<-sub("V6","Study42",colnames(dataout95upper))
dataout95upper[,7]<-theta43Upper;colnames(dataout95upper)<-sub("V7","Study43",colnames(dataout95upper))
dataout95upper[,8]<-theta44Upper;colnames(dataout95upper)<-sub("V8","Study44",colnames(dataout95upper))
dataout95upper[,9]<-theta45Upper;colnames(dataout95upper)<-sub("V9","Study45",colnames(dataout95upper))
dataout95upper[,10]<-theta46Upper;colnames(dataout95upper)<-sub("V10","Study46",colnames(dataout95upper))
dataout95upper[,11]<-theta47Upper;colnames(dataout95upper)<-sub("V11","Study47",colnames(dataout95upper))


dataout95lower<-expand.grid(c(1:21))
dataout95lower[,1]<-theta36lower;colnames(dataout95lower)<-sub("Var1","Study36",colnames(dataout95lower))
dataout95lower[,2]<-theta37lower;colnames(dataout95lower)<-sub("V2","Study37",colnames(dataout95lower))
dataout95lower[,3]<-theta38lower;colnames(dataout95lower)<-sub("V3","Study38",colnames(dataout95lower))
dataout95lower[,4]<-theta40lower;colnames(dataout95lower)<-sub("V4","Study40",colnames(dataout95lower))
dataout95lower[,5]<-theta41lower;colnames(dataout95lower)<-sub("V5","Study41",colnames(dataout95lower))
dataout95lower[,6]<-theta42lower;colnames(dataout95lower)<-sub("V6","Study42",colnames(dataout95lower))
dataout95lower[,7]<-theta43lower;colnames(dataout95lower)<-sub("V7","Study43",colnames(dataout95lower))
dataout95lower[,8]<-theta44lower;colnames(dataout95lower)<-sub("V8","Study44",colnames(dataout95lower))
dataout95lower[,9]<-theta45lower;colnames(dataout95lower)<-sub("V9","Study45",colnames(dataout95lower))
dataout95lower[,10]<-theta46lower;colnames(dataout95lower)<-sub("V10","Study46",colnames(dataout95lower))
dataout95lower[,11]<-theta47lower;colnames(dataout95lower)<-sub("V11","Study47",colnames(dataout95lower))


write.csv(dataoutmean,"C:\\Users\\Ellie\\Documents\\2080\\model_outputMEAN.csv")
write.csv(dataout95upper,"C:\\Users\\Ellie\\Documents\\2080\\model_output95upper.csv")
write.csv(dataout95lower,"C:\\Users\\Ellie\\Documents\\2080\\model_output95lower.csv")

plot(theta36~proportions,ylim=c(0,1),xlim=c(0,1),pch="",
     ylab="Transmission Probability",xlab="Proportion of hosts treated")               
lines(theta36~proportions)
lines(theta37~proportions)
lines(theta38~proportions)
lines(theta40~proportions)
lines(theta41~proportions)
polygon(c(proportions, rev(proportions)),c(theta36Upper,rev(theta36lower)),border=NA, col="aquamarine1")
lines(theta36~proportions)

polygon(c(proportions, rev(proportions)),c(theta37Upper,rev(theta37lower)),border=NA, col="aquamarine1")
lines(theta37~proportions)

polygon(c(proportions, rev(proportions)),c(theta38Upper,rev(theta38lower)),border=NA, col="aquamarine1")
lines(theta38~proportions)

polygon(c(proportions, rev(proportions)),c(theta40Upper,rev(theta40lower)),border=NA, col="aquamarine1")
lines(theta40~proportions)


###################################################################
##################################################################### AUTOMATED
#######################################################################
dataoutmean_i<-expand.grid(c(1:21))
dataout95upper_i<-expand.grid(c(1:21))
dataout95lower_i<-expand.grid(c(1:21))
for (i in 1:206){
  data2<-list(N_st=21,
              N_counts=21,
              para_count = structure(.Data = c(distribsALL[,i],distrib01[,i],distrib05[,i],distrib10[,i],
                                               distrib15[,i],distrib20[,i],distrib25[,i],distrib30[,i],
                                               distrib35[,i],distrib40[,i],distrib45[,i],distrib50[,i],
                                               distrib55[,i],distrib60[,i],distrib65[,i],distrib70[,i],
                                               distrib75[,i],distrib80[,i],distrib85[,i],distrib90[,i],
                                               distrib95[,i]),
                                     .Dim=c(21,21)),###[N_counts,N_st]
              prev = structure(.Data =c(prevdist[,i],prev01[,i],prev05[,i],prev10[,i],prev15[,i],prev20[,i],
                                        prev25[,i],prev30[,i],prev35[,i],prev40[,i],prev45[,i],prev50[,i],
                                        prev55[,i],prev60[,i],prev65[,i],prev70[,i],prev75[,i],prev80[,i],
                                        prev85[,i],prev90[,i],prev95[,i]),
                               .Dim=c(21,21)))

fit1 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080\\modelA2.stan", data=data2,
             iter=1000, chains=2)

#print(fit1)
#data$Label
#proportions<-c(0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
params = extract(fit1);names(params)
dataoutmean_i[,i]<-c(mean(params$theta[,1]),mean(params$theta[,2]),mean(params$theta[,3]),mean(params$theta[,4]),
           mean(params$theta[,5]),mean(params$theta[,6]),mean(params$theta[,7]),mean(params$theta[,8]),
           mean(params$theta[,9]),mean(params$theta[,10]),mean(params$theta[,11]),mean(params$theta[,12]),
           mean(params$theta[,13]),mean(params$theta[,14]),mean(params$theta[,15]),mean(params$theta[,16]),
           mean(params$theta[,17]),mean(params$theta[,18]),mean(params$theta[,19]),mean(params$theta[,20]),
           mean(params$theta[,21]))
dataout95upper_i[,i]<-as.numeric(c(quantile(params$theta[,1],0.975),quantile(params$theta[,2],0.975),quantile(params$theta[,3],0.975),quantile(params$theta[,4],0.975),
                           quantile(params$theta[,5],0.975),quantile(params$theta[,6],0.975),quantile(params$theta[,7],0.975),quantile(params$theta[,8],0.975),
                           quantile(params$theta[,9],0.975),quantile(params$theta[,10],0.975),quantile(params$theta[,11],0.975),quantile(params$theta[,12],0.975),
                           quantile(params$theta[,13],0.975),quantile(params$theta[,14],0.975),quantile(params$theta[,15],0.975),quantile(params$theta[,16],0.975),
                           quantile(params$theta[,17],0.975),quantile(params$theta[,18],0.975),quantile(params$theta[,19],0.975),quantile(params$theta[,20],0.975),
                           quantile(params$theta[,21],0.975)))
dataout95lower_i[,i]<-as.numeric(c(quantile(params$theta[,1],0.025),quantile(params$theta[,2],0.025),quantile(params$theta[,3],0.025),quantile(params$theta[,4],0.025),
                           quantile(params$theta[,5],0.025),quantile(params$theta[,6],0.025),quantile(params$theta[,7],0.025),quantile(params$theta[,8],0.025),
                           quantile(params$theta[,9],0.025),quantile(params$theta[,10],0.025),quantile(params$theta[,11],0.025),quantile(params$theta[,12],0.025),
                           quantile(params$theta[,13],0.025),quantile(params$theta[,14],0.025),quantile(params$theta[,15],0.025),quantile(params$theta[,16],0.025),
                           quantile(params$theta[,17],0.025),quantile(params$theta[,18],0.025),quantile(params$theta[,19],0.025),quantile(params$theta[,20],0.025),
                           quantile(params$theta[,21],0.025)))

}

write.csv(dataoutmean_i,"C:\\Users\\Ellie\\Documents\\2080\\model_outputMEAN_i.csv")
write.csv(dataout95upper_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95upper_i.csv")
write.csv(dataout95lower_i,"C:\\Users\\Ellie\\Documents\\2080\\model_output95lower_i.csv")

plot(theta36~proportions,ylim=c(0,1),xlim=c(0,1),pch="",
     ylab="Transmission Probability",xlab="Proportion of hosts treated")               
lines(theta36~proportions)
lines(theta37~proportions)
lines(theta38~proportions)
lines(theta40~proportions)
lines(theta41~proportions)
polygon(c(proportions, rev(proportions)),c(theta36Upper,rev(theta36lower)),border=NA, col="aquamarine1")
lines(theta36~proportions)

polygon(c(proportions, rev(proportions)),c(theta37Upper,rev(theta37lower)),border=NA, col="aquamarine1")
lines(theta37~proportions)

polygon(c(proportions, rev(proportions)),c(theta38Upper,rev(theta38lower)),border=NA, col="aquamarine1")
lines(theta38~proportions)

polygon(c(proportions, rev(proportions)),c(theta40Upper,rev(theta40lower)),border=NA, col="aquamarine1")
lines(theta40~proportions)
