library(fields)
library(RColorBrewer)
library(nlme)
library(rstan)
library(MASS)
library(boot)
library(coda)
library(R2OpenBUGS)
library(ggplot2)
library("Rlab")
library(contrast)
#install.packages('devtools')
#library(devtools)
#source_url("https://github.com/stan-dev/shinystan/raw/develop/install_shinystan.R")
#install_shinystan()
library(shinyStan)
library(rmngb)


data<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\DataFinalApril2015.csv",header=TRUE)
names(data)

plot(t20 ~ log(k), data=data)
##################
###
####
#####
###### HOST AND PARASITE TAXA
#####
####
###
#########################
dat2 <- data.frame(data$Parasite.taxonomy,data$Host.taxonomy,data$t20,data$"mean..mu.",data$k,data$N)
colnames(dat2) <- c("Parasite_type","Host_type","T20","Mean","k","n")

library(plyr)
dat3 <- ddply(dat2, .(Host_type, Parasite_type), summarise, 
        T20 = mean(T20),
        k = mean(k),
        Mean = mean(Mean))
dat3$CodeHosts <- c(5,5,5,2,2,2,4,4,4,4,6,1,1,1,3)
dat3$CodePara <- c(1,3,4,2,3,4,1,2,3,4,3,2,3,4,3)

DAT4 <- dat3[with(dat3, order(CodeHosts, CodePara)), ]
meansT20 <- c(NA,DAT4$T20[1:3],NA,DAT4$T20[4:6],NA,NA,DAT4$T20[7],NA,DAT4$T20[8:12],NA,DAT4$T20[13:14],NA,NA,DAT4$T20[15],NA)
meansk <- c(NA,DAT4$k[1:3],NA,DAT4$k[4:6],NA,NA,DAT4$k[7],NA,DAT4$k[8:12],NA,DAT4$k[13:14],NA,NA,DAT4$k[15],NA)
meansmean <- c(NA,DAT4$Mean[1:3],NA,DAT4$Mean[4:6],NA,NA,DAT4$Mean[7],NA,DAT4$Mean[8:12],NA,DAT4$Mean[13:14],NA,NA,DAT4$Mean[15],NA)

tapply(data$Parasite.taxonomy,data$Host.taxonomy,summary)
Ngp <- c(0,21,61,16,0,11,23,19,
         0,0,14,0,7,1,4,21,
         3,0,19,9,0,0,4,0)

FinalT20 <- matrix(meansT20,nrow=4,ncol=6)
Finalk <- matrix(meansk,nrow=4,ncol=6)
FinalMean <- matrix(meansmean,nrow=4,ncol=6)
colnames(FinalT20) <- colnames(Finalk) <- colnames(FinalMean) <- hosts <- c("Mammals","Birds","Reptiles","Fish","Amphibians","Inverts")
rownames(FinalT20) <- rownames(Finalk) <- rownames(FinalMean) <- parasites <- c("Acanthocephala","Arthropoda","Nematoda","Platyhelminthes")

par(mfrow=c(1,3))
par(mar=c(12,10,5,5))
image(matrix(FinalT20,nrow=4),axes=F,main="Mean T20 per group",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(hosts), side=2, line=0.3, at=seq(0,1,0.2), las=1, cex=1.2)
mtext(text=c(parasites),side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
image.plot(matrix(FinalT20,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))
text(x=seq(0,1,0.33),y=rep(0,4),c(Ngp[1:4]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(0.2,4),c(Ngp[5:8]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(0.4,4),c(Ngp[9:12]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(0.6,4),c(Ngp[13:16]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(0.8,4),c(Ngp[17:20]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(1.0,4),c(Ngp[21:24]),cex=1.2)

image(matrix(Finalk,nrow=4),axes=F,main="Mean k per group",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(hosts), side=2, line=0.3, at=seq(0,1,0.2), las=1, cex=1.2)
mtext(text=c(parasites), side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
image.plot(matrix(Finalk,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))

image(matrix(FinalMean,nrow=4),axes=F,main="Mean Mean per group",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(hosts), side=2, line=0.3, at=seq(0,1,0.2), las=1, cex=1.2)
mtext(text=c(parasites), side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
image.plot(matrix(FinalMean,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))


##################
###
####
#####
###### Parasite TAXA AND Host MANAGEMENT 
#####
####
###
#########################
dat2 <- data.frame(data$Host.management,data$Host.habitat,data$Parasite.taxonomy,data$t20,data$"mean..mu.",data$k,data$N)
colnames(dat2) <- c("Host_management","Host_habitat","Parasite_type","T20","Mean","k","n")

dat2 <- dat2[with(dat2, order(dat2$T20)), ]

par(mar=c(2,4,5,5))
par(mfrow=c(2,2))
manage <- c("(91) Wild (31)","(6) Seminatural (5)","(21) Captive (5)","(7) Experimental (0)")
dat2nem <- subset(dat2,Parasite_type == "Nematoda")
dat2nemMAT <- c(dat2nem$T20[dat2nem$Host_management == "wild"],
              c(rep(NA,91-sum(ifelse(dat2nem$Host_management == "seminat",1,0))),
                dat2nem$T20[dat2nem$Host_management == "seminat"]),
              c(rep(NA,91-sum(ifelse(dat2nem$Host_management == "captive",1,0))),
                dat2nem$T20[dat2nem$Host_management == "captive"]),
              c(rep(NA,91-sum(ifelse(dat2nem$Host_management == "exp",1,0))),
                dat2nem$T20[dat2nem$Host_management == "exp"]))
dat2nemMAT <- matrix(dat2nemMAT,nrow=91,ncol=4)
image(dat2nemMAT,axes=F,main="Nematodes",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c(manage),side=2, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)

dat2pla <- subset(dat2,Parasite_type == "Platyhelminthes")
sum(ifelse(dat2pla$Host_management == "exp",1,0))
par(mar=c(2,8,5,5))
dat2plaMAT <- c(dat2pla$T20[dat2pla$Host_management == "wild"],
                c(rep(NA,55-sum(ifelse(dat2pla$Host_management == "seminat",1,0))),
                  dat2pla$T20[dat2pla$Host_management == "seminat"]),
                c(rep(NA,55-sum(ifelse(dat2pla$Host_management == "captive",1,0))),
                  dat2pla$T20[dat2pla$Host_management == "captive"]),
                c(rep(NA,55-sum(ifelse(dat2pla$Host_management == "exp",1,0))),
                  dat2pla$T20[dat2pla$Host_management == "exp"]))
dat2plaMAT <- matrix(dat2plaMAT,nrow=55,ncol=4)
image(dat2plaMAT,axes=F,main="Platyhelminthes",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(manage),side=2, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
image.plot(dat2plaMAT, legend.only=T,col=brewer.pal(9, "YlGnBu"))

par(mar=c(2,4,5,5))
hab <- c("(3) Marine (8)","(0) Freshwater/Marine (3)","(3) Estuarine (3)",
         "(1) Freshwater (11)","(16) Semi-aquatic (9)","(102) Terrestrial (31)")
dat2nem <- subset(dat2,Parasite_type == "Nematoda")
dat2nemMAT <- c(c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Marine",1,0))),
                  dat2nem$T20[dat2nem$Host_habitat == "Marine"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Freshwater/Marine",1,0))),
                  dat2nem$T20[dat2nem$Host_habitat == "Freshwater/Marine"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Estuarine",1,0))),
                  dat2nem$T20[dat2nem$Host_habitat == "Estuarine"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Freshwater",1,0))),
                  dat2nem$T20[dat2nem$Host_habitat == "Freshwater"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Semi-aquatic",1,0))),
                  dat2nem$T20[dat2nem$Host_habitat == "Semi-aquatic"]),
                dat2nem$T20[dat2nem$Host_habitat == "Terrestrial"])
dat2nemMAT <- matrix(dat2nemMAT,nrow=102,ncol=6)
image(dat2nemMAT,axes=F,main="",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c(hab),side=2, line=0.3, at=seq(0,1,0.2), las=2, cex=1.2)

par(mar=c(2,8,5,5))
dat2pla <- subset(dat2,Parasite_type == "Platyhelminthes")
dat2plaMAT <- c(c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Marine",1,0))),
                  dat2pla$T20[dat2pla$Host_habitat == "Marine"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Freshwater/Marine",1,0))),
                  dat2pla$T20[dat2pla$Host_habitat == "Freshwater/Marine"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Estuarine",1,0))),
                  dat2pla$T20[dat2pla$Host_habitat == "Estuarine"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Freshwater",1,0))),
                  dat2pla$T20[dat2pla$Host_habitat == "Freshwater"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Semi-aquatic",1,0))),
                  dat2pla$T20[dat2pla$Host_habitat == "Semi-aquatic"]),
                dat2pla$T20[dat2pla$Host_habitat == "Terrestrial"])
dat2plaMAT <- matrix(dat2plaMAT,nrow=31,ncol=6)
image(dat2plaMAT,axes=F,main="",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(hab),side=2, line=0.3, at=seq(0,1,0.2), las=2, cex=1.2)
image.plot(dat2plaMAT, legend.only=T,col=brewer.pal(9, "YlGnBu"))

dat2 <- data.frame(data$Host.management,data$Host.habitat,data$Parasite.taxonomy,data$t20,data$"mean..mu.",data$k,data$N)
colnames(dat2) <- c("Host_management","Host_habitat","Parasite_type","T20","Mean","k","n")

dat2 <- dat2[with(dat2, order(dat2$k)), ]

par(mar=c(2,4,5,5))
par(mfrow=c(2,2))
manage <- c("(91) Wild (31)","(6) Seminatural (5)","(21) Captive (5)","(7) Experimental (0)")
dat2nem <- subset(dat2,Parasite_type == "Nematoda")
dat2nemMAT <- c(dat2nem$k[dat2nem$Host_management == "wild"],
                c(rep(NA,91-sum(ifelse(dat2nem$Host_management == "seminat",1,0))),
                  dat2nem$k[dat2nem$Host_management == "seminat"]),
                c(rep(NA,91-sum(ifelse(dat2nem$Host_management == "captive",1,0))),
                  dat2nem$k[dat2nem$Host_management == "captive"]),
                c(rep(NA,91-sum(ifelse(dat2nem$Host_management == "exp",1,0))),
                  dat2nem$k[dat2nem$Host_management == "exp"]))
dat2nemMAT <- matrix(dat2nemMAT,nrow=91,ncol=4)
image(log(dat2nemMAT),axes=F,main="Nematodes",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c(manage),side=2, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)

dat2pla <- subset(dat2,Parasite_type == "Platyhelminthes")
sum(ifelse(dat2pla$Host_management == "exp",1,0))
par(mar=c(2,8,5,5))
dat2plaMAT <- c(dat2pla$k[dat2pla$Host_management == "wild"],
                c(rep(NA,55-sum(ifelse(dat2pla$Host_management == "seminat",1,0))),
                  dat2pla$k[dat2pla$Host_management == "seminat"]),
                c(rep(NA,55-sum(ifelse(dat2pla$Host_management == "captive",1,0))),
                  dat2pla$k[dat2pla$Host_management == "captive"]),
                c(rep(NA,55-sum(ifelse(dat2pla$Host_management == "exp",1,0))),
                  dat2pla$k[dat2pla$Host_management == "exp"]))
dat2plaMAT <- matrix(dat2plaMAT,nrow=55,ncol=4)
image(log(dat2plaMAT),axes=F,main="Platyhelminthes",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(manage),side=2, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
image.plot(log(dat2plaMAT), legend.only=T,col=brewer.pal(9, "YlGnBu"))

par(mar=c(2,4,5,5))
hab <- c("(3) Marine (8)","(0) Freshwater/Marine (3)","(3) Estuarine (3)",
         "(1) Freshwater (11)","(16) Semi-aquatic (9)","(102) Terrestrial (31)")
dat2nem <- subset(dat2,Parasite_type == "Nematoda")
dat2nemMAT <- c(c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Marine",1,0))),
                  dat2nem$k[dat2nem$Host_habitat == "Marine"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Freshwater/Marine",1,0))),
                  dat2nem$k[dat2nem$Host_habitat == "Freshwater/Marine"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Estuarine",1,0))),
                  dat2nem$k[dat2nem$Host_habitat == "Estuarine"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Freshwater",1,0))),
                  dat2nem$k[dat2nem$Host_habitat == "Freshwater"]),
                c(rep(NA,102-sum(ifelse(dat2nem$Host_habitat == "Semi-aquatic",1,0))),
                  dat2nem$k[dat2nem$Host_habitat == "Semi-aquatic"]),
                dat2nem$k[dat2nem$Host_habitat == "Terrestrial"])
dat2nemMAT <- matrix(dat2nemMAT,nrow=102,ncol=6)
image(log(dat2nemMAT),axes=F,main="",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c(hab),side=2, line=0.3, at=seq(0,1,0.2), las=2, cex=1.2)

par(mar=c(2,8,5,5))
dat2pla <- subset(dat2,Parasite_type == "Platyhelminthes")
dat2plaMAT <- c(c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Marine",1,0))),
                  dat2pla$k[dat2pla$Host_habitat == "Marine"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Freshwater/Marine",1,0))),
                  dat2pla$k[dat2pla$Host_habitat == "Freshwater/Marine"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Estuarine",1,0))),
                  dat2pla$k[dat2pla$Host_habitat == "Estuarine"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Freshwater",1,0))),
                  dat2pla$k[dat2pla$Host_habitat == "Freshwater"]),
                c(rep(NA,31-sum(ifelse(dat2pla$Host_habitat == "Semi-aquatic",1,0))),
                  dat2pla$k[dat2pla$Host_habitat == "Semi-aquatic"]),
                dat2pla$k[dat2pla$Host_habitat == "Terrestrial"])
dat2plaMAT <- matrix(dat2plaMAT,nrow=31,ncol=6)
image(log(dat2plaMAT),axes=F,main="",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(hab),side=2, line=0.3, at=seq(0,1,0.2), las=2, cex=1.2)
image.plot(log(dat2plaMAT), legend.only=T,col=brewer.pal(9, "YlGnBu"))


library(plyr)
dat3 <- ddply(dat2, .(Parasite_type, Host_management), summarise, 
              T20 = mean(T20),
              k = mean(k),
              Mean = mean(Mean))
dat3$CodeHosts <- c(1,2,2,2,2,3,3,3,3,4,4,4)
dat3$Codemanage <- c(1,3,4,2,1,3,4,2,1,3,2,1)

DAT4 <- dat3[with(dat3, order(CodeHosts, Codemanage)), ]
meansT20 <- c(DAT4$T20[1],NA,NA,NA,DAT4$T20[2:12],NA)
meansk <- c(DAT4$k[1],NA,NA,NA,DAT4$k[2:12],NA)
meansmean <- c(DAT4$Mean[1],NA,NA,NA,DAT4$Mean[2:12],NA)

tapply(data$Host.management,data$Parasite.taxonomy,summary)
Ngp <- c(10,0,0,0,29,1,2,1,91,21,6,7,55,5,5,0)

FinalT20 <- matrix(meansT20,nrow=4,ncol=4)
Finalk <- matrix(meansk,nrow=4,ncol=4)
FinalMean <- matrix(meansmean,nrow=4,ncol=4)
colnames(FinalT20) <- colnames(Finalk) <- colnames(FinalMean) <- hosts <- c("Acanthocephala","Arthropoda","Nematoda","Platyhelminthes")
rownames(FinalT20) <- rownames(Finalk) <- rownames(FinalMean) <- manage <- c("Wild","Semi-natural","Captive","Experimental")

par(mfrow=c(1,3))
par(mar=c(12,10,5,5))
image(matrix(FinalT20,nrow=4),axes=F,main="Mean T20 per group",col=brewer.pal(9, "YlGnBu"))
mtext(text=c(hosts), side=2, line=0.3, at=seq(0,1,0.33), las=1, cex=1.2)
mtext(text=c(manage),side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
image.plot(matrix(FinalT20,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))
text(x=seq(0,1,0.33),y=rep(0,4),c(Ngp[1:4]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(0.33,4),c(Ngp[5:8]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(0.66,4),c(Ngp[9:12]),cex=1.2)
text(x=seq(0,1,0.33),y=rep(1.0,4),c(Ngp[13:16]),cex=1.2)

#image(matrix(Finalk,nrow=4),axes=F,main="Mean k per group",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c(hosts), side=2, line=0.3, at=seq(0,1,0.33), las=1, cex=1.2)
#mtext(text=c(manage), side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
#image.plot(matrix(Finalk,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))

#image(matrix(FinalMean,nrow=4),axes=F,main="Mean Mean per group",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c(hosts), side=2, line=0.3, at=seq(0,1,0.33), las=1, cex=1.2)
#mtext(text=c(manage), side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)
#image.plot(matrix(FinalMean,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))


dat3 <- ddply(data, .(Lab.or.field, Host.habitat), summarise, 
              T20 = mean(t20),
              k = mean(k),
              GMean = mean(Gini.Co.efficient))
meansT20 <- c(dat3$T20[1:6],NA,dat3$T20[7],NA,NA,NA,dat3$T20[8])
meansk <- c(dat3$k[1:6],NA,dat3$k[7],NA,NA,NA,dat3$k[8])
meansGmean <- c(dat3$GMean[1:6],NA,dat3$GMean[7],NA,NA,NA,dat3$GMean[8])

image(matrix(meansT20,nrow=6),axes=F,main="T20",col=brewer.pal(9, "YlGnBu"))
mtext(text=c("Field","Lab"), side=2, line=0.3, at=c(0,1), las=2, cex=1.2)
mtext(text=c("Estuarine","Freshwater","Freshwater/Marine",
             "Marine","Semi-aquatic","Terrestrial"), side=1, line=0.3, at=c(0,0.2,0.4,0.6,0.8,1), las=2, cex=1.2)
tapply(data$Lab.or.field,data$Host.habitat,summary)
text(x = seq(0,1,0.2),y=rep(0,6),c(8,15,4,12,25,156))
text(x = seq(0,1,0.2),y=rep(1,6),c(0,1,0,0,0,12))
image.plot(matrix(meansT20,nrow=6), legend.only=T,col=brewer.pal(9, "YlGnBu"))

#image(matrix(meansk,nrow=6),axes=F,main="k",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c("Field","Lab"), side=2, line=0.3, at=c(0,1), las=2, cex=1.2)
#mtext(text=c("Estuarine","Freshwater","Freshwater/Marine",
#             "Marine","Semi-aquatic","Terrestrial"), side=1, line=0.3, at=c(0,0.2,0.4,0.6,0.8,1), las=2, cex=1.2)

#image(matrix(meansGmean,nrow=6),axes=F,main="Gini",col=brewer.pal(9, "YlGnBu"))
#mtext(text=c("Field","Lab"), side=2, line=0.3, at=c(0,1), las=2, cex=1.2)
#mtext(text=c("Estuarine","Freshwater","Freshwater/Marine",
#             "Marine","Semi-aquatic","Terrestrial"), side=1, line=0.3, at=c(0,0.2,0.4,0.6,0.8,1), las=2, cex=1.2)

dat3 <- ddply(data, .(Social.structure, Transmission.route), summarise, 
              T20 = mean(t20),
              k = mean(k),
              GMean = mean(Gini.Co.efficient))
dat3$Codesoc <- c(4,4,3,3,3,3,1,1,1,2,2,2,5,5,5,5)
dat3$Coderoute <- c(2,3,1,2,3,4,1,2,3,1,2,3,1,2,3,4)

DAT4 <- dat3[with(dat3, order(Codesoc, Coderoute)), ]

meansT20 <- c(DAT4$T20[1:3],NA,DAT4$T20[4:6],NA,DAT4$T20[7:10],NA,DAT4$T20[11:12],NA,DAT4$T20[13:16])
meansk <- c(DAT4$k[1:3],NA,DAT4$k[4:6],NA,DAT4$k[7:10],NA,DAT4$k[11:12],NA,DAT4$k[13:16])
meansGmean <- c(DAT4$Gmean[1:3],NA,DAT4$Gmean[4:6],NA,DAT4$Gmean[7:10],NA,DAT4$Gmean[11:12],NA,DAT4$T20[13:16])

image(matrix(meansT20,nrow=4),axes=F,main="T20",col=brewer.pal(9, "YlGnBu"))
mtext(text=c("Solitary","Sol (n)/ Social (b)","Social","Colonial","Other"), side=2, line=0.3, at=seq(0,1,0.25), las=2, cex=1.2)
mtext(text=c("Active","Eats parasite","Eats prey",
             "Other"), side=1, line=0.3, at=seq(0,1,0.33), las=2, cex=1.2)

tapply(data$Social.structure,data$Transmission.route,summary)
text(x = seq(0,1,0.33),y=rep(0,4),c(13,17,38,0))
text(x = seq(0,1,0.33),y=rep(0.25,4),c(3,10,18,0))
text(x = seq(0,1,0.33),y=rep(0.50,4),c(22,43,37,17))
text(x = seq(0,1,0.33),y=rep(0.75,4),c(0,4,2,0))
text(x = seq(0,1,0.33),y=rep(1,4),c(3,3,2,1))
image.plot(matrix(meansT20,nrow=4), legend.only=T,col=brewer.pal(9, "YlGnBu"))
