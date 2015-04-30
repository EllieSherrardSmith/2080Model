###############################################
##
###
#### And Now for R0
###
##
################################################
###using data from CopyOfExplorationR02.r

###Lloyd-Smith et al 2005:
## po <- (1 + R0/k) ^ -k
R0<-seq(0.1,5,0.1)
k<-seq(0.01,10,0.2)
po<-expand.grid(c(1:50))
for (i in 1:length(R0)){
  for (j in 1:length(k)){
    po[j,i]<-(1 + R0[i]/k[j]) ^ -k[j]
  }
}
po[,1:6]

plot(po[,1]~R0,ylim=c(0,1))
for(i in 1:ncol(po)){ lines(po[,i]~R0)}

plot(po[,1]~k,ylim=c(0,1))
for(i in 1:ncol(po)){ lines(po[,i]~k)}

plot(k~R0)
for(i in 1:ncol(po)){ lines(k[,i]~R0)}



##R0<- (k * ( po ^ (-1/k))) - k

par(mfrow=c(1,1))
#if po = 0.2 vary k
po<-c(0.01,seq(0.05,0.95,0.05),0.99)
k<-seq(0.01,10,0.02)

R0<-expand.grid(c(1:length(k)))

for (i in 1:length(k)){
  for(j in 1:length(po)){
  R0[i,j]<- (k[i] * ( po[j] ^ (-1/k[i]))) - k[i]
  }  
}

plot(R0[,1]~k,ylim=c(0,20),pch="")
for(i in 1:ncol(R0)){
lines(R0[,i]~k,col="grey60")
}
text(8,1.5,"0.05 prev")

############################################################################
## From orginal data...and treating most infected through to the least infected

data<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\Exploring R0.csv",header=TRUE)
data[1:10,];summary(data);dim(data)
##Trying for study 1 (Label 36)

proportions<-seq(0,1,length=data[1,4])
plot(proportions~proportions,pch="",ylim=c(0,1),xlim=c(0,1),xaxt="n",par(las=2),
     ylab="Relative reduction in R0 treating most to least infected",xlab="Proportion of hosts treated")
axis(1,at=seq(-0.01,0.99,0.2),labels=seq(0,1,0.2),par(las=1))
##exclude 9-10, 46-48, 56-58, 68, 89, 90, 92-93, 95-96, 99-100, 106, 149, 167-169, 180 
for (j in 181:206){

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
}

for (i in 2:length(R0)){
  test[1]<-0
test[i] <-  (R0[1]-R0[i])/R0[1]
}

proports<-seq(0,1,length=data[j,4])

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


segments(x0=0.05,y0=0,x1=1,y1=0,col="white",lwd=10)
abline(v=0.19,lty=2,col="grey")
