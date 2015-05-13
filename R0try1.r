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

###Anderson and May 1985B:
## po <- 1 - (1 + M/k)^-k

po <- numeric(nrow(data))

for (i in 1:206){
    M=data$mean[i]
    k=data$k[i]

    po[i] <- 1 - (1 + M/k)^-k
}

po

prevREAL <- data$prevalence / data$N

plot(po~prevREAL)

log.binom<-function(p.vec){
  
  a<-p.vec[1]
  b<-p.vec[2]
  c<-p.vec[3]
  #pred1a<- ((exp(a + b * prevREAL)) / (1 + exp(a + b * prevREAL)) ) 
  pred1a<- (a * exp (b * exp(c * prevREAL)))
  prev1<-po
  
  loglik1a<- prev1* log((pred1a)+0.00001)+(1-prev1)*log(1-((pred1a)-0.00001))
  -sum(loglik1a,  na.rm=T)
}
#n.param<-2
#logmod<-optim(c(0,0),log.binom,method="L-BFGS-B",lower=c(-10,-10),upper=c(10,100))
#logmod
gommod<-optim(c(0.75,-2,-0.5),log.binom,method="L-BFGS-B",lower=c(0.5,-5,-10),upper=c(0.99,-1,-0.0001))
gommod
nc<-seq(0,1,0.001)
#pred2<-((exp(logmod$par[1] + logmod$par[2] * nc)) / (1 + exp(logmod$par[1] + logmod$par[2] * nc)) ) 
pred2<-(gommod$par[1] * exp (gommod$par[2] * exp(gommod$par[3] * nc)))
lines(nc,pred2,lwd=2,lty=2,col="red")


############################################################################
## From orginal data...and treating most infected through to the least infected

data<-read.csv("C:\\Users\\Ellie\\Documents\\2080\\Exploring R0.csv",header=TRUE)
data[1:10,];summary(data);dim(data)
##Trying for study 1 (Label 36)
pointsdat<-matrix(nrow=max(data$N),ncol=206)
teststore<-matrix(nrow=max(data$N),ncol=206)

proportions<-seq(0,1,length=data[1,4])
plot(proportions~proportions,pch="",ylim=c(0,1),xlim=c(0,1),xaxt="n",par(las=2),
     ylab="Relative reduction in R0 treating most to least infected",xlab="Proportion of hosts treated")
axis(1,at=seq(-0.01,0.99,0.2),labels=seq(0,1,0.2),par(las=1))


##exclude 48, 89:90, 92:93, 95-96, 99-100, 106, 149, 167-169, 180 
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

#lines(smooth.spline(test~proports),col=transp("aquamarine4"),lwd=2)
#}
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

t01R0<-t05R0<-t10R0<-t15R0<-t20R0<-t25R0<-t30R0<-t50R0<-t80R0<-numeric(206)
for(i in 1:206){
  a<-(teststore[,i][round(pointsdat[,i],2)==0.01]);t01R0[i]<-a[1]
  b<-(teststore[,i][round(pointsdat[,i],2)==0.05]);t05R0[i]<-b[1]
  d<-(teststore[,i][round(pointsdat[,i],2)==0.1]);t10R0[i]<-d[1]
  ee<-(teststore[,i][round(pointsdat[,i],2)==0.15]);t15R0[i]<-ee[1]
  g<-(teststore[,i][round(pointsdat[,i],2)==0.2]);t20R0[i]<-g[1]
  h<-(teststore[,i][round(pointsdat[,i],2)==0.25]);t25R0[i]<-h[1]
  m<-(teststore[,i][round(pointsdat[,i],2)==0.30]);t30R0[i]<-m[1]
  bg<-(teststore[,i][round(pointsdat[,i],2)==0.50]);t50R0[i]<-bg[1]
  bm<-(teststore[,i][round(pointsdat[,i],2)==0.50]);t80R0[i]<-bm[1]
}
t01R0<-ifelse(t01R0<0,0,t01R0)  

test2<-c(0,mean(t01R0,na.rm=TRUE),mean(t05R0,na.rm=TRUE),mean(t10R0,na.rm=TRUE),mean(t15R0,na.rm=TRUE),mean(t20R0,na.rm=TRUE),
         mean(t25R0,na.rm=TRUE),rep(mean(t30R0,na.rm=TRUE),4),rep(mean(t50R0,na.rm=TRUE),6),rep(mean(t80R0,na.rm=TRUE),5))
test2up<-c(0,quantile(t01R0,0.90,na.rm=TRUE),quantile(t05R0,0.90,na.rm=TRUE),quantile(t10R0,0.90,na.rm=TRUE),quantile(t15R0,0.90,na.rm=TRUE),quantile(t20R0,0.90,na.rm=TRUE),
         quantile(t25R0,0.90,na.rm=TRUE),rep(quantile(t30R0,0.90,na.rm=TRUE),4),rep(quantile(t50R0,0.90,na.rm=TRUE),6),rep(quantile(t80R0,0.90,na.rm=TRUE),5))
test2low<-c(0,quantile(t01R0,0.025,na.rm=TRUE),quantile(t05R0,0.025,na.rm=TRUE),quantile(t10R0,0.025,na.rm=TRUE),quantile(t15R0,0.025,na.rm=TRUE),quantile(t20R0,0.025,na.rm=TRUE),
         quantile(t25R0,0.025,na.rm=TRUE),rep(quantile(t30R0,0.025,na.rm=TRUE),4),rep(quantile(t50R0,0.025,na.rm=TRUE),6),rep(quantile(t80R0,0.025,na.rm=TRUE),5))
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

plot(nc,pred2,pty="n",pch="",xaxt="n",yaxt="n",ylab="",xlab="")

polygon(c(nc, rev(nc)),c(pred2u,rev(pred2l)),border=NA, col=transp("darkseagreen1",alpha=0.3))
lines(nc,pred2,lwd=3,lty=1,col=transp("aquamarine4"))

par(new=TRUE)
par(bty="n")
boxplot(t01R0,t05R0,t10R0,t15R0,t20R0,t25R0,t30R0,t50R0,t80R0,na.rm=TRUE,
        ylim=c(0,1),col=transp("aquamarine4",alpha=0.2),cex.lab=1.1,at=seq(-0.5,7.5,1),
        ylab=expression(paste("Relative reduction in   ",   R[0])),
        xlab="Proportion of hosts treated from most to least infected")
axis(1,at=seq(-0.5,7.5,1),labels=c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.50,0.80),cex.lab=1.1)

##############################################################################################
##################################################################################################
##
## Exploring the basic reproductive rate R0 when beta, C and D are allowed to 
## follow different distributions
##
#####################################################################
sampnm <- function(n) qnorm(runif(n,min=pnorm(0),max=pnorm(1)))

##what is beta has a normal distribution
beta_n <- sampnm(100)
hist(beta_n)

##what id beta has a negative binomial distribution
k=seq(0.1,5,0.1)
beta_nb<-matrix(nrow=100,ncol=length(k))
beta<-numeric(100)
for(i in 1:length(k)){
  for(j in 1:length(beta)){
beta[j] <- rnegbin(100, 43, k[i])
  beta_nb[,i] <- beta/max(beta)
}
}
BETAmn<-numeric(length(ncol(beta_nb)))
for (i in 1:ncol(beta_nb)) {
     BETAmn[i]<-mean(beta_nb[,i])
}
C=10
D=20

R0_exp<-beta_n * C * D
hist(R0_exp)
R0_exp2<-matrix(nrow=100,ncol=length(k))
for (i in 1:length(ncol(beta_nb))){
R0_exp2[,i]<-beta_nb[,i] * C * D
}

