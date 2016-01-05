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
library(adegenet)
library(binom)
library(data.table)
library(igraph)
library(fields)
library(RColorBrewer)
library(data.table)
library(gridBase)
library(rootSolve)
library(deSolve)


###
### Data
###
###

dat_hosts <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Traits_Host and Parasite.csv",header=TRUE)
dat_paras <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Traits_Host and Parasite2.csv",header=TRUE)
dat_ecol  <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Traits_Ecological.csv",header=TRUE)
dat_2080  <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Exploring R0.csv",header=TRUE)
dat_combo_temp <- merge(dat_hosts,dat_ecol,by="Label")
dat_combo_temp2 <- merge(dat_combo_temp,dat_2080,by="Label") 
dat_combo <- merge(dat_combo_temp2,dat_paras,by="Label")
names(dat_combo_temp2);dim(dat_combo_temp2)

dat <- subset(dat_combo_temp2,dat_combo_temp2$Parasite.Taxa=="Nematoda")
dat2 <- subset(dat_combo,dat_combo$Parasite.Taxa.x=="Nematoda")
dat <- dat[!is.na(dat$mean_host_length),]
dim(dat)
  
#par(mfrow=c(4,4))


exp_dat <- list(N=77,
                x=dat$mean_host_mass,
                y=dat$T20)

exp_datmin <- list(N=77,
                x=dat$min_host_mass,
                y=dat$T20)

exp_datmax <- list(N=77,
                x=dat$max_host_mass,
                y=dat$T20)

##Trying various functional fits for the data
##C:\\Users\\Ellie\\Documents\\RStudioProjects\\TBI_testing_model\\mouse_to_mouse\\Exploring the model\\logisticfunction for prob trans.stan"
##"C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics1.stan"
##"C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics2.stan"
##"C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics3.stan" ##rubbish
##"C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistic4.stan" ##rubbish
##"C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistic5.stan"

test1 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics1.stan", data=exp_dat,
              iter=1000, chains=4)
print(test1)
params1 = extract(test1);names(params1)
#rstan::traceplot(test1, inc_warmup = FALSE)

test1min <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics1.stan", data=exp_datmin,
              iter=1000, chains=4)
print(test1min)
params1min = extract(test1min);names(params1min)
#rstan::traceplot(test1min, inc_warmup = FALSE)

test1max <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics1.stan", data=exp_datmax,
              iter=1000, chains=4)
print(test1max)
params1max = extract(test1max);names(params1max)
#rstan::traceplot(test1max, inc_warmup = FALSE)

## For mice to oocysts in mosquitos
nc<-seq(0,max(dat$mean_host_mass),1)
pred1 <- (mean(params1$alpha) * nc) / sqrt(1 + nc^ (1/mean(params1$beta)))
pred1min <- (mean(params1min$alpha) * nc) / sqrt(1 + nc^ (1/mean(params1min$beta)))
pred1max <- (mean(params1max$alpha) * nc) / sqrt(1 + nc^ (1/mean(params1max$beta)))

##> logistics1
pred2 <- (mean(params1$alpha) * nc^mean(params1$sigma)) / 
  (mean(params1$delta) + mean(params1$beta) * nc^mean(params1$sigma))
pred2min <- (mean(params1min$alpha) * nc^mean(params1min$sigma)) / 
  (mean(params1min$delta) + mean(params1min$beta) * nc^mean(params1min$sigma))
pred2max <- (mean(params1max$alpha) * nc^mean(params1max$sigma)) / 
  (mean(params1max$delta) + mean(params1max$beta) * nc^mean(params1max$sigma))


##> logistics2
pred3 <- (mean(params1$alpha)/mean(params1$beta)) * exp(-exp(mean(params1$delta) - mean(params1$beta)* nc))
pred3min <- (mean(params1min$alpha)/mean(params1min$beta)) * exp(-exp(mean(params1min$delta) - mean(params1min$beta)* nc))
pred3max <- (mean(params1max$alpha)/mean(params1max$beta)) * exp(-exp(mean(params1max$delta) - mean(params1max$beta)* nc))


##> logistics6
pred6 <- 1 / (1 + exp(-(-mean(params1$alpha) + mean(params1$beta) * nc)))
pred6min <- 1 / (1 + exp(-(-mean(params1min$alpha) + mean(params1min$beta) * nc)))
pred6max <- 1 / (1 + exp(-(-mean(params1max$alpha) + mean(params1max$beta) * nc)))


par(mfrow=c(4,5))
plot(dat$T20 ~ 
       dat$mean_host_mass,cex.lab=1.4,
     ylab="T20",ylim=c(0,1),frame=FALSE,
     xlab="Mean host mass (g)")
points(dat$T20 ~ dat$min_host_mass,pch=20,cex=0.5)
points(dat$T20 ~ dat$max_host_mass,pch=20,cex=0.5)
for(i in 1:nrow(dat)) {
  segments(x0=dat$min_host_mass[i],
           x1=dat$max_host_mas[i],
           y0=dat$T20[i],
           y1=dat$T20[i],lty=2,col="grey")
}
#lines(pred1 ~ nc)
#lines(pred1min ~ nc)
#lines(pred1max ~ nc)

lines(pred2 ~ nc,col="red")
#lines(pred2min ~ nc,col="red")
#lines(pred2max ~ nc,col="red")

plot(dat$T20 ~ 
       dat$mean_host_mass,cex.lab=1.4,
     ylab="T20",ylim=c(0,1),frame=FALSE,
     xlab="Mean host mass (g)",XLIM=c(0,8000))
points(dat$T20 ~ dat$min_host_mass,pch=20,cex=0.5)
points(dat$T20 ~ dat$max_host_mass,pch=20,cex=0.5)
for(i in 1:nrow(dat)) {
  segments(x0=dat$min_host_mass[i],
           x1=dat$max_host_mas[i],
           y0=dat$T20[i],
           y1=dat$T20[i],lty=2,col="grey")
}
plot(dat$T20 ~ 
       dat$mean_host_length,,main="NEMATODA",
     ylab="T20",cex.lab=1.4,ylim=c(0,1),frame=FALSE,
     xlab="Mean host length (cm)")
points(dat$T20 ~ dat$min_host_length,pch=20,cex=0.5)
points(dat$T20 ~ dat$max_host_length,pch=20,cex=0.5)
for(i in 1:nrow(dat)) {
  segments(x0=dat$min_host_length[i],
           x1=dat$max_host_length[i],
           y0=dat$T20[i],
           y1=dat$T20[i],lty=2,col="grey")
}

plot(dat$T20 ~ 
       dat$dLat,ylim=c(0,1),frame=FALSE,
     ylab="T20",cex.lab=1.4,
     xlab="Latitude")
plot(dat$T20 ~ 
       dat$dLong,ylim=c(0,1),frame=FALSE,
     ylab="T20",cex.lab=1.4,
     xlab="Longitude")
labnames <- c(colnames(dat[,45:57]))
for(i in 45:57) {
  plot(dat$T20 ~ 
         dat[,i],ylim=c(0,1),frame=FALSE,
       ylab="T20",cex.lab=1.4,
       xlab=substitute(paste(a), list(a = labnames[i-44])))
}

  plot(dat2$T20 ~ 
         dat2$mean_mm_lengthF,
       ylab="T20",cex.lab=1.4,ylim=c(0,1),frame=FALSE,
       xlab="Mean length female parasite",xlim=c(0,100))
  points(dat2$T20 ~ dat2$min_length,pch=20,cex=0.5)
  points(dat2$T20 ~ dat2$max_length,pch=20,cex=0.5)
for(i in 1:nrow(dat2)) {
  segments(x0=dat2$min_length[i],
           x1=dat2$max_length[i],
           y0=dat2$T20[i],
           y1=dat2$T20[i],lty=2,col="grey")
}

plot(dat2$T20 ~ 
       dat2$mean_mm_lengthF,
     ylab="T20",cex.lab=1.4,ylim=c(0,1),frame=FALSE,
     xlab="Mean length female parasite")
points(dat2$T20 ~ dat2$min_length,pch=20,cex=0.5)
points(dat2$T20 ~ dat2$max_length,pch=20,cex=0.5)
for(i in 1:nrow(dat2)) {
  segments(x0=dat2$min_length[i],
           x1=dat2$max_length[i],
           y0=dat2$T20[i],
           y1=dat2$T20[i],lty=2,col="grey")
}