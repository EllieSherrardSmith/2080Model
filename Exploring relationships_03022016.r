######################################################
##  An analysis of traits data to determine what 
##  impacts on the aggregation of parasites in 
##  host populations. Using T20

######################################################
##  Essential packages

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

######################################################
###                                                ###
###  ####      ##    ##########   ##               ### 
###  ## ##    ####       ##      ####              ### 
###  ##  ##  ##  ##      ##     ##  ##             ###
###  ##  ## ########     ##    ########            ###
###  ## ##  ##    ##     ##    ##    ##            ###
###  ####  ##      ##    ##   ##      ##           ###
###                                                ###
######################################################

dat_hosts <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Traits_Host and Parasite.csv",header=TRUE)
dat_paras <-read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Traits_parasite measurements.csv",header=TRUE)
dat_ecol  <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Traits_Ecological.csv",header=TRUE)
dat_2080  <- read.csv("C:\\Users\\Ellie\\Documents\\Traits\\data_csvs\\Exploring R0.csv",header=TRUE)

dat_combo_temp <- merge(dat_hosts,dat_ecol,by="Label")
dat_combo_temp2 <- merge(dat_combo_temp,dat_2080,by="Label") 
dat_combo <- merge(dat_combo_temp2,dat_paras,by="Label")

dat_nem <- subset(dat_combo,dat_combo$Parasite.Taxa.x=="Nematoda")
dat_pla <- subset(dat_combo,dat_combo$Parasite.Taxa.x=="Platyhelminthes")



######################################################
## The Nematoda

names(dat_nem)

######################################################
## prior exploration
nam   <- c(colnames(dat_nem[,25:43]))
par(mfrow=c(4,5))
for(i in 1:19){
  plot(dat_nem$T20 ~ dat_nem[,24+i],
       xlab=paste(nam[i]))
  abline(lm(dat_nem$T20 ~ dat_nem[,24+i]))
  lm.res <- lm(dat_nem$T20 ~ dat_nem[,24+i])
  print(summary(lm.res)$coefficients[,4] )
}
## bio 5, 10 and 15 are linearly significaiton


nam   <- c(colnames(dat_nem[,44:57]))
par(mfrow=c(3,5))
for(i in 1:14){
  plot(dat_nem$T20 ~ dat_nem[,i+43],
       xlab=paste(nam[i]))
  abline(lm(dat_nem$T20 ~ dat_nem[,i+43]))
  lm.res <- lm(dat_nem$T20 ~ dat_nem[,i+43])
  print(summary(lm.res)$coefficients[,4] )
}


##
###
####
##### The data are skewed such that a ligistic function would fit better than a linear one 

###   First fit functions for y ~ x then work toward a multiple explanatory variables scenario
##


###  Create the data_lists for analysis
##

par(mfrow=c(2,4))

## 1. Parasite traits
  dat_nem_parasitelength <- list(N=58,
                                 x=dat_nem$mean_para_length_F_mm,
                                 y=dat_nem$T20)
## 2. Host traits
  dat_nem_host_mass      <- list(N=58,
                                 x=dat_nem$mean_host_mass,
                                 y=dat_nem$T20)
## 3. Ecological traits
  dat_nem_LONGITUDE      <- list(N=58,
                                 x=(dat_nem$dLong)^2,
                                 y=dat_nem$T20)
  dat_nem_LATITUDE       <- list(N=58,
                                 x=(dat_nem$dLat)^2,
                                 y=dat_nem$T20)
  dat_nem_ai_year        <- list(N=58,
                                 x=dat_nem$ai_yr,
                                 y=dat_nem$T20)  
  dat_nem_bio_5          <- list(N=58,
                                 x=dat_nem$bio_5,
                                 y=dat_nem$T20)
  dat_nem_bio_10         <- list(N=58,
                                 x=dat_nem$bio_10,
                                 y=dat_nem$T20)
  dat_nem_bio_15         <- list(N=58,
                                 x=dat_nem$bio_15,
                                 y=dat_nem$T20)

plotting_funcs <- function(data_list, name_x) {

  
  ## create a variable from which to draw predictive lines
  nc<-seq(0,max(data_list$x),1)
  
  ## fit functions to the data
  ## Function 1: 
  test1 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics simple log function.stan", data=data_list,
                iter=1000, chains=4)
  params1 = extract(test1);names(params1)
  print(test1)
  
  pred1 <- 1 / (1 + exp(mean(params1$alpha) * nc))
  
  ## Function 2:
  test2 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics1.stan", data=data_list,
                iter=1000, chains=4)
  params2 = extract(test2);names(params2)
  print(test2)
  
  pred2 <- (mean(params2$alpha) * nc ^ mean(params2$sigma))/(1 + mean(params2$beta) * nc ^ mean(params2$sigma))
  
  ## Function 3:
  test3 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics1a.stan", data=data_list,
                iter=1000, chains=1)
  params3 = extract(test3);names(params3)
  print(test3)

  pred3 <- (mean(params3$alpha) * nc^mean(params3$sigma))/(mean(params3$delta) + mean(params3$beta) * nc^mean(params3$sigma))


  ## Function 4:
  test4 <- stan(file="C:\\Users\\Ellie\\Documents\\RStudioProjects\\2080Model\\logistics2.stan", data=data_list,
               iter=1000, chains=1)
  params4 = extract(test4);names(params4)
  print(test4)

  pred4 <- (mean(params4$alpha)/mean(params4$beta)) * exp(-exp(mean(params4$delta) - mean(params4$beta) * nc))

  plot(data_list$y ~ 
         data_list$x,,cex.lab=1.4,
       ylab="T20",ylim=c(0,1),frame=FALSE,
       xlab=print(name_x))
  lines(pred1 ~ nc,col="red")
  lines(pred2 ~ nc,col = "red", lty = 2)
  lines(pred3 ~ nc,col = "darkred", lty = 2)
  lines(pred4 ~ nc,col = "darkred", lty = 3)
  
  abline(lm(data_list$y ~ data_list$x),col="blue",lty=2)
  
}

plotting_funcs(dat_nem_parasitelength,name_x="Mean nematode length (mm)")
plotting_funcs(dat_nem_host_mass,name_x="Mean host mass (g) in nematodes")
plotting_funcs(dat_nem_LONGITUDE,name_x="Longitude")
plotting_funcs(dat_nem_LATITUDE,name_x="Latitude")

plotting_funcs(dat_nem_ai_year,name_x="ai_year")
plotting_funcs(dat_nem_bio_5,name_x="bio 5")
plotting_funcs(dat_nem_bio_10,name_x="bio 10")
plotting_funcs(dat_nem_bio_15,name_x="bio 15")
