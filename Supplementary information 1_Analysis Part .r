library(MASS)
library(reshape2)




##################################################################
## A step by step run through of the data analysis for 20:80 #####
##################################################################

#We want to explore whether parasite populations fit the '20:80' distribution
#where 80% of the parasites are contained within 20% of the hosts.
#First, we explore real data on this relationship and then compare these 
#distributions to those that are generated using the mean, N of hosts, and
#a measure of variance.

###############
##Step 1#######
###############

###Real Data; supplied by Darren Shaw September 2013

AclavKennedy84 <- rep(c(seq(0,10,1),20),c(56,29,13,10,9,4,3,2,3,1,2,8))
AgallHodasi69 <- rep(c(seq(0,5,1),13),c(79,11,4,3,2,5,3))
AluciiBrattey88 <- rep(c(seq(0,65,1),119),c(84,57,36,29,29,22,14,14,17,11,2,13,4,4,6,6,4,3,1,2,4,1,2,3,4,4,3,2,3,1,4,1,4,1,1,3,0,0,0,1,2,1,0,0,0,1,1,0,0,0,2,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,1))
BfuscStarJames67 <- rep(seq(0,20,2),c(97,16,4,2,0,0,0,0,1,0,2))
BfuscRedJames67 <- rep(c(seq(0,20,2),40,50),c(114,18,7,10,1,4,1,1,2,0,7,1,1))
ChamulHodasi69 <- rep(seq(0,10,1),c(81,11,4,6,0,1,4,0,0,0,1))
CnodosMarWilliams63 <- rep(seq(0,5,1),c(496,58,18,7,2,4))
CnodosMenWilliams63 <- rep(seq(0,8,1),c(940,283,93,33,18,1,4,6,2))
CsemerNuorteva66 <- rep(seq(0,10,1),c(47,36,35,28,18,10,8,2,3,2,1))
DinterBreyev73 <- rep(seq(0,5,1),c(90,6,3,2,0,1))
HtaranShaw <- rep(c(seq(0,78,1),seq(80,92,1),seq(94,100),102,104,105,106,108,
110,111,112,114,116,117,118,119,120,121,122,124,125,126,128,130,132,134,136,
138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,
176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,
214,218,220,222,224,226,230,232,234,236,238,240,244,246,248,254,260,264,270,
274,280,282,284,290,304,320,322,336,340,360,384,420,432,600,640),
c(30,23,27,29,40,32,37,32,29,38,32,32,34,45,35,32,39,25,40,34,39,36,37,18,19,
27,34,22,38,20,24,22,21,19,25,19,35,21,35,12,33,19,38,21,29,15,26,10,28,11,33,
8,22,8,28,4,25,13,26,5,29,8,30,12,26,7,26,7,27,7,19,4,19,4,17,4,21,4,21,21,2,
18,7,25,2,18,3,19,2,17,2,19,6,1,17,2,17,3,15,23,16,1,8,11,12,1,10,6,8,1,9,1,9,
1,12,9,2,5,7,13,12,9,11,6,8,9,8,6,4,9,6,3,4,5,8,8,1,5,6,6,4,2,3,1,5,1,5,4,2,5,
3,4,1,5,2,2,3,2,3,2,4,5,1,5,2,2,2,4,1,2,1,3,7,4,2,1,1,2,1,1,2,2,1,1,1,1,5,1,1,1,3,1,1,1,1,1))
IricinMilne43 <- rep(c(seq(0,90,1),92,100,116,161,162,168,180),c(59,44,33,33,29,23,17,13,11,14,8,13,5,8,6,6,5,7,5,4,5,11,5,5,5,10,5,4,3,3,5,7,4,4,5,2,4,3,3,4,1,1,2,1,3,2,0,4,3,2,2,1,1,4,1,0,0,0,1,0,2,0,2,1,1,1,2,2,0,0,1,0,0,0,1,0,0,0,2,0,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,1))
PensicJames67 <- rep(seq(0,7,1),c(89,19,8,2,2,1,0,1))
LlusciEvans83 <- rep(seq(0,5,1),c(191,64,18,14,4,2))
LpectBoxshall74 <- rep(seq(0,26,1),c(835,484,295,198,153,87,68,55,29,18,18,7,6,3,3,4,1,1,0,5,0,0,0,0,0,0,1))
PcylinJames67 <- rep(c(seq(0,20,2),50),c(81,29,24,11,5,3,2,1,0,1,1,1))
RechinHodasi69 <- rep(c(seq(0,48,1),140),c(16,10,11,5,4,8,5,2,4,4,1,3,3,1,1,1,4,0,0,1,2,3,0,1,4,0,1,1,0,0,0,2,1,0,0,0,1,0,1,0,0,2,0,1,1,0,0,1,1,1))
RtetraHodasi69 <- rep(c(seq(0,21,1),28,29,33,47,72,75,81),c(45,8,8,2,6,4,3,1,1,4,1,3,1,1,1,1,3,1,1,1,2,1,2,2,1,1,1,1,1))
SjaponLiHsu51 <- rep(c(seq(0,19,1),25,26,42,44,49,109),c(603,28,10,4,3,4,2,3,2,0,1,3,2,2,2,0,1,4,2,1,6,1,2,1,2,1))
TfissiHodasi69 <- rep(c(seq(0,13,1),20,21,32,52),c(73,6,4,5,1,2,2,2,3,0,3,0,0,2,2,1,1,1))
TnodperchChubb64 <- rep(seq(0,5,1),c(211,189,61,23,8,5))
TnodpikeChubb63 <- rep(seq(0,43,1),c(24,4,5,1,2,6,4,5,3,5,6,3,3,1,4,5,3,0,1,3,0,2,0,0,1,1,1,1,2,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1))
CoxyStromberg75 <- rep(seq(0,6,1),c(103,71,43,14,3,2,2))
AdujarJames67 <- rep(seq(0,6,1),c(123,18,9,4,4,0,2))
AcaninLiHsu51 <- rep(c(seq(0,23,1),33,39,43),c(105,27,15,6,12,5,4,1,0,3,3,2,1,1,0,1,2,2,1,0,1,0,2,3,3,1,1))


#First data is summarised to give the appropriate information on
#mean, variance, host and parasite number for the comparison of
#these real distributions with those that are estimated.

data<-sort(AtetraLiHsu51,decreasing=TRUE)# copy the vector of interest from above to here and change number col to fit data)
	data
	sampmean<-mean(data)
	sampmean               ###3.446
		a<-numeric(10000)
		for(i in 1:10000) a[i]<-mean(sample(data,replace=T),na.rm=T)
		quantile(a,c(0.025,0.975))
				###2.589 - 4.545
	sumNpar<-sum(data)
	sumNpar
	N<-length(data)
	N
	var<-var(data)
	var
	k<-(sampmean^2-(var/N))/(var-sampmean)
	k


############################################
###Generate the distribution of parasites ##
###in a given proportion of the hosts     ##  
###Using the real data distributions      ##
############################################    

real2<-data.frame(data)

##distribution of parasites in the sample population
real3<-data.frame(sort(real2[,1],decreasing=TRUE))
real3

samp<-length(real3[,1])
nmk<-numeric(samp)
mk<-matrix(data = NA, nrow = samp, ncol = 1) 


for (j in 1:samp) # a nested loop which counts through the number of hosts

{


            addup<-sum(real3[1:j,1])             # cumulative sum of parasites 
     if (sum(real3)>0){nmk[j]<-addup/sum(real3)} # if there are ANY parasites this gives the cumulative proportion transmission per host

            else {nmk[j]<-0 }                    # if there are no parasites pop stops a problem due to division by 0

}

mk[,1]<-nmk                                      # puts nmk into matrix mk in column i


mk

raw<-data.frame(mk)
props<-c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)#the proportions we want
################################################################################

names<-paste("t",props,sep="")#names for each list of proportions
tvals<-matrix(NA,nrow=dim(raw)[2],ncol=length(props))# a matrix to store the tx values in - rows = studies, cols = tx values


for (j in 1:length(props))# first loop cycles through desired tx values contained in props above
{
	for (i in 1:dim(raw)[2]) #second loop cycles through each column of the data (which is a vector of proportions of parasites in each host from most to least infected)
	{
	vec<-raw[,i]             #make a vector from the column
	vec<-na.omit(vec)        #take out NAs after values reach 1 (it's a ragged array)
	N<-length(vec)           #number of hosts = number of rows
	Nx<-round(props[j]*N)    # gets the row which is x of the way down the column
	if (Nx==0) tvals[i,j]<-0 #if it rounds Nx to 0, there is no such thing as a zeroth row, so make tx 0
	else {tx<-vec[Nx]        # gets the value in that column - the proportion of parasites captured in the x most infected hosts
	tvals[i,j]<-tx           # rowi=study index colj=tx index}

		}
	}

}
t.vals<-cast(tvals)

##########
#data<-sort(Pkell.place[,1:2],decreasing=TRUE)# copy the data matrix here and change number col to fit data)
	data

t.vals[2,]<-t.vals[1,]*sum(data)
t.vals[2,]

t.vals<-cast(t.vals)

######################################################
###Step 2: Generate the estimated distributions from##
###the real data summary information; mean, N, var  ##
######################################################


data<-150    #This is the total number of parasites in study 'Aa'
dmean<-3.75  # e.g. This is the sample mean in study 'Aa'
samp<-40     # e.g. This is the number of hosts (N) in study 'Aa'
k<-0.202     # This is the k estimate from the mean var and N in study 'Aa'


mk<-matrix(data = NA, nrow = samp, ncol = 1000) # creates a matrix with
                                                # rows equal to sample size and columns equal to iteration number
tt20<-numeric(samp)                             # temporary vector for t10 equal to sample size
nmk<-numeric(samp)                              # storage vector also equal to sample size

for( i in 1:10000) # bootstrap starts (runs for 10000) iterations
	{ 

	dist1<-rnegbin (samp, dmean, k)      # creates a negbin distribution with your data
	dizzy<-sort(dist1,decreasing = TRUE) #sorts distribution highest parasite count first

		for (j in 1:samp) # a nested loop which counts through the number of hosts

			{

		            addup<-sum(dizzy[1:j])                      # cumulative sum of parasites 
     			    if (sum(dizzy)>0){nmk[j]<-addup/sum(dizzy)} # if there are ANY parasites this gives the cumulative proportion transmission per host
            		    else {nmk[j]<-0 }                           # if there are no parasites pop stops a problem due to division by 0

			}

	mk[,i]<-nmk # puts nmk into matrix mk in column i

	}

fk<-numeric(samp)

for (i in 1:samp)

	{fk[i]<-mean(mk[i,])} # take average of 10000 iterations for each row (each sample)

raw<-data.frame(fk)
props<-c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)
###The proportion of hosts to which we want to fit parasite populations

names<-paste("t",props,sep="")#names for each list of proportions
tvals<-matrix(NA,nrow=dim(raw)[2],ncol=length(props))# a matrix to store the tx values in - rows = studies, cols = tx values

for (j in 1:length(props))# first loop cycles through desired tx values contained in props above
	{
	for (i in 1:dim(raw)[2]) #second loop cycles through each column of the data (which is a vector of proportions of parasites in each host from most to least infected)
		{
		vec<-raw[,i]          #make a vector from the column
		vec<-na.omit(vec)     #take out NAs after values reach 1 (it's a ragged array)
		N<-length(vec)        #number of hosts = number of rows
		Nx<-round(props[j]*N) #gets the row which is x of the way down the column
		
		if (Nx==0) tvals[i,j]<-0 #if it rounds Nx to 0, there is no such thing as a zeroth row, so make tx 0
		else {tx<-vec[Nx]        # gets the value in that column - the proportion of parasites captured in the x most infected hosts
		tvals[i,j]<-tx           # rowi=study index colj=tx index}
			}
		}
	}

t.vals<-data.frame(tvals) #stick it in a dataframe and name the columns
#t.vals                   #check if you want

#data<-sort(Pkell.place[,1:2],decreasing=TRUE)# copy the data matrix here and change number col to fit data)
	data

t.vals[2,]<-t.vals[1,]*(data)
t.vals[2,]

t.vals2<-cast(t.vals)

####################################################################################
###Data is now ready to be compared                       ##########################
####################################################################################

##The summary data are all stored in "~\CRIPES 4\O5 Sept\Real vs estimate\\Summary data for real and estimated distributions.xlsx
dat<-read.csv("~\\05 Sept\\Real vs estimate\\t20comparison_v2")###t20COMPARISON_v2.csv

#AtetraLiHsu51 ##outlier - deleted
#HtaranShaw ##outlier - deleted

names(dat)

dim(dat)
head(dat)
stackT20<-c(dat$estimated0.2,dat$real0.2)
#stackT30<-c(dat$estimated0.3,dat$real0.3)
stackT40<-c(dat$estimated0.4,dat$real0.4)
#stackT50<-c(dat$estimated0.5,dat$real0.5)
stackT60<-c(dat$estimated0.6,dat$real0.6)
#stackT70<-c(dat$estimated0.7,dat$real0.7)
stackT80<-c(dat$estimated0.8,dat$real0.8)

length(stackT20)
datatype<-c(rep("estimate",30),rep("real",30))
mat<-data.frame(datatype,stackT20,stackT40,
	stackT60,stackT80)
###creating a data frame with the variables I want to investigate...

modt20<-glm(log(mat$stackT20)~mat$datatype,family=gaussian) ###Running a series of models to compare eal and estimated data...
summary.lm(modt20)
par(mfrow=c(2,2))
plot(modt20)
plot(modt20$resid)

modt40<-glm(log(mat$stackT40)~mat$datatype,family=gaussian)
summary.lm(modt40)
par(mfrow=c(2,2))
plot(modt40)

modt80<-glm(log(mat$stackT80)~mat$datatype,family=gaussian)
summary.lm(modt80)
par(mfrow=c(2,2))
plot(modt80)
boxplot(log(mat$stackT20)~mat$datatype,ylab="T20")
boxplot(log(mat$stackT40)~mat$datatype,ylab="T40")
boxplot(log(mat$stackT60)~mat$datatype,ylab="T60")
boxplot(log(mat$stackT80)~mat$datatype,,ylab="T80")

###############################################
####The above is conducted on the 'number'  ### 
####of parasites in a given proportion of   ###
####the host population... repeated with    ###
####proportions.                            ###
###############################################

##e.g.
STACKT20<-c(dat$EST0.2,dat$REAL0.2)###Not sure whether need to use proportion
		                   ###because here just want to compare two numbers...?
Nc<-c(dat$Nhosts,dat$Nhosts)
datatype<-c(rep("estimate",30),rep("real",30))
length(datatype)

mat2<-data.frame(datatype,STACKT20,Nc)
modT<-glm(STACKT20~datatype,family=quasibinomial,data=mat2,
	weights=Nc)
summary(modT)
plot(modT)
plot(modT$resid)
plot(modT$resid~mat2$datatype)

####################################################################
##Step 3: Generate the data set for the ultimate analysis ##########
####################################################################

#As above