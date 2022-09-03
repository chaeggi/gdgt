######## GDGT Monte Carlo Models #############################
rm(list=ls())
install.packages("matrixStats")
library(matrixStats)
setwd("C:/Users/ch/Desktop/Cerrado Data")

######## Load data#############################

brisoGDGT<-read.csv("globalbrisoGDGT.csv")

######## BIT index #############################

################### Model Parameters and set up

#Set seed value for pseudorandom number generation; allows for reproduction of the output values
set.seed(1)
N=nrow(brisoGDGT) #number of sites
S=5 #number of samples per site
variability<- 0.5 #variablity added

#Data for each used compound

Cren<-brisoGDGT[,5] 
BrGDGTIa<-brisoGDGT[,7] 
BrGDGTIIa<-brisoGDGT[,10] 
BrGDGTIIIa<-brisoGDGT[,13] 

#Create vector for the added variability of each modelled compound

Crenvar<-vector(length=N, mode="numeric")
BrGDGTIavar<-vector(length=N, mode="numeric")
BrGDGTIIavar<-vector(length=N, mode="numeric")
BrGDGTIIIavar<-vector(length=N, mode="numeric")

#Create vector for output

BIT<-vector(length=N, mode="numeric")
results<-matrix(nrow = N,ncol = S)
output<-matrix(nrow = N,ncol = 2)

############# Monte Carlo Model


#Run model for S samples per sites

for (i in 1:S){
  for (j in 1:N) { #add variablilty to each compound at each site
    Crenvar[j]<-rnorm(1,Cren[j],Cren[j]*variability)
    BrGDGTIavar[j]<-rnorm(1,   BrGDGTIa[j],   BrGDGTIa[j]*variability)
    BrGDGTIIavar[j]<-rnorm(1,   BrGDGTIIa[j],   BrGDGTIIa[j]*variability)
    BrGDGTIIIavar[j]<-rnorm(1,   BrGDGTIIIa[j],   BrGDGTIIIa[j]*variability)
    BIT[j]<-(BrGDGTIavar[j]+BrGDGTIIavar[j]+BrGDGTIIIavar[j])/(BrGDGTIavar[j]+BrGDGTIIavar[j]+BrGDGTIIIavar[j]+Crenvar[j])
      if (BIT[j] < 0 | BIT[j]> 1) {BIT[j]<-NA } ####The setup allows for a small probability of negative concentrations impossible in nature, which are removed here
  }
  results[,i]<-BIT

}
# Export mean and standard deviation of the data. 
mean<-rowMeans(results)
stdev<-rowSds(results)


output<-data.frame(mean,stdev)


write.csv(output,file=sprintf("C:/Users/ch/Documents/ModelBIT.csv"), row.names=FALSE)
BITor<-vector(length=N, mode="numeric")
BITor<-(BrGDGTIa+BrGDGTIIa+BrGDGTIIIa)/(BrGDGTIa+BrGDGTIIa+BrGDGTIIIa+Cren)
plot(BIT,BITor)

######## Ri/b #############################


#Set seed value for pseudorandom number generation; allows for reproduction of the output values
set.seed(1)


#Define the number of sites for the MC modelling and the number of samples for each site
N=nrow(brisoGDGT) #number of sites
S=5 #number of samples
variability= 0.5 #variablity added

#Data for each used compound
Ia<-brisoGDGT[,7] 
Ib<-brisoGDGT[,8] 
Ic<-brisoGDGT[,9] 
IIa<-brisoGDGT[,10] 
IIb<-brisoGDGT[,11] 
IIc<-brisoGDGT[,12] 
IIIa<-brisoGDGT[,13] 
IIIb<-brisoGDGT[,14] 
IIIc<-brisoGDGT[,15] 


isoGDGT0<-brisoGDGT[,1] 
isoGDGT1<-brisoGDGT[,2] 
isoGDGT2<-brisoGDGT[,3] 
isoGDGT3<-brisoGDGT[,4] 
Cren<-brisoGDGT[,5] 
Crenis<-brisoGDGT[,6] 

#Prepare empty vectors to add variability for each of the modelled comonents
Iavar<-vector(length=N, mode="numeric")
Ibvar<-vector(length=N, mode="numeric")
Icvar<-vector(length=N, mode="numeric")
IIavar<-vector(length=N, mode="numeric")
IIbvar<-vector(length=N, mode="numeric")
IIcvar<-vector(length=N, mode="numeric")
IIIavar<-vector(length=N, mode="numeric")
IIIbvar<-vector(length=N, mode="numeric")
IIIcvar<-vector(length=N, mode="numeric")

isoGDGT0var<-vector(length=N, mode="numeric")
isoGDGT1var<-vector(length=N, mode="numeric")
isoGDGT2var<-vector(length=N, mode="numeric")
isoGDGT3var<-vector(length=N, mode="numeric")
Crenvar<-vector(length=N, mode="numeric")
Crenisvar<-vector(length=N, mode="numeric")

#Create vector for output
Rib<-vector(length=N, mode="numeric")
results<-matrix(nrow = N,ncol = S)
output<-matrix(nrow = N,ncol = 2)




##### Monte Carlo Model

for (i in 1:S){
  for (j in 1:N) {
    
    Iavar[j]<-rnorm(1,Ia[j],Ia[j]*variability)
    Ibvar[j]<-rnorm(1,Ib[j],Ib[j]*variability)
    Ibvar[j]<-rnorm(1,Ic[j],Ic[j]*variability)
    IIavar[j]<-rnorm(1,IIa[j],IIa[j]*variability)
    IIbvar[j]<-rnorm(1,IIb[j],IIb[j]*variability)
    IIcvar[j]<-rnorm(1,IIc[j],IIc[j]*variability)
    IIIavar[j]<-rnorm(1,IIIa[j],IIIa[j]*variability)
    IIIbvar[j]<-rnorm(1,IIIb[j],IIIb[j]*variability)
    IIIcvar[j]<-rnorm(1,IIIc[j],IIIc[j]*variability)
    isoGDGT0var[j]<-rnorm(1,isoGDGT0[j],isoGDGT0[j]*variability)
    isoGDGT1var[j]<-rnorm(1,isoGDGT1[j],isoGDGT1[j]*variability)
    isoGDGT2var[j]<-rnorm(1,isoGDGT2[j],isoGDGT2[j]*variability)
    isoGDGT3var[j]<-rnorm(1,isoGDGT3[j],isoGDGT3[j]*variability)
    Crenvar[j]<-rnorm(1,Cren[j],Cren[j]*variability)
    Crenisvar[j]<-rnorm(1,Crenis[j],Crenis[j]*variability)
    Rib[j]<-(isoGDGT0var[j]+isoGDGT1var[j]+isoGDGT2var[j]+isoGDGT3var[j]+Crenvar[j]+Crenisvar[j])/(Iavar[j]+Ibvar[j]+Icvar[j]+IIavar[j]+IIcvar[j]+IIIavar[j]+IIIbvar[j]+IIIcvar[j])
    if (  Rib[j] < 0) {  Rib[j]<-NA } 
  }
  results[,i]<-Rib
}

######## Export mean and standard deviation of the data. 
mean<-rowMeans(results)
stdev<-rowSds(results)

output<-data.frame(mean,stdev)


write.csv(output,file=sprintf("C:/Users/ch/Documents/ModelRib.csv"), row.names=FALSE)


######## MBT5Me index #############################
brGDGT<-read.csv("globalbrGDGT.csv")

#Set seed value for pseudorandom number generation; allows for reproduction of the output values
set.seed(1)

#Define the number of sites for the MC modelling and the number of samples for each site
N=nrow(brGDGT) #number of sites
S=5 #number of samples
variability= 0.25 #variablity added


#Create vector for each modelled compound
Ia<-brGDGT[,1] 
Ib<-brGDGT[,2] 
Ic<-brGDGT[,3] 
IIa<-brGDGT[,4] 
IIb<-brGDGT[,6] 
IIc<-brGDGT[,8] 
IIIa<-brGDGT[,10] 


#Prepare empty vectors to add variability for each of the modelled comonents
Iavar<-vector(length=N, mode="numeric")
Ibvar<-vector(length=N, mode="numeric")
Icvar<-vector(length=N, mode="numeric")
IIavar<-vector(length=N, mode="numeric")
IIbvar<-vector(length=N, mode="numeric")
IIcvar<-vector(length=N, mode="numeric")
IIIavar<-vector(length=N, mode="numeric")

#Create vector for output
MBT5<-vector(length=N, mode="numeric")
results<-matrix(nrow = N,ncol = S)
output<-matrix(nrow = N,ncol = 2)

#####



for (i in 1:S){
  for (j in 1:N) {
    
    Iavar[j]<-rnorm(1,Ia[j],Ia[j]*variability)
    Ibvar[j]<-rnorm(1,Ib[j],Ib[j]*variability)
    Ibvar[j]<-rnorm(1,Ic[j],Ic[j]*variability)
    IIavar[j]<-rnorm(1,IIa[j],IIa[j]*variability)
    IIbvar[j]<-rnorm(1,IIb[j],IIb[j]*variability)
    IIcvar[j]<-rnorm(1,IIc[j],IIc[j]*variability)
    IIIavar[j]<-rnorm(1,IIIa[j],IIIa[j]*variability)
    MBT5[j]<-(Iavar[j]+Ibvar[j]+Icvar[j])/(Iavar[j]+Ibvar[j]+Icvar[j]+IIavar[j]+IIcvar[j]+IIIavar[j])
    if (MBT5[j] < 0 | MBT5[j]> 1) {MBT5[j]<-NA }
  }
  results[,i]<-MBT5
}

######## Export mean and standard deviation of the data. 
mean<-rowMeans(results)
stdev<-rowSds(results)

output<-data.frame(mean,stdev)


write.csv(output,file=sprintf("C:/Users/ch/Documents/ModelMBT.csv"), row.names=FALSE)

plot(mean,stdev)

