#TOP_Demo.r

rm(list=ls())
library(rrBLUP)
# load the script to obtain the trait predictions of training set and testing set
source("GP.r")
# load the script to construct the top model and obtain the optimal weights
source("TOP.r")

# read genotype 
Riceg<-read.csv("RIL-Genotype.csv",header=TRUE,sep=",",stringsAsFactors=0,check.names=FALSE)
dataS_M<-Riceg[,7:ncol(Riceg)]
data_S<-t(dataS_M)

# read phenotype
ptrait<-read.table("data_trait.txt",header=TRUE,sep="",stringsAsFactors=0,check.names=FALSE)
ptrait<-as.matrix(ptrait)


num_ex<-floor(nrow(ptrait)/2)  #the number of individuals in testing set
N<-c(2,5,10)                   #pool size
m<-50                          #the number of random permutations
b<-0.25                        #the prediction accuracy threshold of entering the TOP model

# obtain the trait predictions of training set and testing set using GBLUP 
A<-A.mat(as.matrix(data_S))    
ord<-sample(1:nrow(ptrait))

dt_ex<-as.matrix(ptrait[ord,])

testop1_exs<-dt_ex[1:num_ex,]
trainop1_exs<-dt_ex[-c(1:num_ex),]

AA<-A[ord,ord]

A_train<-AA[-c(1:num_ex),-c(1:num_ex)]
A_test<-AA

res_M1<-matrix(0,nrow=nrow(dt_ex)-num_ex+11,ncol=ncol(dt_ex))
Tres_M1<-matrix(0,nrow=num_ex+1,ncol=ncol(dt_ex))

for (i in 1:ncol(dt_ex)){
testop1_ex<-dt_ex[1:num_ex,i]
trainop1_ex<-dt_ex[-c(1:num_ex),i]

dt<-trainop1_ex
num_test<-floor(length(dt)/10)
res_M1[,i]<-Pre_ER10(dt,A_train,num_test)

Tdt<-testop1_ex
Tres_M1[,i]<-T_ER10(dt,Tdt,A_test)
}
dis1<-apply(dt_ex[-c(1:num_ex),],2,as.numeric)
Tdis1<-apply(dt_ex[1:num_ex,],2,as.numeric)


res_M<-res_M1              #the trait predictions of training set  
Tres_M<-Tres_M1            #the trait predictions of testing set  
dis<-dis1                  #the real phenotypes of training set 
Tdis<-Tdis1                #the real phenotypes of testing set

# construct TOP model and test the model
Ib<-which(res_M[nrow(res_M),]>=b)
wres<-wres(res_M[,Ib],Tres_M[,Ib],dis[,Ib],Tdis[,Ib],N,m)

write.table(wres$MC,file="Demo_reidentification.txt") #genomic reidentification of pool size N=(2,5,10)

# obtain the optimal weights
W<-T_wres(res_M[,Ib],dis[,Ib])
write.table(W,file="Demo_W.txt")

# the top 10 selection from a pool of 105 individuals 
target<-dis[1,]
Order_Selected<-TrainMres(W,Tres_M[,Ib],res_M[,Ib],dis[,Ib],target[Ib])[1:10]
write.table(Order_Selected,file="Order_Selected.txt")
######################################################################
