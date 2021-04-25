# the script to obtain the trait predictions of training set and testing set using GBLUP 
# Notes:
# dt: phenotype of training set
# A: relationship matrix of training set 
# A_train: relationship matrix of training set 
# A_test: relationship matrix 
# nut: the location of testing set 
# num_test: the number of individuals in testing set   
# trainop: phenotype of training set
# testop: phenotype of testing set


library(rrBLUP)

PR<-function(dt,A,nut){
dt_M<-dt
dt_M[nut]<-NA
data <- data.frame(y=dt_M,gid=rownames(A))
ans <- kin.blup(data=data,geno="gid",pheno="y",K=A)
pl_te<-ans$g[nut]+mean(dt[-nut])
pcor_te<-cor(as.numeric(dt[nut]),pl_te)
Model=list(pcor_te=pcor_te,pl_te=pl_te) 
return(Model)
}

Pre_ER10<-function(dt,A_train,num_test){
A <-A_train
res1<-PR(dt,A,1:num_test)
res2<-PR(dt,A,(num_test+1):(num_test*2))
res3<-PR(dt,A,(num_test*2+1):(num_test*3))
res4<-PR(dt,A,(num_test*3+1):(num_test*4))
res5<-PR(dt,A,(num_test*4+1):(num_test*5))
res6<-PR(dt,A,(num_test*5+1):(num_test*6))
res7<-PR(dt,A,(num_test*6+1):(num_test*7))
res8<-PR(dt,A,(num_test*7+1):(num_test*8))
res9<-PR(dt,A,(num_test*8+1):(num_test*9))
res10<-PR(dt,A,(num_test*9+1):nrow(A))
res_Pcor<-c(res1$pcor_te,res2$pcor_te,res3$pcor_te,res4$pcor_te,res5$pcor_te,res6$pcor_te,res7$pcor_te,res8$pcor_te,res9$pcor_te,res10$pcor_te)
res_Pl<-c(res1$pl_te,res2$pl_te,res3$pl_te,res4$pl_te,res5$pl_te,res6$pl_te,res7$pl_te,res8$pl_te,res9$pl_te,res10$pl_te)
Pcor10<-sum(res_Pcor)/10
res=c(res_Pl,res_Pcor,Pcor10)
return(res)
}

T_ER10<-function(trainop,testop,A_test){
A <- A_test
sop<-c(testop,trainop)
sop[1:length(testop)]<-NA
data <- data.frame(y=sop,gid=rownames(A))
ans <- kin.blup(data=data,geno="gid",pheno="y",K=A)
pl_te<-ans$g[1:length(testop)]+mean(trainop)
pcor_te<-cor(as.numeric(testop),pl_te)
res=c(pl_te,pcor_te)
return(res)
}
