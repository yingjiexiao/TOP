# The script was used to obtain the trait predictions of training set and test set using GBLUP. 
# Notes:
# dt: phenotype of training set
# A: relationship matrix of training set 
# A_train: relationship matrix of training set 
# A_test: relationship matrix 
# nut: the location of test set 
# num_test: the number of individuals in test set   
# trainop: phenotype of training set
# testop: phenotype of test set
# CEVR: the cumulative explained variance ratio
# Pre_ER10: the function was used to obtain the phenotype predictions of training set by 10-fold cross-validation (CV) 
# T_ER10: the function was used to obtain the phenotype predictions of test set

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


T_ER10<-function(trainop,A_test){
A <- A_test
length_test<-nrow(A_test)-length(trainop)
sop<-c(rep(1,length_test),trainop)
sop[1:length_test]<-NA
data <- data.frame(y=sop,gid=rownames(A))
ans <- kin.blup(data=data,geno="gid",pheno="y",K=A)
pl_te<-ans$g[1:length_test]+mean(trainop)
return(pl_te)
}

Prediction_trait<-function(Genotype,Agro,Metabolite,id,pca=TRUE,CEVR=0.8){
genotype_M<-Genotype[,2:ncol(Genotype)]
genotype<-t(genotype_M)

data_Met0<-Metabolite[,2:ncol(Metabolite)]
data_Met_M<-t(data_Met0)
data_Met<-apply(data_Met_M,2,function(x){
x[is.na(x)]=mean(x,na.rm=TRUE)
return(x)})

Agro_M<-Agro[,-1]
Agro1<-apply(Agro_M,2,function(x){
x[is.na(x)]=mean(x,na.rm=TRUE)
return(x)})


if (pca==TRUE){
PCA_Met<-prcomp(data_Met,scale = TRUE)
metabolite<-predict(PCA_Met)[,which(summary(PCA_Met)$
importance[3,]<=CEVR)] 
}else{
metabolite<-data_Met
}
trait_all<-cbind(Agro1,metabolite)
trait_all<-as.matrix(trait_all)

names_trait<-colnames(trait_all) 
test_id<-id[which(id[,2]=="test"),1]
train_id<-id[which(id[,2]=="training"),1]
index_test<-match(test_id,Agro[,1])
index_train<-match(train_id,Agro[,1])
num_test<-length(index_test) 
#the number of individuals in test set

trait_test<-trait_all[index_test,]
trait_train<-trait_all[index_train,]
names_test<-Agro[index_test,1]
#the names of individuals in test set
num_CV_test<-floor(nrow(trait_train)/10)
kinship0<-A.mat(as.matrix(genotype))
kinship<-kinship0[c(index_test,index_train),c(index_test,index_train)] 
K_train<-kinship[-c(1:num_test),-c(1:num_test)]
prediction_train<-matrix(0,nrow=nrow(trait_all)-num_test+11,
ncol=ncol(trait_all))
for (i in 1:ncol(trait_all)){
prediction_train[,i]<-Pre_ER10(trait_train[,i],K_train,num_CV_test)
}
trait_train<-apply(trait_train,2,as.numeric)
prediction_test<-matrix(0,nrow=num_test,ncol=ncol(trait_all))
K_full<-kinship
for (i in 1:ncol(trait_all)){
prediction_test[,i]<-T_ER10(trait_train[,i],K_full)
}
trait_test<-apply(trait_test,2,as.numeric)
model=list(prediction_train=prediction_train,prediction_test=prediction_test,trait_train=trait_train,trait_test=trait_test,names_test=names_test,names_trait=names_trait)
return(model)
}