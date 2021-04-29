# TOP_demo.R

# The script was used to build target-orient prioritization (TOP) model, obtain the optimal trait weights and test the accuracy of the TOP model.

# to run this TOP algrithm, three files are required to be prepared as the format of the demo data:
# 1. demo_geno.csv: a demo genotypic data with 1919 bin markers and 210 individuals from a rice RIL population. The genotypic value needed to be transformed to numeric data (0/1/2).
# 2. demo_pheno.txt: a demo phenotypic data with 210 individuals and 4 traits. The individual identifier (row name) needs to be matched with that of genotypic data.
# 3. demo_training_test_id.txt: a information file with 211 rows and 2 columns. The first column is the individual identifiers and the second indicates the use of training or test set.

# Several parameters and custom functions are demonstrated as followed:
# num_test: the number of individuals in test set
# N: pool size
# m: the number of random permutations
# b: the prediction accuracy threshold of entering the TOP model
# Pre_ER10: the function was used to obtain the phenotype predictions of training set by 10-fold cross-validation (CV) 
# T_ER10: the function was used to obtain the phenotype predictions of test set
# Test_top_res: the function was used to test the Top model and obtain identification rate of different pool size

rm(list=ls())
library(rrBLUP)
# load the script to obtain the trait predictions of training set and test set
source("GP.r")
# load the script to build the top model and obtain the optimal trait weights for calculate the genome-phenome similarity between testing pool and a given target
source("TOP.r")

# read genotype 
Rice_genotype<-read.csv("./data/demo_geno.csv",header=TRUE,sep=",",stringsAsFactors=0,check.names=FALSE)
genotype_M<-Rice_genotype[,2:ncol(Rice_genotype)]
genotype<-t(genotype_M)

# read phenotype
phenotype<-read.table("./data/demo_pheno.txt",header=TRUE,sep="",stringsAsFactors=0,check.names=FALSE)
phenotype<-as.matrix(phenotype)

N<-c(2,5,10)                        #pool size
m<-200                              #the number of random permutations
b<-0.25                             #the prediction accuracy threshold of entering the TOP model

#divide a population into training set and test set 

id<-read.table("./data/demo_training_test_id.txt",header=TRUE,sep="",stringsAsFactors=0,check.names=FALSE)
test_id<-id[which(id[,2]=="test"),1]
train_id<-id[which(id[,2]=="training"),1]
index_test<-match(test_id,rownames(phenotype))
index_train<-match(train_id,rownames(phenotype))
num_test<-length(index_test)        #the number of individuals in test set
phe_test<-phenotype[index_test,]
phe_train<-phenotype[index_train,]

#In training set, GBLUP was used to obtain the phenotype predictions of training set by 10-fold cross-validation (CV).
num_CV_test<-floor(nrow(phe_train)/10)
kinship0<-A.mat(as.matrix(genotype))
kinship<-kinship0[c(index_test,index_train),c(index_test,index_train)] 
K_train<-kinship[-c(1:num_test),-c(1:num_test)]
prediction_train<-matrix(0,nrow=nrow(phenotype)-num_test+11,ncol=ncol(phenotype))
for (i in 1:ncol(phenotype)){
prediction_train[,i]<-Pre_ER10(phe_train[,i],K_train,num_CV_test)
}
phe_train<-apply(phe_train,2,as.numeric)

# construct TOP model and obtain the optimal weights
Ib<-which(prediction_train[nrow(prediction_train),]>=b)
Weight<-Weight_res(prediction_train[,Ib],phe_train[,Ib])
write.table(data.frame(traitNAME=colnames(phenotype),weight=Weight),file="./res/demo_optimal_weight.txt",sep='\t',row.names=F,quote=F)

#obtain the trait predictions of test set using GBLUP 
prediction_test<-matrix(0,nrow=num_test,ncol=ncol(phenotype))
K_full<-kinship
for (i in 1:ncol(phenotype)){
prediction_test[,i]<-T_ER10(phe_train[,i],K_full)
}
phe_test<-apply(phe_test,2,as.numeric)

# test the top model
ide_rate<-Test_top_res(Weight,prediction_train[,Ib],prediction_test[,Ib],phe_train[,Ib],phe_test[,Ib],N,m)

# at the pool of 5 individuals, we output a demo result showing the similarity degree between genomic-predicted values and observed traits by testing each of 5 individuals as a target. 
similarity_matrix=ide_rate$Demo_p
rownames(similarity_matrix)=colnames(similarity_matrix)=paste('ind',1:5,sep='')

#in the similarity matrix of pool size of 5, each column means a given target individual from the pool, the 5 values for this column indicate the genome-phenome similarity between all individuals from the pool and this target individual. If the target individual itself has the highest similarity degree, we regard it as a successful identification. 
write.table(similarity_matrix,file="./res/demo_similarity_matrix_poolsize5.txt",sep='\t',quote=F) 

#at the pool size of 2,5 and 10, we output a demo result of identification rate, which means the proportion of successfull identification by using each individal of a pool as a target.
write.table(data.frame(pool_size=c(2,5,10),identification_rate=ide_rate$MC),file="./res/demo_identification_rate.txt",sep='\t',row.names=F,quote=F) 


