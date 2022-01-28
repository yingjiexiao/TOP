# TOP_demo.R

# The script was used to build target-orient prioritization (TOP) model, obtain the optimal trait weights and test the accuracy of the TOP model.

# to run this TOP algrithm, three files are required to be prepared as the format of the demo data:
# 1. demo_geno.csv: a demo genotypic data with 1919 bin markers and 210 individuals from a rice RIL population. The genotypic value needed to be transformed to numeric data (0/1/2).
# 2. demo_agro.txt: a demo phenotypic data with 210 individuals and 4 traits. The individual identifier needs to be matched with that of genotypic data.
# 3. demo_Metabolites.csv: a demo metabolic data with 210 individuals and 1000 metabolites. The individual identifier needs to be matched with that of genotypic data.
# 4. demo_training_test_id.txt: a information file with 211 rows and 2 columns. The first column is the individual identifiers and the second indicates the use of training or test set.

# Several parameters and custom functions are demonstrated as followed:
# N: pool size
# m: the number of random permutations
# b: the prediction accuracy threshold of entering the TOP model
# Test_top_acc: the function was used to test the Top model and obtain identification rate of different pool size

#The latest update: Jan.27 2022. Wenyu Yang (yangwenyurain@126.com)

rm(list=ls())
library(rrBLUP)
# load the script to obtain the trait predictions of training set and test set
source("GP.r")
# load the script to build the top model and obtain the optimal trait weights for calculate the genome-phenome similarity between testing pool and a given target
source("TOP.r")
library(ggplot2)
# read genotype 
Genotype<-read.csv("./data/demo_geno.csv",
                       header=TRUE,
                       sep=",",
				       stringsAsFactors=0,
				       check.names=FALSE
					   )
# read agronomic trait
Agro<-read.table("./data/demo_agro.txt",
                       header=TRUE,
                       sep="",
				       stringsAsFactors=0,
				       check.names=FALSE
					   )
# read Metabolites
Metabolite<-read.table("./data/demo_Metabolites.csv",
                       header=TRUE,
					   sep=",",
					   stringsAsFactors=0,
					   check.names=FALSE
					   )

# divide a population into training set and test set 
id<-read.table("./data/demo_training_test_id.txt",
                      header=TRUE,
                      sep="",
					  stringsAsFactors=0,
					  check.names=FALSE)

# In training set, GBLUP was used to obtain the predicted trait values of the training set by 10-fold cross-validation (CV).
Pre<-Prediction_trait(Genotype,
                      Agro,
					  Metabolite,
					  id,
					  pca=TRUE,
					  CEVR=0.8
					  )
# Genotype: genotypic data of the population.
# Agro: agronomic data of the population.
# Metabolite: metabolic data of the population.
# id: the names of individuals in the training and the test set.
# pca: if the parameter “pca” is default set to be TRUE, PCA is performed to reduce the dimension of metabolic data. 
# CEVR: cumulative explained variance of principal components. 

prediction_train<-Pre$prediction_train
trait_train<-Pre$trait_train #the observed trait values of the training set

# obtain the predicted trait values of the test set using GBLUP 
prediction_test<-Pre$prediction_test
trait_test<-Pre$trait_test #the observed trait values of the test set

# construct TOP model and obtain the optimal weights


names_test<-Pre$names_test                   #the names of the individuals in test set
names_trait<-Pre$names_trait                 #the names of the traits
target<-trait_train[80,]                     #a given target


Optimal_Weight<-Weight_res(prediction_train,
                           trait_train,
						   names_trait,
						   b=0.25)                                
# prediction_train: the predicted trait values of the training set and the last row stores the average prediction accuracy of ten round.
# trait_train: the observed trait values of the training set.
# names_trait: the names of the traits.
# b: the prediction accuracy threshold of entering the TOP model, the default is 0.25.

# test the top model
Weight<-Optimal_Weight$W_matrix
TOP_acc<-Test_top_acc(prediction_test,
                      trait_test,
					  Weight,
					  names_trait,
					  N=c(2,5,10),
					  m=200)
#prediction_test:the predicted trait values of the test set
#trait_test:the observed trait values of the test set
#Weight:indicate the importance of individual traits in maximizing the similarity
#N:pool size
#m:the number of random permutations

# prediction_test: the predicted trait values of the test set.
# trait_test: the observed trait values of the test set.
# Weight: indicate the importance of individual traits for searching global similarity to a target individual, the default is the optimal weight calculated by the TOP model.
# names_trait: the names of the traits.
# N: pool size for calculating TOP accuracy, the default is c(2,5,10).
# m: the number of random permutations when selecting a pool of individuals from the test set, the default is 200.

# select individuals in the test set according to the similarity function
pre_train_mean<-Optimal_Weight$pre_mean
pre_train_sd<-Optimal_Weight$pre_sd
obs_train_mean<-Optimal_Weight$obs_mean
obs_train_sd<-Optimal_Weight$obs_sd

Top_target_sel<-Top_target(target,
                           prediction_test,
						   names_test,
						   names_trait,
						   pre_train_mean,
						   pre_train_sd,
						   obs_train_mean,
						   obs_train_sd,
						   selection_ratio=0.1,
						   improve_ratio=0.05,
						   improve_trait=4,
                           Weight)					   						   
# target: a vector (1xd, d is the trait number) containing the values of each trait of a given target, often set as the trait observations of an existing commercial cultivar.
# prediction_test: the predicted trait values of the test set.
# names_test: the names of individuals in test set
# names_trait: the names of the traits
# pre_train_mean: the mean of predicted trait values of training set.
# pre_train_sd: the standard deviation of predicted trait values of training set.
# obs_train_mean: the mean of observed trait values of training set.
# obs_train_sd: the standard deviation of observed trait values of training set.
# selection_ratio: identify ratio of individuals in the test set with maximized similarity to the target, the default is 0.1, it indicates the top 10% similarity individuals will be selected. 
# improve_ratio: change ratio of a trait value from the target, the default is 5%, it indicates the improved target at this trait expect to have the phenotype 5% higher than the original one. 
# improve_trait: the trait ID of the target to be improved, if the ID is “NA”, no trait value of target will be modified; if the ID is “4”, the fourth trait in the list of all used traits will be modified. 
# Weight: indicate the importance of individual traits for searching elite candidates, the default is the optimal weights learned from TOP model.


write.table(Weight,file="./res/demo_optimal_weight.txt",sep='\t',row.names=F,quote=F)
# In the similarity matrix of pool size of 5, each column means a given target individual from the pool, the 5 values for this column indicate the genome-phenome similarity between all individuals from the pool and this target individual. If the target individual itself has the highest similarity degree, we regard it as a successful identification. 
write.table(TOP_acc$Demo_p,file="./res/demo_similarity_matrix_poolsize5.txt",sep='\t',quote=F) 

# At the pool size of 2,5 and 10, we output a demo result of identification rate, which means the proportion of successfull identification by using each individal of a pool as a target.
write.table(data.frame(pool_size=c(2,5,10),identification_rate=TOP_acc$Ide_rate),file="./res/demo_identification_rate.txt",sep='\t',row.names=F,quote=F) 

# Let the first individual in the training set be a target, we select individuals in the test set according to the similarity function. 
Top_sel<-cbind(Top_target_sel$names_select,Top_target_sel$values_select)
colnames(Top_sel)<- c("RIL",as.character(Weight[,1]))
write.table(Top_sel,file="./res/demo_target_select.txt",sep='\t',row.names=F,quote=F) 

Top_target_sel_name<-Top_target_sel$names_select
plot_TOP(target,
         trait_test,
		 Top_target_sel_name,
		 improve_trait=4,
		 group_random=5)
