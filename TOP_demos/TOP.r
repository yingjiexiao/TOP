# The script was used to construct TOP model, test TOP model, obtain the optimal weights and select similar individuals for a given target.
# Notes:
# N: pool sizes
# predicted_train0: the predicted values of the traits in the training set  
# predicted_test: the predicted values of the traits in the test set  
# observeded_train: the observeded values of the traits in the training set
# observeded_test: the observeded values of the traits in the test set
# W: the optimal weights
# Pre_test: the normalized predicted values of the traits in the test set   
# Obs_test: the normalized observeded values of the traits in the test set
# n: pool size
# m: the number of random permutations
# target: a given target
###########################TOP model construction and test TOP model
prePro<-function(data){
 Mean<-apply(data,2,mean)
 SD<-apply(data,2,sd)
 value<-apply(data,2,function(x){(x-mean(x))/sd(x)})
 model=list(Mean=Mean,SD=SD,value=value)
 return(model)
}
 pro<-function(w,a_M){
	pi_M<-NULL
	proi<-matrix(0,nrow=dim(a_M)[1],ncol=dim(a_M)[1])
	sf<-matrix(0,nrow=dim(a_M)[1],ncol=dim(a_M)[2])
	for(i in 1:dim(a_M)[1]){
		pi_M<-apply(a_M[,,i],1,function(x){w%*%x})
		proi[i,]<-exp(pi_M)/(sum(exp(pi_M)))
		sf[i,]<-proi[i,]%*%a_M[,,i]
	}
	model=list(proi=proi,sf=sf)
return(model)
}

ide_rate<-function(w,Obs_test,Pre_test,n,m){
rate<-NULL
for(k in 1:m){
set.seed(123)
ord<-sample(1:nrow(Obs_test))
Pre_test_M<-Pre_test[ord,]
Obs_test_M<-Obs_test[ord,]
p<-array(0,dim=c(n,n,dim(Pre_test_M)[1]/n))
for (j in 1:(dim(Pre_test_M)[1]/n)){
pre_test1<-as.matrix(Pre_test_M[(n*(j-1)+1):(n*j),])
obs_test1<-as.matrix(Obs_test_M[(n*(j-1)+1):(n*j),])
Dis<-abs(obs_test1-pre_test1)
m_M<-matrix(1,nrow=dim(obs_test1)[1],ncol=1)
a_M<-array(0,dim=c(dim(pre_test1)[1],dim(pre_test1)[2],dim(pre_test1)[1]))
for(i in 1:dim(pre_test1)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(pre_test1[i,])-obs_test1)
}
p[,,j]<-pro(w,a_M)$proi
}
a<-NULL
for(j in 1:dim(Pre_test_M)[1]/n){
	a[j]<-sum(apply(p[,,j],1,which.max)==c(1:n))
	}
rate[k]<-sum(a)/(floor(dim(Pre_test_M)[1]/n)*n)
}
return(mean(rate))
}

###############################################################Learn optimal weights

Weight_res<-function(predicted_train0,observeded_train,names_trait,b){

Ib<-which(predicted_train0[nrow(predicted_train0),]>=b)
names_trait_M<-names_trait[Ib]
predicted_train<-predicted_train0[1:nrow(observeded_train),]

pre_Pro<-prePro(predicted_train[,Ib])
obs_Pro<-prePro(observeded_train[,Ib])

predicted_train_M<-pre_Pro$value
observeded_train_M<-obs_Pro$value

Pre_train<-as.matrix(predicted_train_M)
Obs_train<-as.matrix(observeded_train_M)
Dis<-abs(Pre_train-Obs_train)

m_M<-matrix(1,nrow=dim(Obs_train)[1],ncol=1)
a_M<-array(0,dim=c(dim(Pre_train)[1],dim(Pre_train)[2],dim(Pre_train)[1]))


for(i in 1:dim(Pre_train)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(Pre_train[i,])-Obs_train)
}

fr<-function(w){
	pi_M<-NULL
	proi<-matrix(0,nrow=dim(a_M)[1],ncol=dim(a_M)[1])
	sf<-matrix(0,nrow=dim(a_M)[1],ncol=dim(a_M)[2])
	for(i in 1:dim(a_M)[1]){
		pi_M<-apply(a_M[,,i],1,function(x){w%*%x})
		proi[i,]<-exp(pi_M)/(sum(exp(pi_M)))
		sf[i,]<-proi[i,]%*%a_M[,,i]
	}
	f<-sum(log(diag(proi)))
return(f)
}
gr<-function(w){
	pi_M<-NULL
	proi<-matrix(0,nrow=dim(a_M)[1],ncol=dim(a_M)[1])
	sf<-matrix(0,nrow=dim(a_M)[1],ncol=dim(a_M)[2])
	for(i in 1:dim(a_M)[1]){
		pi_M<-apply(a_M[,,i],1,function(x){w%*%x})
		proi[i,]<-exp(pi_M)/(sum(exp(pi_M)))
		sf[i,]<-proi[i,]%*%a_M[,,i]
	}
	g<-apply(Dis,2,sum)-apply(sf,2,sum)
return(g)
}
res<-constrOptim(rep(0,ncol(Obs_train)),
fr,gr,method="BFGS",control=list(fnscale=-1,maxit=1500),ui =rbind(rep(1,ncol(Obs_train))),ci = rep(-20,1))
w<-res$par
W<-(-w)
W_matrix<-data.frame(trait=names_trait_M,weight=W)
model=list(W_matrix=W_matrix,pre_mean=pre_Pro$Mean,pre_sd=pre_Pro$SD,obs_mean=obs_Pro$Mean,obs_sd=obs_Pro$SD)
return(model)
}

#################################################################Test the top model
Test_top_acc<-function(predicted_test,observeded_test,W,names_trait,N,m){ 

Ib<-match(W[,1],names_trait)
predicted_test_M<-prePro(predicted_test[,Ib])$value
observeded_test_M<-prePro(observeded_test[,Ib])$value

w<-(-W[,2])
Ide_rate<-NULL
for (i in 1:length(N)){
Ide_rate[i]<-ide_rate(w,predicted_test_M,observeded_test_M,N[i],m)
}
Demo_p<-Demo_P(w,predicted_test_M,observeded_test_M,N[2],m)
model=list(Ide_rate=Ide_rate,Demo_p=Demo_p)
return(model)
}
Demo_P<-function(w,Obs_test,Pre_test,n,m){
for(k in 1:m){
set.seed(123)
ord<-sample(1:nrow(Obs_test))
Pre_test_M<-Pre_test[ord,]
Obs_test_M<-Obs_test[ord,]
p<-array(0,dim=c(n,n,dim(Pre_test_M)[1]/n))
for (j in 1:(dim(Pre_test_M)[1]/n)){
pre_test1<-as.matrix(Pre_test_M[(n*(j-1)+1):(n*j),])
obs_test1<-as.matrix(Obs_test_M[(n*(j-1)+1):(n*j),])
Dis<-abs(obs_test1-pre_test1)
m_M<-matrix(1,nrow=dim(obs_test1)[1],ncol=1)
a_M<-array(0,dim=c(dim(pre_test1)[1],dim(pre_test1)[2],dim(pre_test1)[1]))
for(i in 1:dim(pre_test1)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(pre_test1[i,])-obs_test1)
}
p[,,j]<-pro(w,a_M)$proi
}
}
similarity_matrix<-p[,,1]
rownames(similarity_matrix)=colnames(similarity_matrix)=paste('id',1:5,sep='')
return(similarity_matrix)
}

############################################################Select similar individuals for a given target

Tpro<-function(w,ta_M){
	pi_M<-NULL
	tproi<-matrix(0,nrow=1,ncol=dim(ta_M)[1])
	tpi_M<-apply(ta_M,1,function(x){w%*%x})
	tproi<-exp(tpi_M)/sum(exp(tpi_M))
return(tproi)
}
Top_target<-function(target0,predicted_test,names_test,names_trait,pre_train_mean,pre_train_sd,obs_train_mean,obs_train_sd,selection_ratio,improve_ratio,improve_trait,W){
Ib<-match(W[,1],names_trait)
target1<-target0
if (improve_trait!="NA"){
target1[improve_trait]<-(1+improve_ratio)*target0[improve_trait]
target<-target1[Ib]
}else{
target<-target1[Ib]
}
w<-(-W[,2])
Pre_test<-predicted_test[,Ib]
target_M<-(target-obs_train_mean)/obs_train_sd
predicted_test_M<-(Pre_test-matrix(1,nrow=nrow(Pre_test),ncol=1)%*%pre_train_mean)/matrix(1,nrow=nrow(Pre_test),ncol=1)%*%pre_train_sd
Target<-as.matrix(rep(1,nrow(predicted_test_M)))%*%t(as.matrix(target_M))
ta_M<-abs(Target-predicted_test_M)
similarity0<-Tpro(w,ta_M)
select_order<-match(sort(similarity0,decreasing=TRUE),similarity0)
select<-select_order[1:(nrow(Pre_test)*selection_ratio)]
names_select<-names_test[select]
values_select<-Pre_test[select,]
model=list(names_select=names_select,values_select=values_select)
return(model)
}

plot_TOP<-function(target,trait_test,Top_target_sel_name,improve_trait=4,group_random=5){
value_t<-trait_test[match(Top_target_sel_name,names_test),]
value_r<-matrix(0,nrow=(nrow(value_t)*group_random),ncol=ncol(value_t))
value_r_imp<-NULL
set.seed(123)
for (i in 1:group_random){
value_r[((i-1)*nrow(value_t)+1):(i*nrow(value_t)),]<-trait_test[sample(1:nrow(trait_test))[1:length(Top_target_sel_name)],]
value_r_imp<-c(value_r_imp,value_r[((i-1)*nrow(value_t)+1):(i*nrow(value_t)),improve_trait])
}
Cdata<-c(value_t[,improve_trait],value_r_imp)
CL<-c(rep("TOP",nrow(value_t)),rep("Random",(nrow(value_t)*group_random)))
Select_mode<-factor(CL,levels=c("TOP","Random"))
mydata<-data.frame(CL,Cdata)
p1=ggplot(mydata,aes(x=Cdata,fill=Select_mode))+
  geom_density(alpha=.3)+
  ylab('Density')+
  scale_x_continuous(limits=c(18,30),breaks=seq(18,30,2))+
  scale_y_continuous(limits=c(0,0.4),breaks=seq(0,0.4,0.1))+
  xlab(paste("Trait",improve_trait,sep=""))+scale_fill_manual(values = c("#F8766D","#619CFF"))+
  guides(fill=guide_legend(title=NULL))+
   theme_bw()+
  theme(legend.title=element_blank(),legend.position=c(0.1, 0.95),legend.text = element_text(size=8),
legend.key.size=unit(0.3,'cm'),legend.key = element_blank(),legend.background = element_blank())

 ggsave("./res/demo_density.png")
CMatrix<-rbind(target,value_t,value_r)
CMS<-apply(CMatrix,2,function(x){(x-min(x))/(max(x)-min(x))})
SM0<-as.matrix(rep(1,(nrow(CMS)-1)))%*%as.numeric(CMS[1,])
SM<-CMS[2:nrow(CMS),]-SM0
MSE1<-apply(SM[1:nrow(value_t),-improve_trait],1,function(x){mean(x^2)})
MSE2<-apply(SM[(nrow(value_t)+1):nrow(SM),-improve_trait],1,function(x){mean(x^2)})
Cdata1<-c(MSE1,MSE2)
mydata1<-data.frame(CL,Cdata1)
p2=ggplot(mydata1,aes(x=Select_mode,y=Cdata1,fill=Select_mode))+ geom_boxplot(alpha=.3) +
scale_fill_manual(values = c("#F8766D","#619CFF"))+
  ylab(paste('MSE without Trait',improve_trait,sep="") )+ 
  xlab("")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw()+
  theme(legend.position="none")
  ggsave("./res/demo_MSE_boxplot.png")
 model=list(p1=p1,p2=p2)
 return(model)
}




