# The script was used to construct TOP model, test TOP model, obtain the optimal weights and select similar individuals for a given target.
# Notes:
# N: pool sizes
# res_M: the trait predictions of training set  
# Tres_M: the trait predictions of test set  
# dis: the real phenotypes of training set 
# Tdis: the real phenotypes of test set
# w: the optimal weights
# Tres_M0: the normalized trait predictions of test set   
# Tdis0: the normalized real phenotypes of test set
# n: pool size
# m: the number of random permutations
# target: a given target

###########################TOP model construction and test TOP model
cd<-function(data){
 m<-apply(data,2,mean)
 s<-apply(data,2,sd)
 nop<-apply(data,2,function(x){(x-mean(x))/sd(x)})
 model=list(m=m,s=s,nop=nop)
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

mc<-function(w,Tdis0,Tres_M0,n,m){
sp<-NULL
for(k in 1:m){
so<-sample(1:nrow(Tdis0))
Tres_M<-Tres_M0[so,]
Tdis0N<-Tdis0[so,]
p<-array(0,dim=c(n,n,dim(Tres_M)[1]/n))
for (j in 1:(dim(Tres_M)[1]/n)){
res_Pl<-as.matrix(Tres_M[(n*(j-1)+1):(n*j),])
dt<-as.matrix(Tdis0N[(n*(j-1)+1):(n*j),])
Dis<-abs(dt-res_Pl)
m_M<-matrix(1,nrow=dim(dt)[1],ncol=1)
a_M<-array(0,dim=c(dim(res_Pl)[1],dim(res_Pl)[2],dim(res_Pl)[1]))
for(i in 1:dim(res_Pl)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(res_Pl[i,])-dt)
}
p[,,j]<-pro(w,a_M)$proi
}
a<-NULL
for(j in 1:dim(Tres_M)[1]/n){
	a[j]<-sum(apply(p[,,j],1,which.max)==c(1:n))
	}
sp[k]<-sum(a)/(floor(dim(Tres_M)[1]/n)*n)
}
return(mean(sp))
}


###############################################################Learn optimal weights
Weight_res<-function(res_M,dis){

res_M_M<-res_M[1:nrow(dis),]
dis_M<-dis

dis0<-cd(dis_M)$nop
res_M0<-cd(res_M_M)$nop


res_Pl<-as.matrix(res_M0)
dt<-as.matrix(dis0)
Dis<-abs(dt-res_Pl)
m_M<-matrix(1,nrow=dim(dt)[1],ncol=1)
a_M<-array(0,dim=c(dim(res_Pl)[1],dim(res_Pl)[2],dim(res_Pl)[1]))

for(i in 1:dim(res_Pl)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(res_Pl[i,])-dt)
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

res<-constrOptim(rep(0,ncol(dis)),
fr,gr,method="BFGS",control=list(fnscale=-1,maxit=1500),ui =rbind(rep(1,ncol(dis))),ci = rep(-20,1))
ww<-res$par
W<-(-ww)
return(W)
}
#################################################################Test the top model
Test_top_res<-function(W,res_M,Tres_M,dis,Tdis,N,m){
res_M_M<-res_M[1:nrow(dis),]
dis_M<-dis
Tres_M_M<-Tres_M[1:nrow(Tdis),]
Tdis_M<-Tdis
dis0<-cd(dis_M)$nop
res_M0<-cd(res_M_M)$nop

Tres_M0<-cd(Tres_M_M)$nop
Tdis0<-cd(Tdis_M)$nop

w<-(-W)
MC<-NULL
for (i in 1:length(N)){
MC[i]<-mc(w,Tdis0,Tres_M0,N[i],m)
}
Demo_p<-Demo_P(w,Tdis0,Tres_M0,N[2],m)
model=list(MC=MC,Demo_p=Demo_p)
return(model)
}

Demo_P<-function(w,Tdis0,Tres_M0,n,m){
for(k in 1:m){
so<-sample(1:nrow(Tdis0))
Tres_M<-Tres_M0[so,]
Tdis0N<-Tdis0[so,]
p<-array(0,dim=c(n,n,dim(Tres_M)[1]/n))
for (j in 1:(dim(Tres_M)[1]/n)){
res_Pl<-as.matrix(Tres_M[(n*(j-1)+1):(n*j),])
dt<-as.matrix(Tdis0N[(n*(j-1)+1):(n*j),])
Dis<-abs(dt-res_Pl)
m_M<-matrix(1,nrow=dim(dt)[1],ncol=1)
a_M<-array(0,dim=c(dim(res_Pl)[1],dim(res_Pl)[2],dim(res_Pl)[1]))
for(i in 1:dim(res_Pl)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(res_Pl[i,])-dt)
}
p[,,j]<-pro(w,a_M)$proi
}
}
return(p[,,1])
}

Demo_P<-function(w,Tdis0,Tres_M0,n,m){
for(k in 1:m){
so<-sample(1:nrow(Tdis0))
Tres_M<-Tres_M0[so,]
Tdis0N<-Tdis0[so,]
p<-array(0,dim=c(n,n,dim(Tres_M)[1]/n))
for (j in 1:(dim(Tres_M)[1]/n)){
res_Pl<-as.matrix(Tres_M[(n*(j-1)+1):(n*j),])
dt<-as.matrix(Tdis0N[(n*(j-1)+1):(n*j),])
Dis<-abs(dt-res_Pl)
m_M<-matrix(1,nrow=dim(dt)[1],ncol=1)
a_M<-array(0,dim=c(dim(res_Pl)[1],dim(res_Pl)[2],dim(res_Pl)[1]))
for(i in 1:dim(res_Pl)[1]){
    a_M[,,i]<-abs(m_M%*%as.numeric(res_Pl[i,])-dt)
}
p[,,j]<-pro(w,a_M)$proi
}
}
return(p[,,1])
}
############################################################Select similar individuals for a given target
tpro<-function(w,ta_M){
	pi_M<-NULL
	tproi<-matrix(0,nrow=1,ncol=dim(ta_M)[1])
	tpi_M<-apply(ta_M,1,function(x){w%*%x})
	tproi<-exp(tpi_M)/sum(exp(tpi_M))
return(tproi)
}

TrainMres<-function(ww,Tres_M,res_M,dis,target){
w<-(-ww)
Tres_M_M<-Tres_M[1:nrow(Tdis),]
res_M_M<-res_M[1:nrow(dis),]

dis_M<-dis
target_M<-target


target0<-(target_M-cd(dis_M)$m)/cd(dis_M)$s
Tres_M0<-(Tres_M_M-matrix(1,nrow=nrow(Tres_M_M),ncol=1)%*%cd(res_M_M)$m)/matrix(1,nrow=nrow(Tres_M_M),ncol=1)%*%cd(res_M_M)$s

res_Pl<-as.matrix(rep(1,nrow(Tres_M0)))%*%t(as.matrix(target0))

ta_M<-abs(res_Pl-Tres_M0)
tpa<-tpro(w,ta_M)

otpav<-sort(tpa,decreasing=TRUE)
otpa<-match(sort(tpa,decreasing=TRUE),tpa)
return(otpa)
}