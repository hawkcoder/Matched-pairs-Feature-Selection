
#Given a set of covariates (matrix mtrain) and response vector utrain
#this function fit linear regression and select the best fitting covariate
fwd.lm1=function(mtrain,utrain){ #begin
	#Input parameters: mtrain for protein measurements in training
	#                  vtrain for binary response (0/1 group) in training 

	do.lm <- function(cont.var,bin.var)
	{	
		data.to.fit <- as.data.frame(cbind(bin.var=bin.var,cont.var=cont.var))
		res <-  summary(lm(bin.var ~ cont.var, data = data.to.fit))
		z <- res$coeff[2,4]
    	z
	}

	## Now to perform all the linear regressions:

	logreg.z <- apply(mtrain, 2, FUN=do.lm, utrain)  #This instruction is a substitution of a "while" loop.

	# Get the best proteins in the forward logistic regr and keep the rest of the proteins in rest.proteins

	ind = which.min(logreg.z)
	min.pvalue = min(logreg.z)
	best.proteins=mtrain[,ind[1]]
	best.proteins.numbers=ind[1]

	#run glm on best proteins set
	lm.best =lm(utrain~best.proteins)
	list(best.proteins.numbers=best.proteins.numbers,lm.best=lm.best,min.pvalue=min.pvalue)
}   #end

#This function fits weighted linear regression
#given a covariate f, offset F, response vector v and cluster identifier id (for generating weights)
fwd.lme=function(f,v,F,id){ #begin
  ## mi are individual estimates, sei are their standard errors
	y = v
	Z = model.matrix(~factor(id)-1)
	X = cbind(F,f)

  	loss<-function (par,y){
		tau.b = par[1]
		tau = par[2]
		V = (exp(tau.b)*Z%*%t(Z) + exp(tau)*diag(dim(Z)[1]))
		L = chol(V,pivot=TRUE)
		y = backsolve(L,y,transpose=TRUE)
		X = backsolve(L,X,transpose=TRUE)
		lm.reg = lm(y~offset(X[,1])+X[,2]-1)
		b = list(coef=coef(lm.reg),pvalue=summary(lm.reg)$coef[1,4])
		r = (y-X[,1]-b$coef*X[,2])
		lo = sum((r)^2)
		attr(lo,"fixed") <- b #allow retrieval of beta
	lo
  	}
    	res <- lm(y~offset(X[,1])+X[,2]-1)
	out <- nlminb(c(-10,log(sum(res$residuals^2)/res$df.residual)),loss,lower=c(-10,-10),upper=c(3,3),y=y)
	if (out$convergence == 0) {
		sigmab = sqrt(exp(out$par[1]))
		sigma = sqrt(exp(out$par[2]))
		b = attr(loss(c(out$par[1],out$par[2]),y),"fixed")
		list(coef=c(b$coef[1]),pvalue=b$pvalue,sigma=sigma,sigma.b=sigmab,loss=loss(c(out$par[1],out$par[2]),y))
	}
    	else {
		sigmab = sqrt(exp(out$par[1]))
		sigma = sqrt(exp(out$par[2]))
		b = attr(loss(c(out$par[1],out$par[2]),y),"fixed")
		list(coef=c(b$coef[1]),pvalue=1,sigma=sigma,sigma.b=sigmab,loss=loss(c(out$par[1],out$par[2]),y))
    	}
}

#This function performs WL2Boost given the training dataset (matrix of covariates mtrain, response vector vtrain, cluster identifier id
#										 optional test dataset)
boostingWL2=function(mtrain,vtrain,id,stop,mtest=NULL,vtest=NULL,idtest=NULL){ #begin
	#Input parameters: mtrain/mtest for protein measurements in training/test set
	#                  vtrain/vtest for binary response (-1/1 group) in training/test set 
	#                  stop is the number of iteration
	#                  if mtest/vtest are missing, mtrain/vtrain are used
	if(is.null(mtest)){mtest=mtrain;vtest=vtrain;idtest=id}
	iter=0
	proteins=NULL
	coef=matrix(NA,stop,1)
	pred=rep(0,length(vtrain))
	test.pred=rep(0,length(vtest))
	mat.pred = NULL
	mat.test.pred = NULL
	Ftrain = rep(0,length(vtrain))
	Ftest = rep(0,length(vtest))

	utrain = vtrain
	Ystar = vtrain
	dev = NULL
	aic = NULL
	infl = NULL
	sigma.b = 1
	while(iter<stop){
		iter=iter+1
		res1=fwd.lm1(mtrain,utrain)
		proteins.num = res1$best.proteins.numbers
		f = as.vector(predict(res1$lm.best))
		f.test = as.vector(as.vector(res1$lm.best$coef)%*%rbind(rep(1,length(vtest)),mtest[,as.numeric(proteins.num)]))
		res = fwd.lme(f=f,v=Ystar,F=Ftrain,id)
		coeff = as.vector(res$coef)
		sigma = res$sigma
		sigma.b = res$sigma.b
		Z = model.matrix(~factor(id)-1)
		V = (sigma.b^2)*Z%*%t(Z) + diag(nrow(Z))*sigma^2
		Vinv = solve(V)
		Ftrain = Ftrain + c*coeff*f #pred1
		Ftest = as.vector(Ftest) + as.vector(c*coeff*f.test) #test.pred1
		Ystar = vtrain 
		utrain = Vinv%*%(vtrain-Ftrain)
		proteins[iter]=as.numeric(proteins.num)
		coef[iter]=as.vector(coeff) #res1$glm.best.proteins$coeff
		infl = c(infl,coef[iter]*sd(f))
		mat.pred = cbind(mat.pred,Ftrain)
		mat.test.pred = cbind(mat.test.pred,Ftest)
		cat("\n Iteration: ",iter)
	}
	mat.infl = cbind(1:stop,proteins,infl)
	list(proteins=proteins,coef=coef,vtrain=vtrain,mat.pred=mat.pred,vtest=vtest,mat.test.pred=mat.test.pred,dev=dev,aic=aic,mat.infl=mat.infl)
}

#5-fold cross-validation
boosting.CV.5=function(newdata,stop,seed=12345){
	set.seed(seed)	
	id1=unique(newdata[,1])#newdata[,dim(newdata)[2]]==1]
	id1.5=rep(1:5,ceiling(length(id1)/5))[1:length(id1)]
	id1.5=sample(id1.5)
	for(i in 1:5){
		cat("\n", "crossvalidation run:", i)
		a = sort(c(id1[-which(id1.5==i)]))
		b = sort(c(id1[which(id1.5==i)]))
		m.temp=newdata[which(newdata[,1]%in%a),2:22284]
		v.temp=newdata[which(newdata[,1]%in%a),22285]
		id.temp=newdata[which(newdata[,1]%in%a),1]
		m.tempt=newdata[which(newdata[,1]%in%b),2:22284]
		v.tempt=newdata[which(newdata[,1]%in%b),22285]
		id.tempt=newdata[which(newdata[,1]%in%b),1] 
		res=boostingWL2(m.temp,v.temp,id.temp,stop,m.tempt,v.tempt,id.tempt)
		if(i==1){res1=res}
		if(i==2){res2=res}
		if(i==3){res3=res}
		if(i==4){res4=res}
		if(i==5){res5=res}
	
	}
	list(cv1=res1,cv2=res2,cv3=res3,cv4=res4,cv5=res5)
}

#Data processing and implementation begins

Response=read.csv("time.csv", sep=",",header=T)
Response = as.data.frame(Response)
dim(Response)
time = Response[,5]

new=read.table("mydata.txt",header=T)
time.long = NULL
for (i in 1:35) {
	time.long = c(time.long,rep(time[i],2))
}
relapse = rep(c(-1,1),35)

newdata = as.data.frame(cbind(relapse,t(new[,1:70]),time.long))


#Subsetting the data 
early_late = 1 #Set this to zero to select early relapse cases, one for late relapse cases, otherwise entire data is desired
if (early_late==0) {
	newdata = newdata[which(newdata[,22285]<36),] #Early relapse subset
} else if (early_late==1) {
	newdata = newdata[-which(newdata[,22285]<36),] #Late relapse subset
} else {}

mtrain = newdata[,2:22284]
vtrain = newdata[,1]
id = NULL 
for (i in 1:(length(vtrain)/2)) {
	id = c(id,rep(i,2))
}
newdata = as.data.frame(cbind(id,mtrain,vtrain))
stop = 100
CV.boost = boosting.CV.5(newdata,stop,seed=12345)
niter.cv = min(dim(CV.boost$cv1$mat.test.pred)[2],dim(CV.boost$cv2$mat.test.pred)[2],dim(CV.boost$cv3$mat.test.pred)[2],dim(CV.boost$cv4$mat.test.pred)[2],dim(CV.boost$cv5$mat.test.pred)[2],dim(CV.boost$cv6$mat.test.pred)[2],dim(CV.boost$cv7$mat.test.pred)[2],dim(CV.boost$cv8$mat.test.pred)[2],dim(CV.boost$cv9$mat.test.pred)[2],dim(CV.boost$cv10$mat.test.pred)[2])
t.error = NULL
c.error = NULL
for (i in 1:niter.cv){
	cvpred = c(CV.boost$cv1$mat.test.pred[,i],CV.boost$cv2$mat.test.pred[,i],CV.boost$cv3$mat.test.pred[,i],CV.boost$cv4$mat.test.pred[,i],CV.boost$cv5$mat.test.pred[,i],CV.boost$cv6$mat.test.pred[,i],CV.boost$cv7$mat.test.pred[,i],CV.boost$cv8$mat.test.pred[,i],CV.boost$cv9$mat.test.pred[,i],CV.boost$cv10$mat.test.pred[,i])
	cvtest = c(CV.boost$cv1$vtest,CV.boost$cv2$vtest,CV.boost$cv3$vtest,CV.boost$cv4$vtest,CV.boost$cv5$vtest,CV.boost$cv6$vtest,CV.boost$cv7$vtest,CV.boost$cv8$vtest,CV.boost$cv9$vtest,CV.boost$cv10$vtest)
	post = cvpred[2*(1:(length(cvpred)/2))]
	pre = cvpred[2*(1:(length(cvpred)/2))-1]
	bin.test.pred = NULL
	for (j in 1:length(pre)) {
		if ((max(abs(pre[j]),abs(post[j]))==-pre[j])|(max(abs(pre[j]),abs(post[j]))==post[j])) {a=c(-1,1)}
		else {a=c(1,-1)}  
		bin.test.pred = c(bin.test.pred,a)
	}
	cv.error = length(which(cvtest*bin.test.pred<0))/length(cvtest)
	c.error = c(c.error,cv.error)
} 
plot(1:stop,c.error,type="l",lty=2,col = "blue",ylim=c(0,.5),xlab="Number of Iterations",ylab="Classification Error")
typ.names <- c("CV Error","Training Error") 
legend("topright", legend = typ.names,lty = c(2,1),col = c("blue","red"))

boost.object = boosting(mtrain,vtrain,id,stop)
old.sink.level <- sink("WL2_Aug15_late.txt",append=TRUE) 
on.exit(sink(unsink.to=old.sink.level))
cat("\n Early Relapse is denoted by 0 and Late Relapse is denoted by 1, all otherwise: ", early_late,"\n")
A <- matrix(nrow = stop,ncol=3,NA)
names <- c("# Iterations", "Training Error", "CV Error")
names2 <- c("Protein #", "Protein name", "Beta", "Coef", "Infl")
proteins.num = boost.object$protein
proteins = colnames(mtrain[,boost.object$protein])

coef = round(boost.object$coef,2)
A2 <- matrix(nrow =stop,ncol=4,NA)
A2 <- cbind(proteins.num,proteins,coef,boost.object$mat.infl[,3])
A2 <- as.data.frame(A2)
dimnames(A2) <- list(paste(seq(1:stop)),paste(names2))
print(A2)
#cat("\n The Predictions (Using margin magnitude and paired information): \n")
#print(as.data.frame(boost.object$mat.pred))

t.error.margin = NULL
for (i in 1:stop){
	pred = boost.object$mat.pred[,i]
	vtrain = boost.object$vtrain
	post = pred[2*(1:(length(vtrain)/2))]
	pre = pred[2*(1:(length(vtrain)/2))-1]
	bin.pred = NULL
	for (j in 1:length(pre)) {
		if ((max(abs(pre[j]),abs(post[j]))==-pre[j])|(max(abs(pre[j]),abs(post[j]))==post[j])) {a=c(-1,1)}
		else {a=c(1,-1)} #elseif ((store[i]==pre[i])|(store[i]==-post[i])) 
		bin.pred = c(bin.pred,a)
	}
	train.error.margin = length(which(vtrain*bin.pred<0))/length(vtrain)
	t.error.margin = c(t.error.margin,train.error.margin)
	A[i,] = cbind(i,train.error.margin,c.error[i])
} 
A <- as.data.frame(A)
dimnames(A) <- list(paste(c(1:stop)),paste(names))
print(A)
sink()
lines(1:stop,t.error.margin,lty=1,col = "red")

Response=read.csv("time.csv", sep=",",header=T)
Response = as.data.frame(Response)
dim(Response)
time = Response[,5]

new=read.table("mydata.txt",header=T)
#new=as.data.frame(new)
#id = NULL 
time.long = NULL
for (i in 1:35) {
	#id = c(id,rep(i,2))
	time.long = c(time.long,rep(time[i],2))
}
relapse = rep(c(-1,1),35)
newdata = as.data.frame(cbind(relapse,t(new[,1:70]),time.long))


#Subsetting the data 
early_late = 1 #Set this to zero to select early relapse cases, one for late relapse cases, otherwise entire data is desired
if (early_late==0) {
	newdata = newdata[which(newdata[,22285]<36),] #Early relapse subset
} else if (early_late==1) {
	newdata = newdata[-which(newdata[,22285]<36),] #Late relapse subset
} else {}

mtrain = newdata[,2:22284]
vtrain = newdata[,1]
id = NULL 
for (i in 1:(length(vtrain)/2)) {
	id = c(id,rep(i,2))
}
newdata = as.data.frame(cbind(id,mtrain,vtrain))
stop = 100
CV.boost = boosting.CV.5(newdata,stop,seed=12345)
niter.cv = min(dim(CV.boost$cv1$mat.test.pred)[2],dim(CV.boost$cv2$mat.test.pred)[2],dim(CV.boost$cv3$mat.test.pred)[2],dim(CV.boost$cv4$mat.test.pred)[2],dim(CV.boost$cv5$mat.test.pred)[2],dim(CV.boost$cv6$mat.test.pred)[2],dim(CV.boost$cv7$mat.test.pred)[2],dim(CV.boost$cv8$mat.test.pred)[2],dim(CV.boost$cv9$mat.test.pred)[2],dim(CV.boost$cv10$mat.test.pred)[2])
t.error = NULL
c.error = NULL
for (i in 1:niter.cv){
	cvpred = c(CV.boost$cv1$mat.test.pred[,i],CV.boost$cv2$mat.test.pred[,i],CV.boost$cv3$mat.test.pred[,i],CV.boost$cv4$mat.test.pred[,i],CV.boost$cv5$mat.test.pred[,i],CV.boost$cv6$mat.test.pred[,i],CV.boost$cv7$mat.test.pred[,i],CV.boost$cv8$mat.test.pred[,i],CV.boost$cv9$mat.test.pred[,i],CV.boost$cv10$mat.test.pred[,i])
	cvtest = c(CV.boost$cv1$vtest,CV.boost$cv2$vtest,CV.boost$cv3$vtest,CV.boost$cv4$vtest,CV.boost$cv5$vtest,CV.boost$cv6$vtest,CV.boost$cv7$vtest,CV.boost$cv8$vtest,CV.boost$cv9$vtest,CV.boost$cv10$vtest)
	post = cvpred[2*(1:(length(cvpred)/2))]
	pre = cvpred[2*(1:(length(cvpred)/2))-1]
	bin.test.pred = NULL
	for (j in 1:length(pre)) {
		if ((max(abs(pre[j]),abs(post[j]))==-pre[j])|(max(abs(pre[j]),abs(post[j]))==post[j])) {a=c(-1,1)}
		else {a=c(1,-1)} #elseif ((store[i]==pre[i])|(store[i]==-post[i])) 
		bin.test.pred = c(bin.test.pred,a)
	}
	cv.error = length(which(cvtest*bin.test.pred<0))/length(cvtest)
	c.error = c(c.error,cv.error)
} 
plot(1:stop,c.error,type="l",lty=2,col = "blue",ylim=c(0,.5),xlab="Number of Iterations",ylab="Classification Error")
typ.names <- c("CV Error","Training Error") 
legend("topright", legend = typ.names,lty = c(2,1),col = c("blue","red"))

#stop = which.min(c.error.margin)
boost.object = boosting(mtrain,vtrain,id,stop)
old.sink.level <- sink("WL2_Aug09_late.txt",append=TRUE) 
on.exit(sink(unsink.to=old.sink.level))
cat("\n Early Relapse is denoted by 0 and Late Relapse is denoted by 1, all otherwise: ", early_late,"\n")
A <- matrix(nrow = stop,ncol=3,NA)
names <- c("# Iterations", "Training Error", "CV Error")
names2 <- c("Protein #", "Protein name", "Beta", "Infl")
proteins.num = boost.object$protein
proteins = colnames(mtrain[,boost.object$protein])

coef = round(boost.object$coef,2)
A2 <- matrix(nrow =stop,ncol=4,NA)
A2 <- cbind(proteins.num,proteins,coef,boost.object$mat.infl[,3])
A2 <- as.data.frame(A2)
dimnames(A2) <- list(paste(seq(1:stop)),paste(names2))
print(A2)
#cat("\n The Predictions (Using margin magnitude and paired information): \n")
#print(as.data.frame(boost.object$mat.pred))

t.error.margin = NULL
for (i in 1:stop){
	pred = boost.object$mat.pred[,i]
	vtrain = boost.object$vtrain
	post = pred[2*(1:(length(vtrain)/2))]
	pre = pred[2*(1:(length(vtrain)/2))-1]
	bin.pred = NULL
	for (j in 1:length(pre)) {
		if ((max(abs(pre[j]),abs(post[j]))==-pre[j])|(max(abs(pre[j]),abs(post[j]))==post[j])) {a=c(-1,1)}
		else {a=c(1,-1)} #elseif ((store[i]==pre[i])|(store[i]==-post[i])) 
		bin.pred = c(bin.pred,a)
	}
	train.error.margin = length(which(vtrain*bin.pred<0))/length(vtrain)
	t.error.margin = c(t.error.margin,train.error.margin)
	A[i,] = cbind(i,train.error.margin,c.error[i])
} 
A <- as.data.frame(A)
dimnames(A) <- list(paste(c(1:stop)),paste(names))
print(A)
sink()
lines(1:stop,t.error.margin,lty=1,col = "red")

Response=read.csv("time.csv", sep=",",header=T)
Response = as.data.frame(Response)
dim(Response)
time = Response[,5]

new=read.table("mydata.txt",header=T)
#new=as.data.frame(new)
#id = NULL 
time.long = NULL
for (i in 1:35) {
	#id = c(id,rep(i,2))
	time.long = c(time.long,rep(time[i],2))
}
relapse = rep(c(-1,1),35)
newdata = as.data.frame(cbind(relapse,t(new[,1:70]),time.long))


#Subsetting the data 
early_late = 2 #Set this to zero to select early relapse cases, one for late relapse cases, otherwise entire data is desired
if (early_late==0) {
	newdata = newdata[which(newdata[,22285]<36),] #Early relapse subset
} else if (early_late==1) {
	newdata = newdata[-which(newdata[,22285]<36),] #Late relapse subset
} else {}

mtrain = newdata[,2:22284]
vtrain = newdata[,1]
id = NULL 
for (i in 1:(length(vtrain)/2)) {
	id = c(id,rep(i,2))
}
newdata = as.data.frame(cbind(id,mtrain,vtrain))
stop = 130
CV.boost = boosting.CV.5(newdata,stop,seed=12345)
niter.cv = min(dim(CV.boost$cv1$mat.test.pred)[2],dim(CV.boost$cv2$mat.test.pred)[2],dim(CV.boost$cv3$mat.test.pred)[2],dim(CV.boost$cv4$mat.test.pred)[2],dim(CV.boost$cv5$mat.test.pred)[2],dim(CV.boost$cv6$mat.test.pred)[2],dim(CV.boost$cv7$mat.test.pred)[2],dim(CV.boost$cv8$mat.test.pred)[2],dim(CV.boost$cv9$mat.test.pred)[2],dim(CV.boost$cv10$mat.test.pred)[2])
t.error = NULL
c.error = NULL
for (i in 1:niter.cv){
	cvpred = c(CV.boost$cv1$mat.test.pred[,i],CV.boost$cv2$mat.test.pred[,i],CV.boost$cv3$mat.test.pred[,i],CV.boost$cv4$mat.test.pred[,i],CV.boost$cv5$mat.test.pred[,i],CV.boost$cv6$mat.test.pred[,i],CV.boost$cv7$mat.test.pred[,i],CV.boost$cv8$mat.test.pred[,i],CV.boost$cv9$mat.test.pred[,i],CV.boost$cv10$mat.test.pred[,i])
	cvtest = c(CV.boost$cv1$vtest,CV.boost$cv2$vtest,CV.boost$cv3$vtest,CV.boost$cv4$vtest,CV.boost$cv5$vtest,CV.boost$cv6$vtest,CV.boost$cv7$vtest,CV.boost$cv8$vtest,CV.boost$cv9$vtest,CV.boost$cv10$vtest)
	post = cvpred[2*(1:(length(cvpred)/2))]
	pre = cvpred[2*(1:(length(cvpred)/2))-1]
	bin.test.pred = NULL
	for (j in 1:length(pre)) {
		if ((max(abs(pre[j]),abs(post[j]))==-pre[j])|(max(abs(pre[j]),abs(post[j]))==post[j])) {a=c(-1,1)}
		else {a=c(1,-1)} #elseif ((store[i]==pre[i])|(store[i]==-post[i])) 
		bin.test.pred = c(bin.test.pred,a)
	}
	cv.error = length(which(cvtest*bin.test.pred<0))/length(cvtest)
	c.error = c(c.error,cv.error)
} 
plot(1:stop,c.error,type="l",lty=2,col = "blue",ylim=c(0,.5),xlab="Number of Iterations",ylab="Classification Error")
typ.names <- c("CV Error","Training Error") 
legend("topright", legend = typ.names,lty = c(2,1),col = c("blue","red"))

#stop = which.min(c.error.margin)
boost.object = boosting(mtrain,vtrain,id,stop)
old.sink.level <- sink("WL2_Aug09_all.txt",append=TRUE) 
on.exit(sink(unsink.to=old.sink.level))
cat("\n Early Relapse is denoted by 0 and Late Relapse is denoted by 1, all otherwise: ", early_late,"\n")
A <- matrix(nrow = stop,ncol=3,NA)
names <- c("# Iterations", "Training Error", "CV Error")
names2 <- c("Protein #", "Protein name", "Beta", "Var")
proteins.num = boost.object$protein
proteins = colnames(mtrain[,boost.object$protein])

coef = round(boost.object$coef,2)
A2 <- matrix(nrow =stop,ncol=4,NA)
A2 <- cbind(proteins.num,proteins,coef,boost.object$mat.infl[,3])
A2 <- as.data.frame(A2)
dimnames(A2) <- list(paste(seq(1:stop)),paste(names2))
print(A2)
#cat("\n The Predictions (Using margin magnitude and paired information): \n")
#print(as.data.frame(boost.object$mat.pred))

t.error.margin = NULL
for (i in 1:stop){
	pred = boost.object$mat.pred[,i]
	vtrain = boost.object$vtrain
	post = pred[2*(1:(length(vtrain)/2))]
	pre = pred[2*(1:(length(vtrain)/2))-1]
	bin.pred = NULL
	for (j in 1:length(pre)) {
		if ((max(abs(pre[j]),abs(post[j]))==-pre[j])|(max(abs(pre[j]),abs(post[j]))==post[j])) {a=c(-1,1)}
		else {a=c(1,-1)} #elseif ((store[i]==pre[i])|(store[i]==-post[i])) 
		bin.pred = c(bin.pred,a)
	}
	train.error.margin = length(which(vtrain*bin.pred<0))/length(vtrain)
	t.error.margin = c(t.error.margin,train.error.margin)
	A[i,] = cbind(i,train.error.margin,c.error[i])
} 
A <- as.data.frame(A)
dimnames(A) <- list(paste(c(1:stop)),paste(names))
print(A)
sink()
lines(1:stop,t.error.margin,lty=1,col = "red")
