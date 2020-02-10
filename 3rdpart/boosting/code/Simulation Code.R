library(MASS)

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


#This function performs PQLBoost given the training dataset (matrix of covariates mtrain, response vector vtrain, cluster identifier id
#										 optional test dataset)
boostingPQL=function(mtrain,vtrain,id,stop,mtest=NULL,vtest=NULL,idtest=NULL){ #begin
	#Input parameters: mtrain/mtest for protein measurements in training/test set
	#                  vtrain/vtest for binary response (-1/1 group) in training/test set 
	#                  stop is the number of iteration
	#                  if mtest/vtest are missing, mtrain/vtrain are used
	#                  w is weights
	if(is.null(mtest)){mtest=mtrain;vtest=vtrain;idtest=id}
	iter=0
	proteins=NULL
	coef=matrix(NA,stop,2)
	pred=rep(0,length(vtrain))
	test.pred=rep(0,length(vtest))
	mat.pred = NULL
	mat.test.pred = NULL
	Ftrain = rep(0,length(vtrain))
	Ftest = rep(0,length(vtest))
	dev = NULL
	aic = NULL
	infl = NULL
	while(iter<stop){
		res1=fwd.lmer1(mtrain,(vtrain+1)/2,id,Ftrain)
		coeff = as.vector(res1$glm.best.proteins$coef)
		pred1=as.vector(coeff%*%rbind(rep(1,length(vtrain)),mtrain[,as.numeric(res1$best.proteins.numbers)]))
		test.pred1=as.vector(coeff%*%rbind(rep(1,length(vtest)),mtest[,as.numeric(res1$best.proteins.numbers)]))
		Ftrain = Ftrain + c*pred1
		Ftest = Ftest + c*test.pred1
		iter=iter+1
		proteins[iter]=as.numeric(res1$best.proteins.numbers)
		coef[iter,]=as.vector(coeff)
		infl = c(infl,coef[iter,2]*sd(mtrain[,as.numeric(res1$best.proteins.numbers)]))
		mat.pred = cbind(mat.pred,Ftrain)
		mat.test.pred = cbind(mat.test.pred,Ftest)
		cat("\n Iteration: ",iter,"logLik: ",res1$glm.best.proteins$logLik)
	}
	mat.infl = cbind(1:stop,proteins,infl)
list(proteins=proteins,coef=coef,vtrain=vtrain,mat.pred=mat.pred,vtest=vtest,mat.test.pred=mat.test.pred,dev=dev,aic=aic,mat.infl=mat.infl)
}


#This function fits linear mixed model given the vector of response vtrain, covariate x, offset Ftrain, and id indicating matched observations
# or observations from the same cluster. The model can be fitted using ML or REML
estmixed <- function (vtrain, x, Ftrain, id, method="REML"){
  ## mi are individual estimates, sei are their standard errors
	Z = model.matrix(~factor(id)-1)
	X = cbind(rep(1,length=length(x)),x)
	n = nrow(Z)
	glmfit = glm(cbind(vtrain,1-vtrain)~x,family=binomial,offset=Ftrain)
	mu = fitted.values(glmfit)
	lin.pred = log(mu/(1-mu))
	W = diag(1/(mu*(1-mu)))

	
	v = (1 + exp(-lin.pred))*(1 + exp(lin.pred))
	tran.y = diag(v)%*%(vtrain - mu) + X%*%glmfit$coef
	
  	mll<-function (par, tran.y){
    		## calculate -2 * log likelihood
    		## par[1] is the grand mean`
    		## par[2] is the log of the between-group variance component
		V = exp(par)*Z%*%t(Z) + W
		L = chol(V,pivot=TRUE)
		y = backsolve(L,tran.y,transpose=TRUE)
		X = backsolve(L,X,transpose=TRUE)
		lm.reg = lm(y~X-1)
		b = list(coef=coef(lm.reg),pvalue=summary(lm.reg)$coef[2,4])

		#evaluate log likelihood
		#Vinv = solve(V)
		r = y - X%*%b$coef
		logLik = -n/2*log(2*pi) - sum(log(diag(L))) - sum((r)^2)/2
		attr(logLik,"fixed") <- b #allow retrieval of beta
		logLik
  	}

  	mll.reml<-function (par, tran.y){
		V = exp(par)*Z%*%t(Z) + W
		L = chol(V,pivot=TRUE)
		y = backsolve(L,tran.y,transpose=TRUE)
		X = backsolve(L,X,transpose=TRUE)
		lm.reg = lm(y~X-1)
		b = list(coef=coef(lm.reg),pvalue=summary(lm.reg)$coef[2,4])

		#evaluate log likelihood
		#Vinv = solve(V)
		r = y - X%*%b$coef
		p = length(coef)
		logLik = -((p-n)/2*log(2*pi) - sum(log(diag(L))) - sum((r)^2)/2 - as.numeric(determinant(t(X)%*%X)$mod)/2)
		#logLik = (p-n)/2*log(2*pi) - as.numeric(determinant(V)$mod)/2 - sum((r)^2)/2 - as.numeric(determinant(t(X)%*%X)$mod)/2
		attr(logLik,"fixed") <- b #allow retrieval of beta
		logLik
  	}
  	if(method!="ML" & method!="REML"){
   		stop("Invalid value of method")
  	}
    	res <- resid(lm(tran.y~x,weights=1/diag((W)))) 
	V.est = res%*%t(res)
	C = solve(t(Z)%*%Z)
	sigma = max(log(max(0,mean(diag(C%*%t(Z)%*%(V.est - W)%*%Z%*%C)))),3)
	objfun <- if(method == "ML") mll else mll.reml
	out <- nlminb(start=sigma,objfun,lower=c(-10),upper=c(3),tran.y=tran.y)
	if (out$convergence == 0) {
		tau.var = exp(out$par)
		b = attr(objfun(out$par,tran.y),"fixed")
		list(coef=c(beta.not=b$coef[1],beta.one=b$coef[2]),pvalue=b$pvalue,sigma=sqrt(tau.var),logLik=objfun(out$par,tran.y),out=out)
	}
	else {
		tau.var = 0
		list(coef=c(beta.not=0,beta.one=0),pvalue=1,sigma=sqrt(tau.var),logLik=objfun(out$par,tran.y),out=out)
	}
}

#This function calls the function estmixed to fit linear mixed model to all columns of the covariate matrix mtrain and response vector vtrain
#with an offset Ftrain. One best fitted covariate is selected for update 
fwd.lmer1=function(mtrain,vtrain,id,Ftrain){ #begin
	#Input parameters: mtrain for protein measurements in training
	#                  vtrain for binary response (0/1 group) in training 
	#                  w is weights

	do.glmm <- function(cont.var, bin.var,Ftrain)
	{
		dat = data.frame(bin.var=bin.var,cont.var=cont.var)
		res <- estmixed(bin.var, cont.var, Ftrain, id, method="REML")
		res <- res$pvalue
    	res
	}

	logreg.pvalue <- apply(mtrain, 2, FUN=do.glmm, vtrain, Ftrain)  #This instruction is a substitution of a "while" loop.

	# Get the best proteins in the fitted linear mixed models 

	ind = which.min(logreg.pvalue)
	best.proteins=mtrain[,ind[1]]
	best.proteins.numbers=ind[1]

	#run glm on best proteins set

	glm.best.proteins=estmixed(vtrain, best.proteins, Ftrain, id, method="REML")
	list(best.proteins.numbers=best.proteins.numbers,glm.best.proteins=glm.best.proteins)
}   #end


# Generating the Dataset

mean.response = function(x,u) {
	x = x[,1:10]
	eta = rep(1,dim(x)[1]) + x%*%rep(1,dim(x)[2]) + u
	p = (1 + exp(-eta))^(-1)
}

stop = 300
no.sim = 10
set.seed(12345)
x10 = c(runif(5,-2,-1),runif(5,1,2))}
#Figure 1 uses x10 = c(-1.279096,-1.124227,-1.239018,-1.113875,-1.543519,1.166372,1.325095,1.509224,1.727705,1.989737)
sig = c(.01,1,5,10)
op <- par(mfrow = c(2, 2))  # 2 x 2 pictures on one plot
time1 = NULL
time2 = NULL
time3 = NULL
Miss_Error = NULL
for (k in 1:length(sig)) {
	sigma.u = sig[k]
	vec.t1 = NULL
	vec.t2 = NULL
	vec.t3 = NULL
	mat.test.errorPQL_RE = NULL
	mat.test.errorPQL = NULL
	mat.test.errorWL2 = NULL
	mat.bayes = NULL
	for (m in 1:no.sim) {

#Generate training covariates
	ncov = 50
	rho = 0 #worked very well
	rhox = 1
	ntrain = 100
	x1 = matrix(nrow=ntrain,ncol=ncov)
	x2 = matrix(nrow=ntrain,ncol=ncov)
	x1 = mvrnorm(n = ntrain, mu=rep(0,ncov), Sigma = matrix(1,ncov,ncov)+diag(rep(1,ncov)), tol = 1e-6, empirical = FALSE)
	x1 = cbind(x1,mvrnorm(n=ntrain,mu=rep(0,ncov),Sigma = diag(rep(1,ncov)), tol = 1e-6, empirical = FALSE))
	covmat = matrix(rho,2*ncov,2*ncov)+diag(rep(.1,2*ncov)) # diag(rep(.1,ncov)) worked very well
	ep = mvrnorm(n = ntrain, mu=c(x10,rep(0,(2*ncov-10))), Sigma = covmat, tol = 1e-6, empirical = FALSE)
	x2 = rhox*x1 + ep
	mtrain = NULL
	for (i in 1:ntrain) {
		mtrain = rbind(mtrain,x1[i,],x2[i,])
	}

#Generate test covariates
	ncov = 50
	rho = 0 #worked very well
	rhox = 1
	ntest = 1000
	x1test = mvrnorm(n = ntest, mu=rep(0,ncov), Sigma = matrix(1,ncov,ncov)+diag(rep(1,ncov)), tol = 1e-6, empirical = FALSE)
	x1test = cbind(x1test,mvrnorm(n=ntest,mu=rep(0,ncov),Sigma = diag(rep(1,ncov)), tol = 1e-6, empirical = FALSE))
	ep = mvrnorm(n = ntest, mu=c(x10,rep(0,(2*ncov-10))), Sigma = covmat, tol = 1e-6, empirical = FALSE)
	x2test = rhox*x1test + ep
	mtest = NULL
	for (i in 1:ntest) {
		mtest = rbind(mtest,x1test[i,],x2test[i,])
	}

#Generate training response
	u = rnorm(ntrain,0,sigma.u)
	y1 = 2*rbinom(ntrain,1,mean.response(x1,u)) - 1
	y2 = 2*rbinom(ntrain,1,mean.response(x2,u)) - 1
	id = NULL 
	vtrain = NULL
	vtrainc = NULL
	mtrainc = NULL
	subid = NULL
	for (i in 1:ntrain) {
		id = c(id,rep(i,2))
		vtrain = c(vtrain,y1[i],y2[i])
	}
	train.dat = as.data.frame(cbind(id=id,vtrain=vtrain,mtrain=mtrain))
	train.dat.sub = train.dat[which(id %in% subid),]

#Generate test response
	u = rnorm(ntest,0,sigma.u)
	y1 = 2*rbinom(ntest,1,mean.response(x1test,u)) - 1
	y2 = 2*rbinom(ntest,1,mean.response(x2test,u)) - 1
	y1_bayes = ifelse(1 - pnorm(-(rep(1,dim(x1test[,1:10])[1]) + x1test[,1:10]%*%rep(1,dim(x1test[,1:10])[2]))/sigma.u)>.5,1,-1)
	y2_bayes = ifelse(1 - pnorm(-(rep(1,dim(x2test[,1:10])[1]) + x2test[,1:10]%*%rep(1,dim(x2test[,1:10])[2]))/sigma.u)>.5,1,-1)
	idtest = NULL 
	vtest = NULL
	vtestc = NULL
	mtestc = NULL
	subidtest = NULL
	vtest_bayes = NULL
	for (i in 1:ntest) {
		idtest = c(idtest,rep(i,2))
		vtest = c(vtest,y1[i],y2[i])
		vtest_bayes = c(vtest_bayes,y1_bayes[i],y2_bayes[i])
	}
	test.dat = as.data.frame(cbind(id=id,vtest=vtest,mtest=mtest))
	test.dat.sub = test.dat[which(idtest %in% subidtest),]
	bayes_rate = rep(length(which(vtest*vtest_bayes<0))/length(vtest),stop)

	c = .05
	ptm = proc.time()
	boost.object = boostingPQL(mtrain,vtrain,id,stop,mtest,vtest,idtest)
	t1 = proc.time() - ptm
	vec.t1 = c(vec.t1,t1[3])

	test.error = NULL
	for (i in 1:stop){
		test.pred = boost.object$mat.test.pred[,i]
		error.test = length(which(test.pred*boost.object$vtest<0))/length(test.pred)
		test.error = c(test.error,error.test)
	}
	mat.test.errorPQL = cbind(mat.test.errorPQL,test.error)

	
	ptm = proc.time()
	boost.object = boostingWL2(mtrain,vtrain,id,stop,mtest,vtest,idtest)
	t3 = proc.time() - ptm
	vec.t3 = c(vec.t3,t3[3])

	test.error = NULL
	for (i in 1:stop){
		test.pred = boost.object$mat.test.pred[,i]
		error.test = length(which(test.pred*boost.object$vtest<0))/length(test.pred)
		test.error = c(test.error,error.test)
	}
	mat.test.errorWL2 = cbind(mat.test.errorWL2,test.error)

	mat.bayes = cbind(mat.bayes,bayes_rate)
}
	errorPQL = apply(mat.test.errorPQL,1,mean)
	bayes_rate = apply(mat.bayes,1,mean)
	plot(1:stop,errorPQL,type="l",lty=1,ylim=c(0,.35),col = "red",xlab="Number of Iterations",ylab="Classification Error")
	text(floor(stop/8),.32, bquote(sigma == .(signif(sigma.u,1))))
	lines(1:stop,bayes_rate,type="l",lty=1,col="black")	
	errorWL2 = apply(mat.test.errorWL2,1,mean)
	lines(1:stop,errorWL2,type="l",lty=2,col = "blue")
	Miss_Error = cbind(Miss_Error,errorWL2,errorPQL)
	time1 = c(time1,mean(vec.t1))
	time3 = c(time3,mean(vec.t3))
}
old.sink.level <- sink("Miss_Error05.txt",append=FALSE) 
on.exit(sink(unsink.to=old.sink.level))
print(Miss_Error)
sink()
save.image()
