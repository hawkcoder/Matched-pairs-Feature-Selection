library(MASS)

# *************************
# Test for boosting
# *************************

# TCGA-BRCA data
rm(list = ls())
load("./_data/TCGA_BRCA_case_control_filter4.rdata")
load("./_data/data_clin.rdata")

## GEO data
# rm(list = ls())
# load("./_data/GSE70947_case_control_filter4.rdata")
# 

data_case = data_case4
data_ctrl = data_ctrl4
rm(data_case4, data_ctrl4)

data_case = (data_case - min(data_case)) / (max(data_case) - min(data_case))
data_ctrl = (data_ctrl - min(data_ctrl)) / (max(data_ctrl) - min(data_ctrl))


data_case10 = data_case[c(1:10),]
data_case50 = data_case[c(1:50),]
data_case100 = data_case[c(1:100),]
data_case1000 = data_case[c(1:1000),]
data_case1500 = data_case[c(1:1500),]

data_ctrl10 = data_ctrl[c(1:10),]
data_ctrl50 = data_ctrl[c(1:50),]
data_ctrl100 = data_ctrl[c(1:100),]
data_ctrl1000 = data_ctrl[c(1:1000),]
data_ctrl1500 = data_ctrl[c(1:1500),]

# 
# table(data_clin$vital_status)
# 
# data_train_case = t(data_case10)
# data_train_ctrl = t(data_ctrl10)
# data_train = rbind(data_train_case, data_train_ctrl)
# rownames(data_train) = 1:nrow(data_train)
# 
# dim(data_train_case)
# dim(data_train_ctrl)
# dim(data_train)
# 
# data_train[, 1:5]
# 
# v_index = c(1:nrow(data_train_case),1:nrow(data_train_ctrl))
# v_target = c(rep(1, nrow(data_train_case)), rep(0, nrow(data_train_ctrl)) )


# v_target = data_clin$vital_status
# v_target = c(v_target, v_target)
# v_target = ifelse( v_target == "alive", 1, -1)

# v_target_stage = data_clin$tumor_stage
# v_target_stage = c(v_target_stage, v_target_stage)
# v_target_stage[v_target_stage %in% c("stage I", "stage II")] = 0
# v_target_stage[v_target_stage %in% c("stage III", "stage IV")] = 1



# # 使用预测死亡数据来
# c = .05
# ptm = proc.time()
# boost.PQL = boostingPQL(data_train, v_target, v_index, 2500)
# t1 = proc.time() - ptm
# save(boost.PQL, file = "./_data/boosting.PQL.3000.rdata")

# *************************
# Boost PQL
# *************************
time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
# v_numBS = c(100, 100, 200, 200, 200)
v_numBS = c(100, 200, 500, 2000, 3000)
rm(i, top, datacase, datactrl,data_train_case, data_train_ctrl, boost.WL2, df_result)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  niter = v_numBS[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  # tide data for run
  data_train_case = t(datacase)
  data_train_ctrl = t(datactrl)
  data_train = rbind(data_train_case, data_train_ctrl)
  rownames(data_train) = 1:nrow(data_train)
  
  # get index and predict target 
  v_index = c(1:nrow(data_train_case),1:nrow(data_train_ctrl))
  v_target = c(rep(1, nrow(data_train_case)), rep(0, nrow(data_train_ctrl)) )
  
  # run 
  c = .05
  ptm = proc.time()
  boost.WL2 = boostingPQL(data_train, v_target, v_index, niter)
  #boost.WL2 = boostingWL2(data_train, v_target, v_index, niter)
  time_pttest[i,] = c(top, proc.time() - ptm)
  # save(boost.WL2, file = paste0("./_data/output2/boosting.PQL.top",top,".",niter,".rdata"))
  
  # get variable importance
  v_protein = colnames(data_train)
  boost.WL2.infl = as.data.frame(boost.WL2$mat.infl)
  boost.WL2.infl$features = v_protein[boost.WL2.infl$proteins]
  
  # 排序
  boost.WL2.infl = boost.WL2.infl[order(-boost.WL2.infl$infl),]
  
  # 去掉重复的 features
  v_duplicated = duplicated(boost.WL2.infl$features)
  boost.WL2.infl2 = boost.WL2.infl[!duplicated(boost.WL2.infl$features),]
  
  # 合并成一个
  df_result = data.frame(features= boost.WL2.infl2$features, importance = boost.WL2.infl2$infl, stringsAsFactors = F)
  df_resul2 = data.frame(features= v_protein[!(v_protein %in% boost.WL2.infl2$features)], importance = 0, stringsAsFactors = F)
  df_result = rbind(df_result, df_resul2)
  
  # save 
  write.csv(df_result, file = paste0("./data/output/7_boostPQL_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/7_boostPQL_time.csv", row.names = F)



# *************************
# Boosting WL2
# *************************

time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
v_numBS = c(100, 200, 500, 2000, 3000)
rm(i, top, datacase, datactrl,data_train_case, data_train_ctrl, boost.WL2, df_result)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  niter = v_numBS[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  # tide data for run
  data_train_case = t(datacase)
  data_train_ctrl = t(datactrl)
  data_train = rbind(data_train_case, data_train_ctrl)
  rownames(data_train) = 1:nrow(data_train)
  
  # get index and predict target 
  v_index = c(1:nrow(data_train_case),1:nrow(data_train_ctrl))
  v_target = c(rep(1, nrow(data_train_case)), rep(0, nrow(data_train_ctrl)) )
  
  # run 
  c = .05
  ptm = proc.time()
  boost.WL2 = boostingWL2(data_train, v_target, v_index, niter)
  time_pttest[i,] = c(top, proc.time() - ptm)
  # save(boost.WL2, file = paste0("./data/output/boosting.WL2.top",top,".",niter,".rdata"))
  
  # get variable importance
  v_protein = colnames(data_train)
  boost.WL2.infl = as.data.frame(boost.WL2$mat.infl)
  boost.WL2.infl$features = v_protein[boost.WL2.infl$proteins]
  
  # 排序
  boost.WL2.infl = boost.WL2.infl[order(-boost.WL2.infl$infl),]
  
  # 去掉重复的 features
  v_duplicated = duplicated(boost.WL2.infl$features)
  boost.WL2.infl2 = boost.WL2.infl[!duplicated(boost.WL2.infl$features),]
  
  # 合并成一个
  df_result = data.frame(features= boost.WL2.infl2$features, importance = boost.WL2.infl2$infl, stringsAsFactors = F)
  df_resul2 = data.frame(features= v_protein[!(v_protein %in% boost.WL2.infl2$features)], importance = 0, stringsAsFactors = F)
  df_result = rbind(df_result, df_resul2)
  
  # save 
  write.csv(df_result, file = paste0("./data/output/8_boostWL2_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/8_boostWL2_time.csv", row.names = F)





# load("./_data/boosting.WL2.3000.rdata")
# View(boost.WL2$mat.infl)
# 
# v_protein = colnames(data_train)
# boost.WL2.infl = as.data.frame(boost.WL2$mat.infl)
# boost.WL2.infl$features = v_protein[boost.WL2.infl$proteins]
# 
# # 排序
# boost.WL2.infl = boost.WL2.infl[order(-boost.WL2.infl$infl),]
# 
# # 去掉重复的 features
# v_duplicated = duplicated(boost.WL2.infl$features)
# boost.WL2.infl2 = boost.WL2.infl[!duplicated(boost.WL2.infl$features),]
# 
# # 合并成一个
# df_result = data.frame(features= boost.WL2.infl2$features, importance = boost.WL2.infl2$infl, stringsAsFactors = F)
# df_resul2 = data.frame(features= v_protein[!(v_protein %in% boost.WL2.infl2$features)], importance = 0, stringsAsFactors = F)
# df_result = rbind(df_result, df_resul2)
# 
# sum(df_result$features %in% colnames(data_train))
# 
# # 保存
# result_boost_wl2 = df_result
# save(result_boost_wl2, file = "./_data/result_boost_wl2.rdata")
# write.csv(df_result, file = "./_data/result_boost_wl2.csv", row.names = F)


# *************************
# Functions
# *************************
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
    infl = c(infl, coef[iter]*sd(f))
    mat.pred = cbind(mat.pred,Ftrain)
    mat.test.pred = cbind(mat.test.pred,Ftest)
    cat("\n Iteration: ",iter)
  }
  mat.infl = cbind(1:stop,proteins,infl)
  list(proteins=proteins,coef=coef,vtrain=vtrain,
       mat.pred=mat.pred,vtest=vtest,
       mat.test.pred=mat.test.pred,dev=dev,aic=aic,mat.infl=mat.infl)
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
    mat.pred = cbind(mat.pred, Ftrain)
    mat.test.pred = cbind(mat.test.pred,Ftest)
    cat("\n Iteration: ",iter,"logLik: ",res1$glm.best.proteins$logLik)
  }
  mat.infl = cbind(1:stop,proteins,infl)
  
  list(proteins=proteins,coef=coef,
       vtrain=vtrain,mat.pred=mat.pred,
       vtest=vtest,mat.test.pred=mat.test.pred,
       dev=dev,aic=aic,mat.infl=mat.infl)
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

