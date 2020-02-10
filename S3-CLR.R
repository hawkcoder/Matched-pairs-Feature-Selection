# ****************************************************
# Email: lsen@infervision.com
# Data: Wed Dec 20 09:59:38 2017
# Description: Conditional Logistic Regression
#                
# *****************************************************

# load("./_data/TCGA_BRCA_case_control_filter.rdata")
library(RPCLR)
library(clogitL1)
library(Hmisc)
library(stringr)

# *************************
# CLR: (1) RPCLR, Random penalized conditional logistic regression
# *************************
library(RPCLR)

# mtry = 3
# numBS = 4000 #ncol(data_case2) * 2
# 
# result_rpclr = clr_rpclr(data_case2, data_ctrl2, mtry, numBS)
# 

clr_rpclr = function(data_case, data_ctrl, mtry=10, numBS = 2000)
{
  
  if(any(dim(data_case) != dim(data_ctrl))) {
    stop("Error: dim of data_case and data_ctrl are not equal!")
  }
  
  data_case = t(data_case)
  data_ctrl = t(data_ctrl)
  

  data = rbind(data_ctrl, data_case)
  rownames(data) = seq(1, nrow(data))
  out = c(rep(0, nrow(data_ctrl)), rep(1, nrow(data_case)))
  strat = c(seq(1,nrow(data_ctrl)), seq(1, nrow(data_case)))
  
  # make normalization
  for(i in 1:ncol(data)){
    data[, i] = scale(data[,i], center = T, scale = T)
  }
  
  # Get Variable Importance
  MyResults = GetVarImp(data, out, strat, mtry = mtry, numBS = numBS)
  
  
  df_result <- data.frame(features= colnames(data), importance=MyResults)
  df_result = df_result[order(-df_result$importance),]
  return(df_result)
  
}



# *************************
# CLR: (2) Penalized conditional and unconditonal logistic regression
# *************************
library(clogitL1)

# flattenCorrMatrix
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  v_row = rownames(cormat)[row(cormat)[ut]]
  v_col = rownames(cormat)[col(cormat)[ut]]
  v_new_features = paste(v_row, v_col, sep = "_")
  
  return(data.frame( row = v_row, column = v_col, new_feature = v_new_features, cor =(cormat)[ut], pvalue = pmat[ut] ))
}


## Penalized conditional and unconditonal logistic regression
clogitPenalty =  function(datacase, datactrl, numLambda = 1000, alpha = 0.5, save_prefix = "feature_top_10"){
  
  data_case2 = as.data.frame(t(datacase))
  data_ctrl2 = as.data.frame(t(datactrl))
  
  data_head = data.frame(strata = 1:nrow(data_case2), y = 1)
  data_case2 = cbind(data_head, data_case2)
  # data_case2[1:10, 1:5]
  
  data_head2 = data.frame(strata = 1:nrow(data_ctrl2), y = 0)
  data_ctrl2 = cbind(data_head2, data_ctrl2)
  # data_ctrl2[1:10, 1:5]
  
  data = rbind(data_case2, data_ctrl2)
  rownames(data) = 1:nrow(data)
  
  # numLambda = 1000
  # alpha = 0.5
  
  data_case = data[which(data$y == 1), -c(1,2)]
  data_ctrl = data[which(data$y == 0), -c(1,2)]
  
  v_features = colnames(data)[-c(1,2)]
  
  
  # Step1: 只用Main Feature去做penalty conditional Regression
  mod1 = clogitL1(y=data$y, 
                   strata = data$strata,
                   x=as.matrix(data[,v_features]),
                   numLambda = numLambda, 
                   minLambdaRatio = 0.000001, 
                   alpha = alpha)
  cv.mod1 = cv.clogitL1(mod1)

  # par(mfrow = c(1,2))
  # plot(mod10, logX = T)
  # plot(cv.mod10)
  
  # 选择最优解即lambda 及对应的beta
  mod1.smy = summary(mod1)
  cv.mod1.smy = summary(cv.mod1)
  n_best = which.min(abs(cv.mod1.smy$lambda_minCV - cv.mod1$lambda))
  
  coef_best = cv.mod1$beta[n_best,]
  names(coef_best) = colnames(data)[-c(1,2)]
  v_feature_selected = coef_best[coef_best != 0]        #这就是 EN1步骤选择出来的Xmain特征
  
  
  # df_result_pclrn10 = data.frame(features = names(coef_best), coef=coef_best)
  # df_result_pclrn10 = df_result_pclrn10[order(-abs(df_result_pclrn10$coef)), ]
  # df_result_pclrn10
  #### But we need add some interaction of feature in int 
  
  # STEP2：去掉筛选掉的feature,并与two-way interaction feature合并
  data2 = data[,c("strata", "y", names(v_feature_selected))]

  # 做Correlation，选择p < 0.0005
  mod1_cor_case = rcorr(as.matrix(data_case[,names(v_feature_selected)]))
  mod1_cor_case = flattenCorrMatrix(mod1_cor_case$r, mod1_cor_case$P)
  mod1_cor_case = mod1_cor_case[which(mod1_cor_case$pvalue < 0.005),]
  # mod1_cor_case
  
  mod1_cor_ctrl = rcorr(as.matrix(data_ctrl[,names(v_feature_selected)]))
  mod1_cor_ctrl = flattenCorrMatrix(mod1_cor_ctrl$r, mod1_cor_ctrl$P)
  mod1_cor_ctrl = mod1_cor_ctrl[which(mod1_cor_ctrl$pvalue < 0.005),]
  #mod1_cor_ctrl
  
  ## 选择两者都有的interaction features， 这是选择出有交叉的特征
  v_mod1_inter = intersect(mod1_cor_case$new_feature, mod1_cor_ctrl$new_feature)
  v_mod1_inter
  
  ## 创建cor的矩阵
  n_trol = nrow(data_ctrl)
  n_case = nrow(data_case)
  data_cor_feature = data.frame(matrix(NA, nrow = n_trol+n_case, ncol = length(v_mod1_inter)))
  colnames(data_cor_feature) = v_mod1_inter
  
  for (newfeature in v_mod1_inter) {
    # newfeature = v_mod1_inter[1]
    cor_case = mod1_cor_case$cor[which(mod1_cor_case$new_feature == newfeature)]
    cor_ctrl = mod1_cor_ctrl$cor[which(mod1_cor_ctrl$new_feature == newfeature)]
    data_cor_feature[newfeature] = c(rep(cor_case, n_case), rep(cor_ctrl, n_trol))
  }
  
  # 与原始数据合并在一起
  data3 = cbind(data2, data_cor_feature)
  
  #### STEP3：开始训练
  mod2 = clogitL1(y=data3$y, 
                     strata = data3$strata,
                     x=as.matrix(data3[,-c(1,2)]),
                     numLambda = numLambda, 
                     minLambdaRatio = 0.000001, 
                     alpha = alpha)
  cv.mod2 = cv.clogitL1(mod2)
  
  # par(mfrow = c(1,2))
  # plot(mod2, logX = T)
  # plot(cv.mod2)
  
  mod1.smy2 = summary(mod2)
  cv.mod1.smy2 = summary(cv.mod2)
  
  n_best2 = which.min(abs(cv.mod1.smy2$lambda_minCV - cv.mod2$lambda))
  coef_best2 = cv.mod2$beta[n_best2,]
  names(coef_best2) = colnames(data3)[-c(1,2)]
  v_feature_selected2 = coef_best2[coef_best2 != 0] ## 选择从不非0的特征
  
  ## 建立结果矩阵
  df_result_pclrn10 = data.frame(features = names(v_feature_selected2), coef=v_feature_selected2, stringsAsFactors = F)
  df_result_pclrn10 = df_result_pclrn10[order(-abs(df_result_pclrn10$coef)), ]
  #df_result_pclrn10
  
  
  
  v_order_feature_list = unlist(strsplit(df_result_pclrn10$features, split = "_"))
  v_order_feature_list = v_order_feature_list[!duplicated(v_order_feature_list)]
  #v_order_feature_list
  
  #将之前的数据合并在一起
  v_not_in_order_feature_list = v_features[ !(v_features %in% v_order_feature_list) ]
  
  df_result_pclrn10_all = data.frame(features = c(v_order_feature_list, v_not_in_order_feature_list), 
                                     value = c(rep(1, length(v_order_feature_list)), rep(0, length(v_not_in_order_feature_list))), 
                                     stringsAsFactors = F)
  
  ## 保存
  #df_result_pclrn10
  #df_result_pclrn10_all
  # save(df_result_pclrn10, df_result_pclrn10_all, file = paste0("./data/output/5_clogitL1_",save_prefix, ".rdata"))
  # write.csv(df_result_pclrn10, file = paste0("./_data/output2/5_clogitL1_",save_prefix, ".ceof.csv"), row.names = F)
  write.csv(df_result_pclrn10_all, file = paste0("./data/output/5_clogitL1_",save_prefix, ".csv"), row.names = F)
  
}
