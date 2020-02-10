# ****************************************************
# Email: lsen@infervision.com
# Data: Wed Jan 10 20:15:45 2018
# Description: null
#                
# *****************************************************
rm(list = ls())
load("./_data/TCGA_BRCA_case_control_filter4.rdata")
load("./_data/data_clin.rdata")

data_case = data_case4
data_ctrl = data_ctrl4
rm(data_case4, data_ctrl4)

# save data
write.csv(data_case, file = "./data/input/input_data_case4.csv")
write.csv(data_ctrl, file = "./data/input/input_data_ctrl4.csv")
write.csv(data_clin, file = "./data/input/clinical_data.csv")

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

# save data 
write.csv(data_case10, file = "./data/input/input_data_case4_10.csv")
write.csv(data_ctrl10, file = "./data/input/input_data_ctrl4_10.csv")

write.csv(data_case50, file = "./data/input/input_data_case4_50.csv")
write.csv(data_ctrl50, file = "./data/input/input_data_ctrl4_50.csv")

write.csv(data_case100, file = "./data/input/input_data_case4_100.csv")
write.csv(data_ctrl100, file = "./data/input/input_data_ctrl4_100.csv")

write.csv(data_case1000, file = "./data/input/input_data_case4_1000.csv")
write.csv(data_ctrl1000, file = "./data/input/input_data_ctrl4_1000.csv")

write.csv(data_case1500, file = "./data/input/input_data_case4_1500.csv")
write.csv(data_ctrl1500, file = "./data/input/input_data_ctrl4_1500.csv")

# save data  normlization to 0~1
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

write.csv(data_case, file = "./data/input/input_data_case4_norm.csv")
write.csv(data_ctrl, file = "./data/input/input_data_ctrl4_norm.csv")

write.csv(data_case10, file = "./data/input/input_data_case4_norm_10.csv")
write.csv(data_ctrl10, file = "./data/input/input_data_ctrl4_norm_10.csv")

write.csv(data_case50, file = "./data/input/input_data_case4_norm_50.csv")
write.csv(data_ctrl50, file = "./data/input/input_data_ctrl4_norm_50.csv")

write.csv(data_case100, file = "./data/input/input_data_case4_norm_100.csv")
write.csv(data_ctrl100, file = "./data/input/input_data_ctrl4_norm_100.csv")

write.csv(data_case1000, file = "./data/input/input_data_case4_norm_1000.csv")
write.csv(data_ctrl1000, file = "./data/input/input_data_ctrl4_norm_1000.csv")

write.csv(data_case1500, file = "./data/input/input_data_case4_norm_1500.csv")
write.csv(data_ctrl1500, file = "./data/input/input_data_ctrl4_norm_1500.csv")



### （1）ttest系列：Paired t-test
source("./S1-ttest.R")

time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
rm(i, top, datacase, datactrl,result_pttest)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  ptm = proc.time()
  result_pttest = test_pttest(datacase, datactrl)
  time_pttest[i,] = c(top, proc.time() - ptm)
  write.csv(result_pttest, file = paste0("./data/output/1_pttest_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/1_pttest_time.csv", row.names = F)


### （2）ttest系列：Modified Paired t-test （Tan et al.）
source("./S1-ttest.R")

time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
rm(i, top, datacase, datactrl,result_pttest)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  ptm = proc.time()
  result_pttest = test_mpttest(datacase, datactrl)
  time_pttest[i,] = c(top, proc.time() - ptm)
  write.csv(result_pttest, file = paste0("./data/output/2_mpttest_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/2_mpttest_time.csv", row.names = F)



### （3）ttest系列：Fold-Change Paired t-test (Cao et al.)
source("./S1-ttest.R")

time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
rm(i, top, datacase, datactrl,result_pttest)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  ptm = proc.time()
  result_pttest = test_fcpttest(datacase, datactrl)
  time_pttest[i,] = c(top, proc.time() - ptm)
  write.csv(result_pttest, file = paste0("./data/output/3_fcpttest_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/3_fcpttest_time.csv", row.names = F)


### （4）CLR系列：RPCLR, Random penalized conditional logistic regression
source("./S3-CLR.R")
#load("./data/result_rpclr.rdata")

mtry = 3
#numBS = 4000 #ncol(data_case) * 2

time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
v_numBS = c(100, 500, 1000, 2000, 3000)
rm(i, top, datacase, datactrl,result_pttest)
for(i in 1:5){
  #  i = 5
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  ptm = proc.time()
  result_pttest = clr_rpclr(datacase, datactrl, mtry, v_numBS[i])
  time_pttest[i,] = c(top, proc.time() - ptm)
  write.csv(result_pttest, file = paste0("./data/output/4_rpclr_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/4_rpclr_time.csv", row.names = F)


### （5）CLR系列： （待续）
source("./S3-CLR.R")

time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
# v_numBS = c(3000, 3000, 3000, 3000, 3000)
v_numBS = c(100, 500, 1000, 2000, 3000)
rm(i, top, datacase, datactrl)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  datactrl = get(paste0("data_ctrl",top))
  
  ptm = proc.time()
  clogitPenalty(datacase, datactrl, numLambda = v_numBS[i], alpha = 0.5, save_prefix = top)
  time_pttest[i,] = c(top, proc.time() - ptm)
  #write.csv(result_pttest, file = paste0("./data/output/5_clogitL1_",top,".csv"), row.names = F)
}
time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/5_clogitL1_time.csv", row.names = F)


### （6）CLR系列：Baysian CLR
# 见 Python notebook

### （7）Boosting系列： boosting-PQL
# 见 S22-boosting.R

### （8）Boosting系列： boosting-WL2
# 见 S22-boosing.R
