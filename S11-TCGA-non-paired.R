# ****************************************************
# Email: lsen@infervision.com
# Data: Mon Jan 22 22:57:18 2018
# Description: null
#                
# *****************************************************
rm(list = ls())

library(mRMRe)
library(RWeka)

# set.thread.count(10)

rm(list = ls())
load("./_data/TCGA_BRCA_case_control_filter4.rdata")
load("./_data/data_clin.rdata")

data_case = data_case4
data_ctrl = data_ctrl4
rm(data_case4, data_ctrl4)

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

write.csv(data_case, file = "./_data/input/input_data_case4_norm.csv")
write.csv(data_ctrl, file = "./_data/input/input_data_ctrl4_norm.csv")

write.csv(data_case10, file = "./_data/input/input_data_case4_norm_10.csv")
write.csv(data_ctrl10, file = "./_data/input/input_data_ctrl4_norm_10.csv")

write.csv(data_case50, file = "./_data/input/input_data_case4_norm_50.csv")
write.csv(data_ctrl50, file = "./_data/input/input_data_ctrl4_norm_50.csv")

write.csv(data_case100, file = "./_data/input/input_data_case4_norm_100.csv")
write.csv(data_ctrl100, file = "./_data/input/input_data_ctrl4_norm_100.csv")

write.csv(data_case1000, file = "./_data/input/input_data_case4_norm_1000.csv")
write.csv(data_ctrl1000, file = "./_data/input/input_data_ctrl4_norm_1000.csv")

write.csv(data_case1500, file = "./_data/input/input_data_case4_norm_1500.csv")
write.csv(data_ctrl1500, file = "./_data/input/input_data_ctrl4_norm_1500.csv")


# *************************
# 10- MRMD 
# *************************

# 写arff文件
data10 = as.data.frame(t(cbind(data_case10, data_ctrl10)))
data10$class = c(rep(1, ncol(data_case10)), rep(0, ncol(data_ctrl10)))
write.arff(data10, file = "./_data/input/arff/MRMD-master/input_data_case4_norm_10.arff")

data50 = as.data.frame(t(cbind(data_case50, data_ctrl50)))
data50$class = c(rep(1, ncol(data_case50)), rep(0, ncol(data_ctrl50)))
write.arff(data50, file = "./_data/input/arff/MRMD-master/input_data_case4_norm_50.arff")

data100 = as.data.frame(t(cbind(data_case100, data_ctrl100)))
data100$class = c(rep(1, ncol(data_case100)), rep(0, ncol(data_ctrl100)))
write.arff(data100, file = "./_data/input/arff/MRMD-master/input_data_case4_norm_100.arff")

data1000 = as.data.frame(t(cbind(data_case1000, data_ctrl1000)))
data1000$class = c(rep(1, ncol(data_case1000)), rep(0, ncol(data_ctrl1000)))
write.arff(data1000, file = "./_data/input/arff/MRMD-master/input_data_case4_norm_1000.arff")

data1500 = as.data.frame(t(cbind(data_case1500, data_ctrl1500)))
data1500$class = c(rep(1, ncol(data_case1500)), rep(0, ncol(data_ctrl1500)))
write.arff(data1500, file = "./_data/input/arff/MRMD-master/input_data_case4_norm_1500.arff")

# 读取处理的文件
v_top = c(10, 50, 100, 1000, 1500)
rm(i, top, df_result, datacase)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  
  df_result = read.table(file = paste0("/home/lsen/Data/CODE/feature_selection/_data/input/arff/MRMD-master/10_nonparied_MRMD_",top,".txt")
                         ,sep = "\t", header = T, stringsAsFactors = F)
  v_fea_index = sapply(df_result$FeaName, function(ss){as.numeric(substr(ss, 4, nchar(ss)))}, simplify = T, USE.NAMES = F)
  v_fea_index = v_fea_index + 1
  
  df_result$FeaName = rownames(datacase)[v_fea_index]
  colnames(df_result) = c("rank", "features", "score")
  write.csv(df_result, file = paste0("./data/output/10_nonparied_MRMD_",top,".csv"), row.names = F)

}



# *************************
# 9- mRMR algorithm
# *************************
time_pttest = data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(time_pttest) = c("Top_n_feature","user.self","sys.self", "elapsed", "user.child", "sys.child")

v_top = c(10, 50, 100, 1000, 1500)
rm(i, top, datacase, datactrl,result_pttest)
for(i in 1:5){
  #  i = 1
  top = v_top[i]
  datacase = get(paste0("data_case",top))
  
  data = mRMR.data(data.frame(t(datacase)))
  
  
  ptm = proc.time()
  ss = mRMR.ensemble(solution_count = 1, data=data, target_indices=1, feature_count = top)
  time_pttest[i,] = c(top, proc.time() - ptm)
  
  
  v_filters = as.vector(ss@filters$`1`)
  v_scores = as.vector(ss@scores$`1`)
  v_features = ss@feature_names[v_filters]
  
  
  ## 建立结果矩阵
  df_result = data.frame(features = v_features, score=ss@scores$`1`, rank = seq(top, 1, -1), stringsAsFactors = F)
  df_result = df_result[order(df_result$rank), ]
  
  df_result$features = stringr::str_replace(df_result$features, "\\.", "|")
  
  write.csv(df_result, file = paste0("./data/output/9_nonparied_mRMR_",top,".csv"), row.names = F)
}

time_pttest = round(time_pttest, 4)
write.csv(time_pttest, file = "./data/output/9_nonparied_mRMR_time.csv", row.names = F)



