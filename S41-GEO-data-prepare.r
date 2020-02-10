# ****************************************************
# Email: lsen@infervision.com
# Data: Wed Jan 17 21:53:45 2018
# Description: null
#                
# *****************************************************

library(stringr)

load("./_data/data_case_ctrol_GSE70947.rdata")


# Step 1: 取log
# GEO 数据已经取了log
data_case2 = data_case
data_ctrl2 = data_ctrl

# Step 2: 删掉mean值都小于1的feature
mean_case = rowMeans(data_case2)
mean_ctrl = rowMeans(data_ctrl2)

i_mean_case = which(mean_case < 1) 
i_mean_ctrl = which(mean_ctrl < 1) 
# 没有mean特别小的数据

max(data_case)
min(data_case)

max(data_ctrl)
min(data_ctrl)

# Step 3: 删除var特别小的数据
v_var_case = sapply(1:nrow(data_case), function(i){var(as.numeric(data_case[i,]))})
v_var_ctrl = sapply(1:nrow(data_ctrl), function(i){var(as.numeric(data_ctrl[i,]))})

hist(v_var_case)
hist(v_var_ctrl)

sum(v_var_case > 0.8)
sum(v_var_ctrl > 0.8)

v_gene_var_case = which(v_var_case >= 0.8)
v_gene_var_ctrl = which(v_var_ctrl >= 0.8)

v_gene_var = unique(v_gene_var_case, v_gene_var_ctrl)

data_case2 = data_case[v_gene_var, ]
data_ctrl2 = data_ctrl[v_gene_var, ]


# Step 4: 选择一些FC > 0.5的
mean_case = rowMeans(data_case2)
mean_ctrl = rowMeans(data_ctrl2)
mean_FC = mean_case - mean_ctrl

hist(mean_FC)

i_mean_FC = which(abs(mean_FC) > 0.8) 

data_case2 = data_case2[i_mean_FC, ]
data_ctrl2 = data_ctrl2[i_mean_FC, ]


# Step 4: 做T-test检验
v_pvalue = vector(mode = 'numeric', length = nrow(data_case2))
for (i in 1:nrow(data_case2)) {
  v_pvalue[i] = t.test(data_case2[i,],  data_ctrl2[i,])$p.value
}

#hist(-log10(v_pvalue), breaks = 100)

# 去掉p-value > 0.005的
sum(v_pvalue < 0.005)

data_case2 = data_case2[which(v_pvalue < 0.005), ]
data_ctrl2 = data_ctrl2[which(v_pvalue < 0.005), ]


# Step 5: 验证数据
data_case2 = as.data.frame(data_case2)
data_ctrl2 = as.data.frame(data_ctrl2)

# (1) 数据维度是否一致
dim(data_case2)
dim(data_ctrl2)

# (2) sample ID 是否重复
length(unique(colnames(data_case2))) == ncol(data_case2)
length(unique(colnames(data_ctrl2))) == ncol(data_ctrl2)

# (3) case/contrl的sample ID是否一致
sum(colnames(data_case2) == colnames(data_ctrl2)) == ncol(data_case2)

# 肯定不一致，the GSE name of tumor or normal sample are not same in GEO database 

# (4) gene ID 是否重复
length(unique(rownames(data_case2))) == nrow(data_case2)
length(unique(rownames(data_ctrl2))) == nrow(data_ctrl2)

# (5) case/contrl的gene ID是否一致
sum(rownames(data_case2) == rownames(data_ctrl2)) == nrow(data_case2)

# save data
save(data_case2, data_ctrl2, data_clinical, file = "./_data/input2/GSE70947_case_control_filter.rdata")

# save to csv data
write.csv(data_case2, file = "./_data/input2/GSE70947_case_filter.csv")
write.csv(data_ctrl2, file = "./_data/input2/GSE70947_ctrl_filter.csv")


##### shuffet 一下gene的顺序
set.seed(7)
v_sample_id = sample(nrow(data_case2))

data_case4 = data_case2[v_sample_id, ]
data_ctrl4 = data_ctrl2[v_sample_id, ]

# 保存Filter后的数据
save(data_case4, data_ctrl4, data_clinical, file = "./_data/GSE70947_case_control_filter4.rdata")

# save to csv data
write.csv(data_case4, file = "./_data/GSE70947_case_filter4.csv")
write.csv(data_ctrl4, file = "./_data/GSE70947_ctrl_filter4.csv")
