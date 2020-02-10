# ****************************************************
# Email: lsen@infervision.com
# Data: Wed Jan 10 19:59:10 2018
# Description: null
#                
# *****************************************************
rm(list = ls())
load("./_data/TCGA_BRCA_case_control_filter.rdata")


set.seed(7)
v_sample_id = sample(nrow(data_case4))

data_case4 = data_case2[v_sample_id, ]
data_ctrl4 = data_ctrl2[v_sample_id, ]

# （1）保存Filter后的数据
save(data_case4, data_ctrl4, file = "./_data/TCGA_BRCA_case_control_filter4.rdata")
# save to csv data
write.csv(data_case4, file = "./_data/TCGA_BRCA_case_filter4.csv")
write.csv(data_ctrl4, file = "./_data/TCGA_BRCA_ctrl_filter4.csv")


###########GEO data shuffte
rm(list = ls())
load("./_data/input2/GSE70947_case_control_filter.rdata")


set.seed(7)
v_sample_id = sample(nrow(data_case2))

data_case4 = data_case2[v_sample_id, ]
data_ctrl4 = data_ctrl2[v_sample_id, ]

# （1）保存Filter后的数据
save(data_case4, data_ctrl4, data_clinical, file = "./_data/
     ")

# save to csv data
write.csv(data_case4, file = "./_data/GSE70947_case_filter4.csv")
write.csv(data_ctrl4, file = "./_data/GSE70947_ctrl_filter4.csv")
