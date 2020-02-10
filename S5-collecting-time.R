# ****************************************************
# Email: lsen@infervision.com
# Data: Mon Jan 29 21:58:58 2018
# Description: null
#                
# *****************************************************

df_time_TCGA = data.frame(matrix(NA, nrow = 5, ncol = 10), stringsAsFactors = F)

v_method_name = c("pttest", "mpttest", "fcpttest", "RP-CLR", "PCU-CLR", "BVS-CLR", "1-Step PQLBoost","WL2Boost", "mRMR", "MRMD")

colnames(df_time_TCGA) = v_method_name
rownames(df_time_TCGA) = c("gene10", "gene50", "gene100", "gene1000", "gene1500")


df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/1_pttest_time.csv")
df_time_TCGA$pttest = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/2_mpttest_time.csv")
df_time_TCGA$mpttest = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/3_fcpttest_time.csv")
df_time_TCGA$fcpttest = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/4_rpclr_time.csv")
df_time_TCGA$`RP-CLR` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/5_clogitL1_time.csv")
df_time_TCGA$`PCU-CLR` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/6_bayesian_time.csv")
df_time_TCGA$`BVS-CLR` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/7_boostPQL_time.csv")
df_time_TCGA$`1-Step PQLBoost` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/8_boostWL2_time.csv")
df_time_TCGA$WL2Boost = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/9_nonparied_mRMR_time.csv")
df_time_TCGA$mRMR = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/TCGA-BRCA/output/10_nonparied_MRMD_time.csv")
df_time_TCGA$MRMD = df_tmp$elapsed


mycolor = rainbow(10, alpha = 1)

df_time_TCGA[is.na(df_time_TCGA)] = Inf


par(xpd=F, mar=c(5,4,4,13))
par(xpd=T, mar=c(5,4,4,13))

plot(x= NA,xlim = c(1,5),ylim = c(0,5000), xlab = "Five gene groups", xaxt='n',
     ylab="Runing times (s)", main = "TCGA dataset running time comparison")
axis(1, at=c(1:5), label = c("10-genes", "50-genes", "100-genes", "1000-genes", "1500-genes"))

for (i in 1:10) {
  # i = 7
  points(x = c(1:5), y=df_time_TCGA[,i], pch=i, col=mycolor[i], lwd=2)
  lines(x = c(1:5), y=df_time_TCGA[,i], pch=i, col=mycolor[i], lwd=2)
}

legend(5.2, 5000, legend = v_method_name, pch = c(1:10), col=mycolor, lwd=2)


##### GEO 
df_time_GEO = data.frame(matrix(NA, nrow = 5, ncol = 10), stringsAsFactors = F)

v_method_name = c("pttest", "mpttest", "fcpttest", "RP-CLR", "PCU-CLR", "BVS-CLR", "1-Step PQLBoost","WL2Boost", "mRMR", "MRMD")

colnames(df_time_GEO) = v_method_name
rownames(df_time_GEO) = c("gene10", "gene50", "gene100", "gene1000", "gene1500")


df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/1_pttest_time.csv")
df_time_GEO$pttest = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/2_mpttest_time.csv")
df_time_GEO$mpttest = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/3_fcpttest_time.csv")
df_time_GEO$fcpttest = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/4_rpclr_time.csv")
df_time_GEO$`RP-CLR` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/5_clogitL1_time.csv")
df_time_GEO$`PCU-CLR` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/6_bayesian_time.csv")
df_time_GEO$`BVS-CLR` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/7_boostPQL_time.csv")
df_time_GEO$`1-Step PQLBoost` = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/8_boostWL2_time.csv")
df_time_GEO$WL2Boost = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/9_nonparied_mRMR_time.csv")
df_time_GEO$mRMR = df_tmp$elapsed

df_tmp = read.csv("/home/lsen/Data/CODE/feature_selection/_data/Feature_Selection/GEO-GSE70947/output/10_nonparied_MRMD_time.csv")
df_time_GEO$MRMD = df_tmp$elapsed

library(ggplot2)

mycolor = rainbow(10, alpha = 1)

df_time_GEO[is.na(df_time_GEO)] = Inf


par(xpd=F, mar=c(5,4,4,13))
# par(xpd=T, mar=c(5,4,4,13))

plot(x= NA,xlim = c(1,5),ylim = c(0,5000), xlab = "Five gene groups", xaxt='n',
     ylab="Runing times (s)", main = "GEO dataset running time comparison")
axis(1, at=c(1:5), label = c("10-genes", "50-genes", "100-genes", "1000-genes", "1500-genes"))

for (i in 1:10) {
  points(x = c(1:5), y=df_time_GEO[,i], pch=i, col=mycolor[i], lwd=2)
  lines(x = c(1:5), y=df_time_GEO[,i], pch=i, col=mycolor[i], lwd=2)
}

legend(5.2, 5000, legend = v_method_name, pch = c(1:10), col=mycolor, lwd=2)




