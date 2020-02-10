# ****************************************************
# Email: lsen@infervision.com
# Data: Sun Jan 14 14:35:22 2018
# Description: null
#                
# *****************************************************
  
require(Biobase)
library(GEOquery)

gse2 <- getGEO('GSE70947',GSEMatrix=TRUE, destdir = "./GEO/")
save(gse2, file = "./_data/GSE70947")

gse = gse$GSE70947_series_matrix.txt.gz

# gse2 <- getGEO('GSE70947')
# save(gse, file = "./_data/GSE70947")

# 提取三类数据
data_clin = pData(gse)                             # 临床信息
data_expr = as.data.frame(assayData(gse)$exprs)    # 基因表达信息
data_gene = featureData(gse)@data                  # 基因probe信息





gse <- getGEO(filename="/home/lsen/GSE70947_family.soft.gz")
save(gse, file = "./_data/GSE70947.rdata")
head(Meta(gse))

gsm_list = GSMList(gse)
gpl_list = GPLList(gse)

names(gsm_list)
names(gpl_list)

# make sure all the data 
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)

gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL13607'},GSMList(gse))
length(gsmlist)

Table(gsmlist[[1]])[1:5,]
Columns(gsmlist[[1]])[1:5,]


# get the probeset ordering
probesets <- Table(GPLList(gse)[[1]])$ID

# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',
                       lapply(gsmlist,function(x){ 
                         tab <- Table(x) 
                         mymatch <- match(probesets,tab$ID_REF)
                         return(tab$VALUE[mymatch])}))

data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
# data.matrix <- log2(data.matrix)
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)

# 获得数据
data_expr = data.frame(data.matrix)                # 基因表达信息
data_gene = gpl_list$GPL13607@dataTable@table      # 临床信息
data_clin = featureData(gse)@data                  # 基因probe信息

# 筛选基因
# (1) 去掉Control的probe
table(data_gene$ControlType)
data_gene2 = data_gene[which(data_gene$ControlType == 0),]

# (2) 去掉lincRNA
library(stringr)
data_gene2$accessions = str_trim(data_gene2$accessions)
data_gene3 = data_gene2[which(data_gene2$SPOT_ID != "lincRNA"),]

length(unique(data_gene3$GeneName))

### 用上面筛选出来的gene去filter data_expr
data_expr2 = data_expr[rownames(data_expr) %in% data_gene3$ID, ]
View(data_expr2[1:1000, 1:20])

### 先去掉重复基因，然后再去掉NA缺失的情况
# 重复的probe取平均值
# 1-找出有重复的probe
data_gene3$GeneName = str_trim(data_gene3$GeneName)

v_copy = duplicated(data_gene3$GeneName)
v_copy_gene = unique(data_gene3$GeneName[v_copy])

v_all_gene = unique(data_gene3$GeneName)
v_all_gsm = colnames(data_expr2)

data_expr3 = data.frame(matrix(NA, nrow = length(v_all_gene), ncol = length(v_all_gsm)), stringsAsFactors = F)
colnames(data_expr3) = v_all_gsm
rownames(data_expr3) = v_all_gene

# 将所有的基因第一次出现的那个放到data.frame中
for( i in 1:length(v_all_gene)){
  # i = 2
  igene = v_all_gene[i]
  v_id = data_gene3$ID[which(data_gene3$GeneName == igene)][1]
  
  data_expr3[i, ] = data_expr2[which(rownames(data_expr2) == v_id),]
} # 很地效的方式，换种方式吧

# 将重复的基因重新更新
for( i in 1:length(v_copy_gene)){
  #i = 2
  igene = v_copy_gene[i]
  v_id = data_gene3$ID[which(data_gene3$GeneName == igene)]
  
  
  tmp_expr = data_expr2[which(rownames(data_expr2) %in% v_id),]
  
  tmp_mean = colMeans(tmp_expr, na.rm = T)
  tmp_mean[is.nan(tmp_mean)] = NA
  
  data_expr3[which(v_all_gene == igene), ] = tmp_mean
}

save.image(file = "S4-GEO-2018018.rdata")

### 统计NA的缺失率,并补充缺失的数据
# 1- 去掉缺失较大的基因(去掉缺失率>1%的数据)
sum(is.na(data_expr3[1,]))

v_sum_na = sapply(1:nrow(data_expr3), function(i){sum(is.na(data_expr3[i,]))}, simplify = T)
v_sum_na_rate = v_sum_na / (1.0*ncol(data_expr3))

hist(v_sum_na_rate)
sum(v_sum_na_rate < 0.01)

data_expr4 = data_expr3[ which(v_sum_na_rate < 0.01),]

# 2- 去掉基因缺失率较大的Samples
v_sum_na2 = sapply(1:ncol(data_expr4), function(i){sum(is.na(data_expr4[,i]))}, simplify = T)
v_sum_na_rate2 = v_sum_na2 / (1.0*nrow(data_expr4))

hist(v_sum_na_rate2)
ncol(data_expr4)
sum(v_sum_na_rate2 < 0.01)

data_expr5 = data_expr4[,which(v_sum_na_rate2 < 0.01)]

v_sample1 = colnames(data_expr5) # 剩下291个数据了

# 3- 划分成Tumor/Normal, 然后分别对Normal/Tumor去用均值补充NA数据
table(data_clinical$source_name_ch1)
colnames(data_clinical)

data_clinical$title = str_trim(as.character(data_clinical$title))
data_clinical$geo_accession = str_trim(as.character(data_clinical$geo_accession))

v_title = str_split(data_clinical$title, "-", simplify = T)
data_clinical$sampleID = v_title[,1]
data_clinical$tumorOrNormal = v_title[,2]


# 选择出v—sample1中的数据
data_clinical2 = data_clinical[which(data_clinical$geo_accession %in% v_sample1),]


v_sampelID = unique(data_clinical2$sampleID)
v_sampleTable = table(data_clinical2$sampleID)
v_sampleTable = v_sampleTable[v_sampleTable ==2]
v_sampleID = names(v_sampleTable)

# 选择出最终的数据集
data_clinical2 = data_clinical2[which(data_clinical2$sampleID %in% v_sampleID),]
data_expr6 = data_expr5[, which(colnames(data_expr5) %in% data_clinical2$geo_accession)]

# 接下来划分出case/ctrol
data_clinical2 = data_clinical2[order(data_clinical2$sampleID),]
v_sample_case = data_clinical2$geo_accession[which(data_clinical2$tumorOrNormal == "tumor")]
v_sample_ctrol = data_clinical2$geo_accession[which(data_clinical2$tumorOrNormal == "normal")]

data_case = data_expr6[,v_sample_case]
data_ctrl = data_expr6[,v_sample_ctrol]

save(data_ctrl, data_case, data_clinical2, file = "./_data/data_case_ctrol_GSE70947.rdata")



# 2- 用均值法补气NA数据
# 对case而言
v_gene_mean_case = rowMeans(data_case, na.rm = T)

v_have_na = sapply(1:nrow(data_case), function(i){sum(is.na(data_case[i,]))}, simplify = T)
v_have_na2 = which(v_have_na > 0)
for(ind in v_have_na2){
  data_case[ind, is.na(data_case[ind,])] = v_gene_mean_case[ind]
}

sum(is.na(data_case))

# 对control而言
v_gene_mean_ctrl = rowMeans(data_ctrl, na.rm = T)
v_have_na = sapply(1:nrow(data_ctrl), function(i){sum(is.na(data_ctrl[i,]))}, simplify = T)
v_have_na2 = which(v_have_na > 0)

for(ind in v_have_na2){
  data_ctrl[ind, is.na(data_ctrl[ind,])] = v_gene_mean_ctrl[ind]
}

sum(is.na(data_ctrl))

## 整理成和TCGA一样的数据
save(data_ctrl, data_case, data_clinical, file = "./_data/data_case_ctrol_GSE70947.rdata")
save.image(file = "S4-GEO-20180117.rdata")
