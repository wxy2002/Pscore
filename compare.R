expression = read.csv('TCGA-PRAD.star_tpm.tsv', sep = '\t', row.names = 1)
library(stringi)
expression$ID = rownames(expression)
ensg2symbol = read.csv('gencode.v36.annotation.gtf.gene.probemap', sep = '\t')
# 根据expression$ID，找到ensg2symbol中对应的symbol
expression = merge(expression, ensg2symbol[, c('id', 'gene')], by.x = 'ID', by.y = 'id')
# 去除重复基因，保留第一个基因
expression = expression[!duplicated(expression$gene), ]
# 转置为行为样本，列为基因
rownames(expression) = expression$gene
expression = expression[, -which(colnames(expression) == 'ID')]
expression = expression[, -which(colnames(expression) == 'gene')]
expression = as.data.frame(t(expression))
expression$SampleID = rownames(expression)
# 仅保留最后三个字符为01A的样本
expression = expression[which(stri_sub(expression$SampleID, -3, -1) == '01A'), ]
# 提取前12个字符作为样本ID
expression$SampleID = stri_sub(expression$SampleID, 1, 12)
# 把SampleID中的.改为-
expression$SampleID = gsub('\\.', '-', expression$SampleID)

clin = read.csv('survival_PRAD_survival.txt', sep = '\t')
expression = merge(expression, clin[, c('X_PATIENT', 'PFI', 'PFI.time')], by.x = 'SampleID', by.y = 'X_PATIENT')
# 去重
expression = expression[!duplicated(expression$SampleID), ]




# A novel prognostic model for prostate cancer based on androgen biosynthetic and catabolic pathways
# 公式
# (-0.1028 * AFF3) + (0.2922 * B4GALNT4) + (-0.511 * CD38) + (0.0072 * CHRNA2) +
# (0.038318 * CST2) + (0.2041 * ADGRF5) + (0.1326 * KLK14) + (0.1493 * LRRC31) + 
# (-0.1262 * MT1F) + (-0.0116 * MT1G) + (-0.1639 * SFTPA2) + (-0.1697 * SLC7A4) + (-0.4706 * TDRD1)

score = expression$AFF3 * -0.1028 + expression$B4GALNT4 * 0.2922 + expression$CD38 * -0.511 + expression$CHRNA2 * 0.0072 +
  expression$CST2 * 0.038318 + expression$ADGRF5 * 0.2041 + expression$KLK14 * 0.1326 + expression$LRRC31 * 0.1493 +
  expression$MT1F * -0.1262 + expression$MT1G * -0.0116 + expression$SFTPA2 * -0.1639 + expression$SLC7A4 * -0.1697 +
  expression$TDRD1 * -0.4706
expression$score = score
# 计算C指数
library(survival)
library(tidyverse)
fit = coxph(Surv(PFI.time, PFI) ~ score, data = expression)
summary(fit)
df = expression %>%
  select(SampleID, score)
write.csv(df, 'TCGA-PRAD.score.36439479.csv', row.names = FALSE)



# mtPCDI: a machine learning-based prognostic model for prostate cancer recurrence
# 筛选出的基因：PIK3R1、MSRB3、PTGIS、TPM1、HSPB8、EYA4、PAK3、GGCT
library(tidyverse)
gene_expression <- expression %>%
  select(SampleID, PIK3R1, MSRB3, PTGIS, TPM1, HSPB8, EYA4, PAK3, GGCT, PFI, PFI.time)
library(survival)
library(gbm)
set.seed(42)
ind <- sample(1:nrow(gene_expression), size = 0.7*nrow(gene_expression))
train <- gene_expression[ind,]
test <- gene_expression[-ind,]
gbm1 <- gbm(
  formula = Surv(PFI.time, PFI) ~ . - SampleID - PFI.time - PFI,
  distribution = "coxph",
  data = train,
  n.trees = 100,
  interaction.depth = 3,
  shrinkage = 0.01,
  n.minobsinnode = 10,
  verbose = FALSE
)
predicted_risk <- predict(gbm1, newdata = gene_expression, n.trees = 100)
gene_expression$mtPCDI_score <- predicted_risk
fit_mtp = coxph(Surv(PFI.time, PFI) ~ mtPCDI_score, data = gene_expression)
summary(fit_mtp)
# 保存结果
write.csv(gene_expression %>% select(SampleID, mtPCDI_score), 'TCGA-PRAD.mtPCDI.score.csv', row.names = FALSE)

# CAPRA Score
clin_all = read.csv('TCGA-PRAD.clinical.tsv', sep = '\t', row.names = 1)
# 仅保留需要的行
library(tidyverse)
clin_all <- clin_all %>%
  select(age_at_diagnosis.diagnoses, primary_gleason_grade.diagnoses,
         secondary_gleason_grade.diagnoses, ajcc_clinical_t.diagnoses, submitter_id)
# 空白值转为NA
clin_all[clin_all == ''] <- NA
# 去除有NA的行
clin_all <- na.omit(clin_all)
clin_all$rowname = rownames(clin_all)
clin_all = clin_all[which(stri_sub(clin_all$rowname, -3, -1) == '01A'), ]

clin_all$primary_gleason_grade.diagnoses = as.numeric(stri_sub(clin_all$primary_gleason_grade.diagnoses, -1, -1))
clin_all$secondary_gleason_grade.diagnoses = as.numeric(stri_sub(clin_all$secondary_gleason_grade.diagnoses, -1, -1))
colnames(clin_all)
age = clin_all$age_at_diagnosis.diagnoses
# 小于50为0，否则为1
age_score = ifelse(age < 50, 0, 1)
gleason_p = clin_all$primary_gleason_grade.diagnoses
gleason_s = clin_all$secondary_gleason_grade.diagnoses
gleason_score = ifelse((gleason_p<4 & gleason_s<4), 0,
                       ifelse(gleason_p >= 4, 3, 1))
clin_stage = clin_all$ajcc_clinical_t.diagnoses
# T1或T2为0，其余为1
clin_score = ifelse(clin_stage %in% c('T1', 'T2'), 0, 1)
clin_all$capra_score1 = age_score + gleason_score + clin_score

clin_all_2 = read.csv('TCGA.PRAD.sampleMap_PRAD_clinicalMatrix', sep = '\t', row.names = 1)
clin_all_2 <- clin_all_2 %>%
  select(psa_value, X_PATIENT)
clin_all_2$rowname = rownames(clin_all_2)
clin_all_2 = clin_all_2[which(stri_sub(clin_all_2$rowname, -2, -1) == '01'), ]
# 空白值转为NA
clin_all_2[clin_all_2 == ''] <- NA
# 去除有NA的行
clin_all_2 <- na.omit(clin_all_2)
psa = clin_all_2$psa_value
psa_score = ifelse(psa <=6, 0, ifelse(psa <=10, 1, ifelse(psa <=20, 2, ifelse(psa <=30, 3, 4))))
clin_all_2$psa_score = psa_score

clin_all = merge(clin_all, clin_all_2[, c('X_PATIENT', 'psa_score')], by.x = 'submitter_id', by.y = 'X_PATIENT')
clin_all$capra_score = clin_all$capra_score1 + clin_all$psa_score
expression = merge(expression, clin_all[, c('submitter_id', 'capra_score1', 'capra_score')], by.x = 'SampleID', by.y = 'submitter_id')
fit2 = coxph(Surv(PFI.time, PFI) ~ capra_score, data = expression)
summary(fit2)

