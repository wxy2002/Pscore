# 读入数据
library(tinyarray)
library(dplyr)
library(limma)
library(ggplot2)
library(edgeR)
proj <- "TCGA-PRAD"
dat <- read.table("TCGA-PRAD.star_counts.tsv.gz", check.names= F, row.names = 1, header = T)
gtf_v22 <- rtracklayer::import('gencode.v22.annotation.gtf.gz') %>% as.data.frame()
gtf_mrna_v22 <- dplyr::select(gtf_v22, c("gene_id", "gene_type", "gene_name")) %>%
  subset(., gene_type == "protein_coding") %>%
  unique()
range(dat)
dat <- as.matrix(2^dat - 1)
exp <- round(dat)
exp <- trans_exp_new(exp)
exp[1:4,1:4]
exp <- exp[rownames(exp) %in% gtf_mrna_v22$gene_name, ]
nrow(exp)
exp <- exp[apply(exp, 1, function(x) sum(x > 0) > 0.8*ncol(exp)), ]
nrow(exp)
Group <- make_tcga_group(exp)
table(Group)
save(exp, Group, proj, file = paste0(proj,".Rdata"))

# TCGA数据集：limma进行差异分析
load("TCGA-PRAD.Rdata")
table(Group)
logFC_t <- 1
pvalue_t <- 0.05
dge <- edgeR::DGEList(counts=exp)
dge <- edgeR::calcNormFactors(dge)
design <- model.matrix(~ Group)
v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)
fit <- eBayes(fit)
DEG <- topTable(fit, coef=2, n=Inf)
DEG <- na.omit(DEG)
k1 = (DEG$adj.P.Val < pvalue_t) & (DEG$logFC < -logFC_t)
k2 = (DEG$adj.P.Val < pvalue_t) & (DEG$logFC > logFC_t)
DEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(DEG$change)
head(DEG)

exp[1:4,1:4]
dat = log2(cpm(exp)+1)
pca.plot = draw_pca(dat,Group);pca.plot
save(pca.plot,file = paste0(proj,"_pcaplot.Rdata"))

cg = rownames(DEG)[DEG$change !="NOT"]
v = draw_volcano(DEG, pkg = 3, logFC_cutoff = logFC_t, pvalue_cutoff = pvalue_t)
m = data.frame(m1 = colnames(dat), m2 = Group)
m = arrange(m,Group)
h = draw_heatmap(dat[cg, m$m1], m$m2, cluster_cols = F)
cg <- data.frame(Name = cg)
write.csv(DEG, 'Degs-TCGA.csv')

# GEO差异表达分析
logFC_t <- 1
pvalue_t <- 0.05
data_geo1 <- read.csv('GSE38241.top.table (2).tsv', sep = '\t')
k1 = (data_geo1$adj.P.Val < pvalue_t)&(data_geo1$logFC < -logFC_t)
k2 = (data_geo1$adj.P.Val < pvalue_t)&(data_geo1$logFC > logFC_t)
data_geo1$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
data_geo1$Gene.symbol <- data.frame(sapply(data_geo1$Gene.symbol, function(x)unlist(strsplit(x,"///"))[1]), stringsAsFactors=F)[,1]
data_geo1 <- na.omit(data_geo1)
table(data_geo1$change)
data_geo2 <- read.csv('GSE46602.top.table (1).tsv', sep = '\t')
k1 = (data_geo2$adj.P.Val < pvalue_t)&(data_geo2$logFC < -logFC_t)
k2 = (data_geo2$adj.P.Val < pvalue_t)&(data_geo2$logFC > logFC_t)
data_geo2$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
data_geo2$Gene.symbol <- data.frame(sapply(data_geo2$Gene.symbol, function(x)unlist(strsplit(x,"///"))[1]), stringsAsFactors=F)[,1]
data_geo2 <- na.omit(data_geo2)
table(data_geo2$change)

data_TCGA <- read.csv('Degs-TCGA.csv')
data_high <- read.csv('AUC75.csv')
genes <- list(GSE38241 = data_geo1[data_geo1$change == "UP", 'Gene.symbol'],
              GSE46602 = data_geo2[data_geo2$change == "UP", 'Gene.symbol'],
              `TCGA-PRAD` = data_TCGA[data_TCGA$change == "DOWN", 'X'],
              HighAUC = data_high$Gene)
genes2 <- list(GSE38241 = data_geo1[data_geo1$change == "DOWN", 'Gene.symbol'],
              GSE46602 = data_geo2[data_geo2$change == "DOWN", 'Gene.symbol'],
              `TCGA-PRAD` = data_TCGA[data_TCGA$change == "UP", 'X'],
              HighAUC = data_high$Gene)
p1 <- draw_venn(genes,"Up-regulated Genes")
p2 <- draw_venn(genes2,"Down-regulated Genes")
library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1)

name <- intersect(data_geo1[data_geo1$change == "DOWN", 'Gene.symbol'], data_geo2[data_geo2$change == "DOWN", 'Gene.symbol'])
name <- intersect(name, data_TCGA[data_TCGA$change == "UP", 'X'])
name <- intersect(name, data_high$Gene)
cg1 <- data.frame(name = name)
name <- intersect(data_geo1[data_geo1$change == "UP", 'Gene.symbol'], data_geo2[data_geo2$change == "UP", 'Gene.symbol'])
name <- intersect(name, data_TCGA[data_TCGA$change == "DOWN", 'X'])
name <- intersect(name, data_high$Gene)
cg2 <- data.frame(name = name)
cg <- rbind(cg1, cg2)

# 预后数据处理
data <- data.frame(dat, check.names= F)
data <- rbind(data, Group)
data <- data[cg$name, ]
data <- t(data) %>%
  data.frame()
clin <- read.csv('survival_PRAD_survival.txt', sep = '\t') %>%
  select(sample, X_PATIENT, PFI, PFI.time)
library(tidyverse)
library(survival)
s <- str_sub(rownames(data), start = 1, end = -2)
data$sample <- s
data <- left_join(data, clin, by = 'sample')
write.csv(data, 'Gene_cpm.csv')

train <- data %>%
  select(-sample)
train <- train[!duplicated(train$X_PATIENT), ] 
rownames(train) <- train$X_PATIENT
x <- train %>%
  select(-PFI, -PFI.time, -X_PATIENT) %>%
  as.matrix()
y <- train[, c('PFI', 'PFI.time')]
save(exp, x, y, train, file = paste0(proj,"-PFI.Rdata"))

# Lasso回归筛选
library(glmnet)
mod <- glmnet(x, Surv(y$PFI.time, y$PFI), 
              type.measure = "deviance", family = "cox")
plot(mod, xvar = "lambda", label = T, lwd=2)
set.seed(42)
lasso_fit <- cv.glmnet(x, Surv(y$PFI.time, y$PFI), family = 'cox', 
                       type.measure = 'deviance', nfolds = 10)
plot(lasso_fit)
lambda <- lasso_fit$lambda.min
model_lasso <- glmnet(x, Surv(y$PFI.time, y$PFI), family = 'cox',
                      type.measure = 'deviance', lambda = lambda)
gene_name <- rownames(model_lasso$beta)[as.numeric(model_lasso$beta)!=0]
gene_name
paste(gene_name, collapse = ' + ')

cox_model_fin <- coxph(Surv(PFI.time, PFI) ~ 
                         CHTF18 + CDCA5 + TRPM4 + SERPINB5 + ARHGEF40 + LIX1, data = train)
save(cox_model_fin, file = 'model.rdata')
summary(cox_model_fin)
rms::vif(cox_model_fin)
library(survminer)
ggforest(cox_model_fin,
         main = "Hazard ratio", data = train)
train$Gene_Score <- predict(cox_model_fin, train)
library(survminer)
name = 'Gene_Score'
train[, name] <- ifelse(train[, name] >= median(train$Gene_Score), 'High', 'Low')
fit <- survfit(Surv(PFI.time, PFI) ~ train[, name], data = train)
fit
p <- ggsurvplot(fit, data = train,
                conf.int = TRUE, risk.table = TRUE,
                legend.labs = c("High Score", "Low Score"))
data.survdiff <- survdiff(Surv(PFI.time, PFI) ~ train[, name], data = train)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")
p$plot <- p$plot + 
  annotate("text", x = 0, y = 0.1, 
           label = ifelse(p.val < 0.0001, 
                          paste0(" p < 0.0001","\n HR (95% CI) = ",ci), 
                          paste0(" P value = ",sprintf("%.3f",p.val),"\n HR (95% CI) = ",ci)), 
           size = 5, color = "black", hjust = 0) +
  theme(text = element_text(size = 15))
p

clin <- read.csv('Score_train.csv', row.names = 1)
gene_pre <- read.csv('gene_express_data.csv')
cox_pre <- coxph(Surv(PFI.time, PFI) ~ 
                   CHTF18 + CDCA5 + TRPM4 + SERPINB5 + ARHGEF40 + LIX1, data = gene_pre)
summary(cox_pre)
cox.zph(cox_pre)
sqrt(rms::vif(cox_pre)) < 2
gene_pre$Score_gene <- predict(cox_pre, gene_pre)
gene_pre <- inner_join(gene_pre, clin, by = 'X_PATIENT')
summary(coxph(Surv(PFI.time.x, PFI.x) ~ Score + Score_Clin + Score_gene, data = gene_pre))

wsi_Score <- read.csv("Score_TCGA_wsi.csv")
p <- predict(step_wise_model, train)
train$GeneScore <- p
train <- left_join(train, wsi_Score, by = 'X_PATIENT')
c <- coxph(Surv(PFI.time, PFI) ~ GeneScore, data = train)
summary(c)

# AUC75（显著的）
SERPINB5 + CHTF18 + ARHGEF40 + CDCA5 + LIX1
C_G = 0.732, C_G/W = 0.749, C_W = 0.691
# AUC80（显著的）
SRD5A2 + CHTF18 + ARHGEF40 + LIX1 + PTGS1
C_G = 0.727, C_G/W=0.746, C_W = 0.691
# AUC85
LIX1
C_G = 0.642, C_G/W=0.708, C_W = 0.691

list <- c('SERPINB5', 'CHTF18', 'ARHGEF40', 'CDCA5', 'LIX1')
p1 <- data_TCGA %>%
  filter(X %in% list)
p2 <- data_geo1 %>%
  filter(Gene.symbol %in% list) %>%
  distinct(Gene.symbol, .keep_all = T)
p3 <- data_geo2 %>%
  filter(Gene.symbol %in% list) %>%
  distinct(Gene.symbol, .keep_all = T)
