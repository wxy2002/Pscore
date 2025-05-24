library(tinyarray)
library(dplyr)
library(limma)
library(ggplot2)
library(edgeR)
library(ggpubr)
df <- read.csv('Score_TCGA.csv')
colnames(df) <- c('Pscore', 'Wscore', 'Gscore', 'Cscore')
df2 <- read.csv('Score_PLCO.csv')
colnames(df2) <- c('Pscore', 'Wscore', 'Gscore', 'Cscore')
df <- rbind(df, df2)
ggscatter(df, x = "Pscore", y = "Cscore",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray")) +
  stat_cor(method = "spearman")
cor.test(df$Pscore, df$Gscore, alternative = "two.side",
         method = "spearman", conf.level = 0.95, exact=FALSE)

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
group <- read.csv('Score_train2.csv', header = T, row.names = 1)
Group <- make_tcga_group(exp)
exp <- exp[, Group == 'tumor']
name <- colnames(exp)
name <- substr(name, 1, nchar(name) - 4)
colnames(exp) <- name
exp <- exp[, group$X_PATIENT]
group <- group[match(colnames(exp), group$X_PATIENT), ]
save(exp, group, proj, file = paste0(proj,".Rdata"))
exp2 <- as.data.frame(exp)
write.csv(exp2, file = "express.txt", row.names = T, quote = F)

#-------------------------
find_H <- function(name){
  load("TCGA-PRAD.Rdata")
  # name = 'Score'
  group[, name] <- factor(group[, name], 
                          levels = c('Low', 'High'))
  dge <- edgeR::DGEList(counts=exp)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~ group[, name])
  v <- voom(dge, design, normalize="quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  DEG <- topTable(fit, coef=2, n=Inf)
  DEG <- na.omit(DEG)
  library(tidyverse)
  library(data.table)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  DEG$SYMBOL <- rownames(DEG)
  gene <- bitr(rownames(DEG),
               fromType="SYMBOL",
               toType="ENTREZID",
               OrgDb="org.Hs.eg.db") %>%
    distinct(SYMBOL,.keep_all=TRUE)
  data_all <- DEG %>% 
    inner_join(gene, by="SYMBOL")
  data_all_sort <- data_all %>% 
    arrange(desc(logFC))
  head(data_all_sort)
  geneList <- data_all_sort$logFC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
  head(geneList)
  kegg_gmt <- read.gmt("F:/数据/富集基因集/MSigDB/c6.all.v2024.1.Hs.entrez.gmt")
  gsea <- GSEA(geneList,
               TERM2GENE = kegg_gmt) #GSEA分析
  head(gsea)
  return(gsea)
}
find_GO <- function(name){
  load("TCGA-PRAD.Rdata")
  # name = 'Score'
  group[, name] <- factor(group[, name], 
                          levels = c('Low', 'High'))
  dge <- edgeR::DGEList(counts=exp)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~ group[, name])
  v <- voom(dge, design, normalize="quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  DEG <- topTable(fit, coef=2, n=Inf)
  DEG <- na.omit(DEG)
  library(tidyverse)
  library(data.table)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  DEG$SYMBOL <- rownames(DEG)
  gene <- bitr(rownames(DEG),
               fromType="SYMBOL",
               toType="ENSEMBL",
               OrgDb="org.Hs.eg.db") %>%
    distinct(SYMBOL,.keep_all=TRUE)
  data_all <- DEG %>% 
    inner_join(gene, by="SYMBOL")
  data_all_sort <- data_all %>% 
    arrange(desc(logFC))
  head(data_all_sort)
  geneList <- data_all_sort$logFC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENSEMBL #给上面提取的foldchange加上对应上ENTREZID
  head(geneList)
  res <- gseGO(
    geneList,    # 根据logFC排序的基因集
    ont = "ALL",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
    OrgDb = org.Hs.eg.db,    # 使用人的OrgDb
    keyType = "ENSEMBL",    # 基因id类型
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",    # p值校正方法
  )
  return(res)
}
find_KEGG <- function(name){
  load("TCGA-PRAD.Rdata")
  name = 'Score'
  group[, name] <- factor(group[, name], 
                          levels = c('Low', 'High'))
  dge <- edgeR::DGEList(counts=exp)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~ group[, name])
  v <- voom(dge, design, normalize="quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  DEG <- topTable(fit, coef=2, n=Inf)
  DEG <- na.omit(DEG)
  library(tidyverse)
  library(data.table)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  DEG$SYMBOL <- rownames(DEG)
  gene <- bitr(rownames(DEG),
               fromType="SYMBOL",
               toType="ENTREZID",
               OrgDb="org.Hs.eg.db") %>%
    distinct(SYMBOL,.keep_all=TRUE)
  data_all <- DEG %>% 
    inner_join(gene, by="SYMBOL")
  data_all_sort <- data_all %>% 
    arrange(desc(logFC))
  head(data_all_sort)
  geneList <- data_all_sort$logFC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
  head(geneList)
  kegg_gmt <- read.gmt("F:/数据/富集参考基因集/MSigDB/c2.cp.kegg_medicus.v2024.1.Hs.entrez.gmt")
  gsea <- GSEA(geneList,
               TERM2GENE = kegg_gmt) #GSEA分析
  head(gsea)
  return(gsea)
}
wsi <- find_KEGG('Score')
gene <- find_KEGG('Score_Gene')
clin <- find_KEGG('Score_Clin2')
wsi <- wsi %>%
  filter(p.adjust < 0.05 & abs(NES) > 1 & qvalue < 0.25)
gene <- gene %>%
  filter(p.adjust < 0.05 & abs(NES) > 1 & qvalue < 0.25)
clin <- clin %>%
  filter(p.adjust < 0.05 & abs(NES) > 1 & qvalue < 0.25)
venn_list <- list(
  wsi = wsi$ID,
  gene = gene$ID,
  clin = clin$ID
)
library(GseaVis)
library(RColorBrewer)
wsi2 <- wsi[!(wsi$ID %in% gene$ID | wsi$ID %in% clin$ID), ] %>%
  mutate(Description = substr(Description, 10, nchar(Description)),
         dataset = 'WSI')
gene2 <- gene[!(gene$ID %in% wsi$ID | gene$ID %in% clin$ID), ] %>%
  mutate(Description = substr(Description, 10, nchar(Description)),
         dataset = 'Gene')
clin2 <- clin[!(clin$ID %in% wsi$ID | clin$ID %in% gene$ID), ] %>%
  mutate(Description = substr(Description, 10, nchar(Description)),
         dataset = 'Clin')
all <- rbind(wsi2, gene2, clin2) %>%
  select(ID, Description, NES, pvalue, dataset)
display.brewer.pal(6, 'Set1')
id <- wsi2$ID
gseaNb(object = wsi,
       geneSetID = id,
       newGsea = F,
       subPlot = 2,
       addPval = T, pvalX = 0.01, pvalY = -0.5,
       curveCol = brewer.pal(4, 'Set1'))
id <- gene2$ID
gseaNb(object = gene,
       geneSetID = id,
       newGsea = F,
       subPlot = 2,
       addPval = T, pvalX = 0.01, pvalY = -0.5,
       curveCol = brewer.pal(4, 'Set1'))
id <- clin2$ID
gseaNb(object = clin,
       geneSetID = id,
       newGsea = F,
       subPlot = 2,
       addPval = T, pvalX = 0.01, pvalY = -0.5,
       curveCol = brewer.pal(5, 'Set1'))

library(ggvenn)
ggvenn(
  data = venn_list,         # 数据列表
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)
# 筛选出wsi、gene、clin中独特的部分
# 去除wsi2 ID列中前9个字符

# ids <- gsea_res@result$ID[1:5]
gseadist(gsea_res,
         IDs = gsea_res2$ID,
         type = "density" # boxplot
) +
  theme(legend.direction = "vertical")
ridgeplot(gsea_res,
          showCategory = 10,
          fill = "pvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
) +
  theme(axis.text.y = element_text(size=10))

#-------------------------
# 基于exp，绘制某个基因不同分组表达箱线图
# 证明在前面的就是对照组
gene <- "UGT2B15"
load("TCGA-PRAD.Rdata")
exp <- as.data.frame(exp) %>%
  t() %>%
  as.data.frame()
exp$gene <- exp[, gene]
exp$group <- group$Score
exp$group <- factor(exp$group, levels = c('Low', 'High'))
exp$gene <- as.numeric(exp$gene)
exp$gene <- log(exp$gene + 1)
library(ggpubr)
p <- ggplot(exp, aes(x=group, y=gene)) +
  geom_boxplot(aes(fill=group), outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  labs(title = gene)
p

find_DEG <- function(name){
  load("TCGA-PRAD.Rdata")
  group[, name] <- factor(group[, name], 
                          levels = c('Low', 'High'))
  dge <- edgeR::DGEList(counts=exp)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~ group[, name])
  v <- voom(dge, design, normalize="quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  DEG <- topTable(fit, coef=2, n=Inf)
  DEG <- na.omit(DEG)
  logFC_t <- 1
  pvalue_t <- 0.05
  k1 = (DEG$adj.P.Val < pvalue_t) & (DEG$logFC < -logFC_t)
  k2 = (DEG$adj.P.Val < pvalue_t) & (DEG$logFC > logFC_t)
  DEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
  return(DEG)
}
wsi_DEG <- find_DEG('Score') %>%
  filter(change != 'NOT') %>%
  mutate(ID = rownames(.))
gene_DEG <- find_DEG('Score_Gene') %>%
  filter(change != 'NOT') %>%
  mutate(ID = rownames(.))
clin_DEG <- find_DEG('Score_Clin2') %>%
  filter(change != 'NOT') %>%
  mutate(ID = rownames(.))
venn_list <- list(
  wsi = wsi_DEG$ID,
  gene = gene_DEG$ID,
  clin = clin_DEG$ID
)

#-------------------------------
data <- read.table('TCGA-PRAD.star_fpkm.tsv.gz', header = T,
                   check.names = F, row.names = 1)
data <- data[, substr(colnames(data), nchar(colnames(data)), nchar(colnames(data))) == 'A']
colnames(data) <- substr(colnames(data), 1, nchar(colnames(data)) - 1)
data <- data[, substr(colnames(data), nchar(colnames(data))-1, nchar(colnames(data))) < '10']
colnames(data) <- substr(colnames(data), 1, nchar(colnames(data)) - 3)
data <- data[, !duplicated(colnames(data))]
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(msigdbr)
library(stringi)
data$ENSEMBL <- rownames(data)
data$ENSEMBL <- stri_sub(data$ENSEMBL, 1, 15)
data <- data[!duplicated(data$ENSEMBL), ]
gene_entrezid <- bitr(geneID = data$ENSEMBL, 
                      fromType = "ENSEMBL", 
                      toType = c("SYMBOL"),
                      OrgDb = "org.Hs.eg.db")
data2 <- read.csv('express.csv')
# gene_entrezid <- gene_entrezid[!duplicated(gene_entrezid$ENSEMBL), ]
data <- inner_join(data, gene_entrezid, by = "ENSEMBL")
data <- data[data$SYMBOL %in% data2$symbol, ]

library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
data3 <- data %>%
  select(SYMBOL, everything()) %>%
  select(-ENSEMBL)
data3 <- data3[!duplicated(data3$SYMBOL), ]
rownames(data3) <- data3$SYMBOL
data3 <- data3[,-1]
results <- cibersort(sig_matrix, data3, perm = 1000, QN = F)
head(results[,1:4],n=12)

load("TCGA-PRAD.Rdata")
res <- results[, 1:22] %>%
  data.frame() %>%
  mutate(X_PATIENT = rownames(.))
res <- inner_join(res, group, by = "X_PATIENT")
res1 <- res %>%
  pivot_longer(cols = colnames(.)[1:22],
               names_to = "cell.type",
               values_to = 'value')

library(ggpubr)
ggplot(res1, aes(cell.type, value, fill = Score)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  stat_compare_means(aes(group = Score), label = "p.signif", size=3,
                     method = "wilcox.test")
ggplot(res1, aes(cell.type, value, fill = Score_Gene)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  stat_compare_means(aes(group = Score_Gene), label = "p.signif", size=3,
                     method = "wilcox.test")
ggplot(res1, aes(cell.type, value, fill = Score_Clin2)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  stat_compare_means(aes(group = Score_Clin2), label = "p.signif", size=3,
                     method = "wilcox.test")
write.csv(res1, 'CIBERSORT_data.csv')
