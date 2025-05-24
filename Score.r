library(tidyverse)
library(survival)
library(Hmisc)
library(pROC)
library(caret)
library(survminer)
library(caret)
library(compareGroups)

time <- function(data, cph){
  library(timeROC)
  p <- cph
  ROC <- timeROC(T = data$PFI.time, #生存时间
                 delta = data$PFI,   #生存状态
                 marker = p, #计算timeROC的变量
                 cause = 1,
                 weighting = "marginal",
                 times = c(1 * 365, 3 * 365, 5 * 365, 10 * 365),
                 iid=TRUE)
  return(ROC)
}
mult <- function(p_all, p_WSI, p_Gene, p_Clin, time){
  library(ggplot2)
  label <- case_when(
    time == 365 ~ 1,
    time == 1095 ~ 2,
    time == 1825 ~ 3
  )
  AUC <- c(p_all$AUC[label], p_WSI$AUC[label], p_Gene$AUC[label], p_Clin$AUC[label])
  data_all <- data.frame(
    FP = p_all$FP[, label],
    TP = p_all$TP[, label]
  )
  data_WSI <- data.frame(
    FP = p_WSI$FP[, label],
    TP = p_WSI$TP[, label]
  )
  data_Gene <- data.frame(
    FP = p_Gene$FP[, label],
    TP = p_Gene$TP[, label]
  )
  data_Clin <- data.frame(
    FP = p_Clin$FP[, label],
    TP = p_Clin$TP[, label]
  )
  ggplot() + 
    geom_line(data = data_all, aes(x = FP, y = TP, color = '1'), size = 1) +
    geom_line(data = data_WSI, aes(x = FP, y = TP, color = '2'), size = 1) +
    geom_line(data = data_Gene, aes(x = FP, y = TP, color = '3'), size = 1) +
    geom_line(data = data_Clin, aes(x = FP, y = TP, color = '4'), size = 1) +
    geom_line(aes(x = c(0,1), y = c(0,1)), color = "grey")+
    theme_bw() + 
    labs(x = "1 - Specificity",
         y = "Sensitivity") +
    scale_color_manual(name = NULL,
                       breaks = c('1', '2', '3', '4'),
                       values = c('1' = "#92C5DE", '2' = "#F4A582", '3' = "#66C2A5", '4' = 'pink'),
                       labels = paste0(c('All Model', 'WSI Model', 'Gene Model', 'Clin Model'), " : ",
                                       format(round(AUC, 2), nsmall = 2))) +
    scale_x_continuous(expand = c(0.005,0.005)) +
    scale_y_continuous(expand = c(0.005,0.005)) +
    coord_fixed() + 
    ggtitle(paste0(time, ' days AUC of Models'))
}

name <- read.csv('GeneScore（先筛选）/gene_express_data.csv')
train <- read.csv('GeneScore（先筛选）/gene_express_data.csv') %>%
  select(CPLX1, FCGR2A, FOXS1, MGAT4B, SUGP2, TPX2, U2AF1L4, ASPN, FTH1, MXD3, RBIS,
         ABCB6, ALDOA, CPLX1, FAM193B, GRIPAP1, MC1R, TRAIP, PRCC, VPS36, 
         SERPINB5, CHTF18, ARHGEF40, CDCA5, LIX1, TRPM4, SRD5A2, PFI.time, PFI, OS, OS.time,
         DSS, DSS.time, X_PATIENT)
train <- train[!duplicated(train$X_PATIENT), ]
# data <- read.csv('GeneScore/train.csv') %>%
#   select(EVX2, CTF1, PFDN6, S100A6, PFI.time, PFI, X_PATIENT)
clin <- read.table("临床资料/2024更新（和2019相同）/TCGA-PRAD.clinical.tsv",
                   sep = '\t', header = T)
clin2 <- read.table("D:/科研/WSI-2024/PRAD/临床资料/2019更新/TCGA.PRAD.sampleMap_PRAD_clinicalMatrix",
                    header = T, sep = '\t') %>%
  select(X_PATIENT, psa_value, age_at_initial_pathologic_diagnosis)
clin2$logPSA <- log(clin2$psa_value + 1)
clin2$age <- clin2$age_at_initial_pathologic_diagnosis
gleason <- clin %>%
  select(submitter_id, primary_gleason_grade.diagnoses, 
         secondary_gleason_grade.diagnoses, ajcc_pathologic_n.diagnoses, 
         ajcc_pathologic_t.diagnoses, ajcc_clinical_t.diagnoses, ajcc_clinical_m.diagnoses)
colnames(gleason) <- c('X_PATIENT', 'Primary', 'Secondary', 'pN', 'pT', 'cT', 'cM')
gleason <- gleason %>%
  distinct(X_PATIENT, .keep_all = T)
clin2 <- clin2 %>%
  distinct(X_PATIENT, .keep_all = T)
gleason <- inner_join(gleason, clin2, by = 'X_PATIENT')
gleason$Primary <- as.numeric(str_split(gleason$Primary, ' ', simplify = T)[, 2])
gleason$Secondary <- as.numeric(str_split(gleason$Secondary, ' ', simplify = T)[, 2])
gleason$Gleason <- gleason$Primary + gleason$Secondary
gleason$pN <- factor(gleason$pN)
gleason <- gleason %>% mutate(
  NCCN = case_when(
    cT <= 'T2a' & Gleason <= 6 & psa_value < 10 ~ 1,
    cT >= 'T3a' | Gleason >= 8 | (psa_value > 20) ~ 3,
    TRUE ~ 2
  )
)
gleason <- gleason %>% mutate(
  `Clincal T` = case_when(
    cT < 'T2' ~ 'T1',
    cT < 'T3' ~ 'T2',
    cT < 'T4' ~ 'T3',
    cT == 'T4' ~ 'T4',
    TRUE ~ 'Unknown'
  )
) %>%
  mutate(cT = case_when(
  cT < 'T2' ~ 1,
  cT < 'T3' ~ 2,
  cT < 'T4' ~ 3,
  cT == 'T4' ~ 4,
  TRUE ~ 0
))
head(gleason)

clin3 <- read.table('D:/科研/WSI-2024/PRAD/临床资料/clinical.tsv', sep = '\t',
                    fill = T, header = T)

Score_WSI <- read.csv("Score_TCGA_wsi.csv")
Score_train <- train
Score_train <- left_join(Score_train, gleason, by = 'X_PATIENT')
Score_train <- left_join(Score_train, Score_WSI, by = 'X_PATIENT')
Score_train2 <- Score_train %>%
  mutate(Primary = as.character(Primary),
         Secondary = as.character(Secondary),
         Gleason = as.character(Gleason),
         PFI = as.character(PFI),
         OS = as.character(OS),
         PCSM = as.character(DSS),
         NCCN = as.character(NCCN))
p <- descrTable(~ age + psa_value + Primary + Secondary + Gleason + 
                  `Clincal T` + PFI + OS + PCSM + OS.time + NCCN, 
                data = Score_train2,
                method = c(age = 2, psa_value = 2,
                           OS.time = 2))
export2xls(p, file='TCGA_Clin_Table.xlsx')

cph_G <- coxph(Surv(PFI.time, PFI) ~ 
                 CHTF18 + CDCA5 + TRPM4 + SERPINB5 + ARHGEF40 + LIX1, data = Score_train)
library(rms)
vif(cph_G)
summary(cph_G)
cox.zph(cph_G)
ggcoxzph(cox.zph(cph_G))
ggcoxdiagnostics(cph_G, 
                 type = "deviance",
                 ox.scale = "linear.predictions")
# sqrt(rms::vif(cph_G)) < 2
cph_C <- coxph(Surv(PFI.time, PFI) ~ Primary + Secondary + cT, data = Score_train)
summary(cph_C)
vif(cph_C)
cox.zph(cph_C)
# sqrt(rms::vif(cph_C)) < 2
cph_C2 <- coxph(Surv(PFI.time, PFI) ~ Gleason + cT, data = Score_train)
summary(cph_C2)
cph_W <- coxph(Surv(PFI.time, PFI) ~ Score, data = Score_train)
summary(cph_W)

Score_train$Score_Gene <- 2.09636*Score_train$CHTF18 - 
  0.11524*Score_train$CDCA5 - 1.00421*Score_train$TRPM4 + 
  0.02448*Score_train$SERPINB5 + 0.51506*Score_train$ARHGEF40 - 
  0.87403*Score_train$LIX1
Score_train$Score_Clin <- 0.8227*Score_train$Primary +
  0.7457*Score_train$Secondary + 0.4466*Score_train$cT
Score_train$Score_Clin2 <- 0.7759*Score_train$Gleason + 0.4513*Score_train$cT
cph_all <- coxph(Surv(PFI.time, PFI) ~ Score_Gene + Score + Score_Clin, data = Score_train)
summary(cph_all)
vif(cph_all)
ggplot(Score_train, aes(x = Score)) + 
  geom_density(alpha = 0.3,
               fill = 'pink', lwd = 1)
cox.zph(cph_all)
ggcoxzph(cox.zph(cph_all))
cph_all2 <- coxph(Surv(PFI.time, PFI) ~ Score_Gene + Score + Score_Clin2, data = Score_train)
summary(cph_all2)
cox.zph(cph_all2)
ggcoxzph(cox.zph(cph_all))
# sqrt(rms::vif(cph_all))

cph_GW <- coxph(Surv(PFI.time, PFI) ~ Score_Gene + Score , data = Score_train)
summary(cph_GW)
cph_WC <- coxph(Surv(PFI.time, PFI) ~ Score + Score_Clin , data = Score_train)
summary(cph_WC)
cph_GC <- coxph(Surv(PFI.time, PFI) ~ Score_Gene + Score_Clin , data = Score_train)
summary(cph_GC)

Score_train$Score_All <- 0.4320 * Score_train$Score_Gene + 
  1.0939 * Score_train$Score + 0.7207 * Score_train$Score_Clin
Score_train$Score_All2 <- 0.4332 * Score_train$Score_Gene + 
  1.1055 * Score_train$Score + 0.7224 * Score_train$Score_Clin2
summary(coxph(Surv(PFI.time, PFI) ~ Score_All2, data = Score_train))

# Score_train$Score_Gene <- scale(Score_train$Score_Gene)
# Score_train$Score_Clin <- scale(Score_train$Score_Clin)
write.csv(Score_train, 'Score_train.csv')

save(cph_all, cph_all2, cph_C, cph_G, cph_W, cph_GW, 
     cph_WC, cph_GC, cph_C2, file = 'D:/科研/WSI-2024/PLCO/model.rdata')
summary(cph_all)

clin_var <- c('age_at_initial_pathologic_diagnosis', 'psa_value', 'Gleason', 'cT', 'Score_All2')
clin_var <- data.frame(clin_var)
clin_var <- clin_var %>%
  mutate(n = NA, HR = NA, `95% CI` = NA, p = NA)
Score_train$Gleason_2 <- ifelse(Score_train$Gleason >= 7, 1, 0)
for (i in 1:length(clin_var[, 1])){
  if (clin_var[i, 'clin_var'] == 'age_at_initial_pathologic_diagnosis'){
    reg <- coxph(Surv(PFI.time, PFI) ~ age_at_initial_pathologic_diagnosis, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'psa_value'){
    reg <- coxph(Surv(PFI.time, PFI) ~ log2(psa_value + 1), data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'Gleason'){
    reg <- coxph(Surv(PFI.time, PFI) ~ Gleason, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'cT'){
    reg <- coxph(Surv(PFI.time, PFI) ~ cT, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'Score_All2'){
    reg <- coxph(Surv(PFI.time, PFI) ~ Score_All, data = Score_train)
  }
  sumreg <- summary(reg)
  clin_var[i, 'n'] <- round(sumreg$n, 4)
  clin_var[i, 'HR'] <- round(sumreg$coef[2], 4)
  clin_var[i, '95% CI'] <- paste0(round(sumreg$conf.int[3], 4), '-', round(sumreg$conf.int[4], 4))
  clin_var[i, 'p'] <- round(sumreg$sctest[3], 4)
}
clin_var
reg <- coxph(Surv(PFI.time, PFI) ~ 
               Score_All + age_at_initial_pathologic_diagnosis + log2(psa_value + 1) + Gleason_2 + cT,
             data = Score_train)
sumreg <- summary(reg)
p <- round(sumreg$coef, 4) %>% as.data.frame()
p$`95% CI` <- paste0(round(sumreg$conf.int[, 3], 4), '-', round(sumreg$conf.int[, 4], 4))
p$n <- round(sumreg$n, 4)
p$p <- round(sumreg$coef[, 5], 4)
p <- p %>%
  select(n, `exp(coef)`, `95% CI`, p)

clin_var <- c('age_at_initial_pathologic_diagnosis', 'psa_value', 'Gleason', 'cT', 'pN', 'Primary', 'Secondary')
clin_var <- data.frame(clin_var)
clin_var <- clin_var %>%
  mutate(n = NA, HR = NA, `95% CI` = NA, p = NA)
Score_train$pN <- ifelse(Score_train$pN=='N1', 1, ifelse(Score_train$pN=='N0', 0, NA))
Score_train$Gleason_2 <- ifelse(Score_train$Gleason >= 7, 1, 0)
for (i in 1:length(clin_var[, 1])){
  if (clin_var[i, 'clin_var'] == 'age_at_initial_pathologic_diagnosis'){
    reg <- coxph(Surv(PFI.time, PFI) ~ age_at_initial_pathologic_diagnosis, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'psa_value'){
    reg <- coxph(Surv(PFI.time, PFI) ~ log2(psa_value + 1), data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'Gleason'){
    reg <- coxph(Surv(PFI.time, PFI) ~ Gleason, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'cT'){
    reg <- coxph(Surv(PFI.time, PFI) ~ cT, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'pN'){
    reg <- coxph(Surv(PFI.time, PFI) ~ pN, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'Primary'){
    reg <- coxph(Surv(PFI.time, PFI) ~ Primary, data = Score_train)
  } else if (clin_var[i, 'clin_var'] == 'Secondary'){
    reg <- coxph(Surv(PFI.time, PFI) ~ Secondary, data = Score_train)
  }
  sumreg <- summary(reg)
  clin_var[i, 'n'] <- round(sumreg$n, 4)
  clin_var[i, 'HR'] <- round(sumreg$coef[2], 4)
  clin_var[i, '95% CI'] <- paste0(round(sumreg$conf.int[3], 4), '-', round(sumreg$conf.int[4], 4))
  clin_var[i, 'p'] <- round(sumreg$sctest[3], 4)
}
clin_var
reg <- coxph(Surv(PFI.time, PFI) ~ 
               age_at_initial_pathologic_diagnosis + log2(psa_value + 1) + Gleason + cT + pN,
             data = Score_train)
sumreg <- summary(reg)
p <- round(sumreg$coef, 4) %>% as.data.frame()
p$`95% CI` <- paste0(round(sumreg$conf.int[, 3], 4), '-', round(sumreg$conf.int[, 4], 4))
p$n <- round(sumreg$n, 4)
p$p <- round(sumreg$coef[, 5], 4)
p <- p %>%
  select(n, `exp(coef)`, `95% CI`, p)
p

Score_train$S_Score_All <- scale(Score_train$Score_All)
Score_train$S_Score_Clin <- scale(Score_train$Score_Clin)
Score_train$S_Score_Gene <- scale(Score_train$Score_Gene)
Score_train$S_Score <- scale(Score_train$Score)
coxph(Surv(PFI.time, PFI) ~ S_Score_All , data = Score_train)
coxph(Surv(PFI.time, PFI) ~ S_Score_Clin , data = Score_train)
coxph(Surv(PFI.time, PFI) ~ S_Score_Gene , data = Score_train)
coxph(Surv(PFI.time, PFI) ~ S_Score, data = Score_train)
coxph(Surv(PFI.time, PFI) ~ S_Score + S_Score_Gene + S_Score_Clin, data = Score_train)
reg <- coxph(Surv(PFI.time, PFI) ~ NCCN, data = Score_train)
summary(reg)

anova(cph_all, cph_G)
# anova(cph_all, cph_WC)
anova(cph_all, cph_W)
anova(cph_all, cph_C)
library(compareC)
compareC(Score_train$PFI.time,
         Score_train$PFI,
         Score_train$Score_All2,
         Score_train$Score)
compareC(Score_train$PFI.time,
         Score_train$PFI,
         Score_train$Score_All2,
         Score_train$Score_Gene)
compareC(Score_train$PFI.time,
         Score_train$PFI,
         Score_train$Score_All2,
         Score_train$Score_Clin)
compareC(Score_train$PFI.time,
         Score_train$PFI,
         Score_train$Score_All2,
         Score_train$Score_Clin2)
compareC(Score_train$PFI.time,
         Score_train$PFI,
         Score_train$Score_All2,
         Score_train$NCCN)

p_all <- time(Score_train, Score_train$Score_All)
p_Score <- time(Score_train, Score_train$Score)
p_Gene <- time(Score_train, Score_train$Score_Gene)
p_gleason <- time(Score_train, Score_train$Score_Clin)
p_all2 <- time(Score_train, Score_train$Score_All2)
p_gleason2 <- time(Score_train, Score_train$Score_Clin)
compare(p_all, p_Score)
compare(p_all, p_Gene)
compare(p_all, p_gleason)
mult(p_all, p_Score, p_Gene, p_gleason, 365 * 3)

#-----------------------------
Score_train %>%
  select(Score_All2, Score, Score_Gene, Score_Clin2) %>%
  write.csv('D:/科研/WSI-2024/GSEA/Score_TCGA.csv', row.names = F)
ggscatter(Score_train, x = "Score_Clin", y = "Score_Gene",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray")) +
  stat_cor(method = "pearson") +
  xlab('Cscore') + ylab('Gscore')
ggplot(Score_train, aes(x = Score_All2, y = Score_Gene)) +
  geom_point() +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_bw()

# cutoff包--------------------
library(cutoff)
# res3 %>% arrange(as.numeric(p.adjust))
cutpoint <- c()
name = ''
Score_train <- Score_train %>%
  filter(Score_All2 <= 5.396036)
for (i in 1:6){
  i = 2
  if (i == 1) name = 'Score_All'
  if (i == 2) name = 'Score_All2'
  if (i == 3) name = 'Score'
  if (i == 4) name = 'Score_Gene'
  if (i == 5) name = 'Score_Clin'
  if (i == 6) name = 'Score_Clin2'
  if (i >= 4){
    y = 0.05
  }else{
    y = 0
  }
  p <- cox(data=Score_train,
           time = 'PFI.time',#时间变量
           y='PFI', #生存状态
           x=name,
           cut.numb=1,
           n.per=0,#新分组x少数组最小比例
           y.per=0,
           round = 5)
  p <- p %>% 
    filter(as.numeric(pvalue) < 0.05) %>%
    arrange(desc(hr))
  cutpoint <- c(cutpoint, p[1, 'cut1'])
}
cutpoint <- data.frame(
  Score_All = cutpoint[1],
  Score_All2 = cutpoint[2],
  Score = cutpoint[3],
  Score_Gene = cutpoint[4],
  Score_Clin = cutpoint[5],
  Score_Clin2 = cutpoint[6]
)
write.csv(cutpoint, 'D:/科研/WSI-2024/PLCO/cutpoint.csv')

cutoff<-surv_cutpoint(Score_train, #数据集
                      time="PFI.time",#“ ”里写数据集时间变量
                      event="PFI",##“ ”里数据集结局变量名称
                      variables=c('Score_All2'))
summary(cutoff) #输出结果


#----------------------
cutpoint <- read.csv('D:/科研/WSI-2024/PLCO/cutpoint.csv', row.names = 1)
cutpoint$Score_All2 = 4.288867
Score_train2 <- Score_train %>%
  mutate(
    Score_All = ifelse(Score_All >= cutpoint$Score_All, 'High', 'Low'),
    Score_All2 = ifelse(Score_All2 >= cutpoint$Score_All2, 'High', 'Low'),
    Score = ifelse(Score >= cutpoint$Score, 'High', 'Low'),
    Score_Gene = ifelse(Score_Gene >= cutpoint$Score_Gene, 'High', 'Low'),
    Score_Clin = ifelse(Score_Clin >= cutpoint$Score_Clin, 'High', 'Low'),
    Score_Clin2 = ifelse(Score_Clin2 >= cutpoint$Score_Clin2, 'High', 'Low')
  ) %>%
  select(X_PATIENT, Score_All2, Score, Score_Gene, Score_Clin2)
write.csv(Score_train2, 'D:/科研/WSI-2024/GSEA/Score_train2.csv')

draw_KM <- function(Score_train, res.cut, name){
  Score_train[, name] <- ifelse(Score_train[, name] >= res.cut, 'High', 'Low')
  fit <- survfit(Surv(PFI.time, PFI) ~ Score_train[, name], data = Score_train)
  fit
  p <- ggsurvplot(fit, data = Score_train,
                  conf.int = TRUE, risk.table = TRUE,
                  legend.labs = c("High Score", "Low Score"))
  p$plot <- p$plot +
    ylab('Progression Free Interval Probability')
  data.survdiff <- survdiff(Surv(PFI.time, PFI) ~ Score_train[, name], data = Score_train)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")
  p$plot <- p$plot + 
    annotate("text", x = 0, y = 0.1, 
             label = ifelse(p.val < 0.0001, 
                            paste0(" p < 0.0001","\n HR (95% CI) = ",ci), 
                            paste0(" p value = ",sprintf("%.4f",p.val),"\n HR (95% CI) = ",ci)), 
             size = 5, color = "black", hjust = 0)+
    theme(text = element_text(size = 15))
  return (p)
}
pic <- list()
for (i in 1:ncol(cutpoint)){
  i = 2
  if (i == 1) name = 'Score_All'
  if (i == 2) name = 'Score_All2'
  if (i == 3) name = 'Score'
  if (i == 4) name = 'Score_Gene'
  if (i == 5) name = 'Score_Clin'
  if (i == 6) name = 'Score_Clin2'
  pic[[i]] <- draw_KM(Score_train, cutpoint[, i], name)
}
pic[[2]]
library(cowplot)
plot_grid(pic[[1]]$plot, pic[[2]]$plot, pic[[3]]$plot, 
          pic[[4]]$plot, pic[[5]]$plot, pic[[6]]$plot,
          labels = c("A", "B", "C",
                     'D', 'E', 'F'),
          ncol = 3, nrow = 2)

data_s2 <- Score_train %>%
  filter(NCCN == 3)
summary(coxph(Surv(PFI.time, PFI) ~ Score_All2 + Gleason,
              data = data_s2))
name = 'Score_All2'
res.cut <- 5.518809
data_s2[, name] <- ifelse(data_s2[, name] >= res.cut, 'High', 'Low')
fit <- survfit(Surv(PFI.time, PFI) ~ data_s2[, name], data = data_s2)
fit
p <- ggsurvplot(fit, data = data_s2,
                conf.int = TRUE, risk.table = TRUE,
                legend.labs = c("High Score", "Low Score"))
p$plot <- p$plot +
  ylab('Prostate Cancer Specific Survival')
data.survdiff <- survdiff(Surv(PFI.time, PFI) ~ data_s2[, name], data = data_s2)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")
p$plot <- p$plot + 
  annotate("text", x = 0, y = 0.1, 
           label = ifelse(p.val < 0.0001, 
                          paste0(" p < 0.0001","\n HR (95% CI) = ",ci), 
                          paste0(" p value = ",sprintf("%.4f",p.val),"\n HR (95% CI) = ",ci)),
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
p

# 绘制柱状图，要求每个样本一条柱子，按照从大到小排列，且填充颜色为PFI，柱子之间紧密排列
# 将Score_train2按照Score_All从大到小排序
Score_train2 <- Score_train %>%
  arrange(desc(Score_All2))
Score_train2$num <- seq(1, nrow(Score_train2), 1)
ggplot(Score_train2, aes(x = num, y = Score_All2)) +
  geom_bar(stat = "identity", aes(fill = as.character(PFI))) +
  labs(x = "Sample ID", y = "Score") +
  scale_fill_manual(values = c('0' = "grey", '1' = "red"),
                    labels = c('0' = "No Progress", '1' = "Progress"),
                    name = "Status") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = Score_train2[which.min(abs(rev(cummax(rev(Score_train2$Score_All2)))-cutpoint$Score_All2)),
                                       'num'], 
             colour = "black", size = 1, linetype = 'dashed') +
  theme_bw() + xlab('') + ylab('Risk Score') +
  theme(axis.text.x = element_blank())
