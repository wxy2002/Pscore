##--------------------------------------------------
library(tidyverse)

clin <- read.csv('./package-plcoi-1792.2025-01-16/Prostate/pros_data_t20241011.csv')
path_data <- read.csv('./package-plcoi-1792.2025-01-16/Prostate/Pathology\ Image\ Linkage/pros_image_data_t20241011.csv')
clin <- clin %>%
  filter(pros_cancer == 1) %>%
  filter(pros_has_deliv_heslide_img == 1) %>%
  subset(!is.na(pros_cancer_diagdays)) %>%
  mutate(OS.time = dth_days - pros_cancer_diagdays) %>%
  mutate(OS = is_dead) %>%
  mutate(time_data = mortality_exitdays) %>%
  subset(!is.na(pros_stage_t)) %>%
  subset(!is.na(pros_dx_psa)) %>%
  subset(!is.na(age))
clin$OS.time <- ifelse(clin$OS == 0, clin$time_data, clin$OS.time)
clin$PCSM <- ifelse(clin$d_cause_of_death == 1, 1, 0)
clin$PCSM1 <- clin$PCSM
clin$PCSM2 <- clin$PCSM
clin$PCSM1[is.na(clin$PCSM)] <- 0
clin$PCSM2[is.na(clin$PCSM)] <- 2
path_data <- path_data %>%
  filter(plco_id %in% clin$plco_id)
clin$gleason = clin$pros_gleason

##--------------------------------------------------
Score <- read.csv('Score_PLCO_wsi.csv')
Score$X_PATIENT <- paste0(Score$X_PATIENT, '.svs')
path_s <- path_data %>%
  filter(pros_heslide_img_name %in% Score$X_PATIENT)
data_s <- clin %>%
  filter(plco_id %in% path_s$plco_id)
load('model.rdata')
colnames(Score) <- c('pros_heslide_img_name', 'Score')
Score <- left_join(Score, path_s, by = 'pros_heslide_img_name')
p <- Score %>%
  group_by(plco_id) %>%
  summarize(Max = max(Score, na.rm = T))
colnames(p) <- c('plco_id', 'Score')
data_s <- left_join(data_s, p, by = 'plco_id')

pre <- read.csv('PLCO_Gene.csv')
pre$Score_Gene <- 2.09636*pre$CHTF18 - 
  0.11524*pre$CDCA5 - 1.00421*pre$TRPM4 + 
  0.02448*pre$SERPINB5 + 0.51506*pre$ARHGEF40 - 
  0.87403*pre$LIX1
pre <- pre %>%
  select(PATIENT, Score_Gene)
for (i in 1:length(pre[, 1])){
  pre[i, 'pros_heslide_img_name'] = paste(pre[i, 1], '.svs', sep = '')
}
pre <- left_join(pre, path_s, by = 'pros_heslide_img_name')
p <- pre %>%
  group_by(plco_id) %>%
  summarize(Max = max(Score_Gene))
colnames(p) <- c('plco_id', 'Score_Gene')
data_s <- left_join(data_s, p, by = 'plco_id')

library(survival)
data_s$Gleason <- data_s$gleason
data_s$Gleason <- ifelse(data_s$Gleason <= 6, 6, data_s$Gleason)
data_s <- data_s %>% 
  mutate(cT = case_when(
    pros_stage_t < '200' ~ 1,
    pros_stage_t < '300' ~ 2,
    pros_stage_t < '400' ~ 3,
    pros_stage_t >= '400' ~ 4,
    TRUE ~ 0
  )) %>%
  mutate(NCCN = case_when(
    pros_stage_t <= '210' & Gleason <= 6 & pros_dx_psa < 10 ~ 1,
    
    pros_stage_t >= '310' | Gleason >= 8 | (pros_dx_psa > 20) ~ 3,
    TRUE ~ 2
  ))
data_s$Score_Clin2 <- 0.7759*data_s$Gleason + 0.4513*data_s$cT

data_s$Score_All <- 0.4332 * data_s$Score_Gene + 
  1.1055 * data_s$Score + 0.7224 * data_s$Score_Clin2
data_s$Score_All_Scale <- scale(data_s$Score_All)
data_s <- data_s %>% 
  mutate(cT_2 = case_when(
    cT >= 2 ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(Gleason_2 = case_when(
    Gleason >= 7 ~ 1,
    TRUE ~ 0
  ))
reg <- coxph(Surv(OS.time, PCSM1) ~ Score_Clin2, data = data_s)
summary(reg)
reg <- coxph(Surv(OS.time, PCSM1) ~ Score_Gene, data = data_s)
summary(reg)
reg <- coxph(Surv(OS.time, PCSM1) ~ Score, data = data_s)
summary(reg)
reg <- coxph(Surv(OS.time, PCSM1) ~ Score_All, data = data_s)
summary(reg)
reg <- coxph(Surv(OS.time, PCSM1) ~ NCCN, data = data_s)
summary(reg)
cutpoint <- read.csv('cutpoint.csv', header = T, row.names = 1)

data_s %>%
  select(Score_All, Score, Score_Gene, Score_Clin2) %>%
  write.csv('D:/科研/WSI-2024/GSEA/Score_PLCO.csv', row.names = F)
Score_train <- read.csv('D:/科研/WSI-2024/GSEA/Score_TCGA.csv')
Score_test <- read.csv('D:/科研/WSI-2024/GSEA/Score_PLCO.csv')
colnames(Score_train) <- c('Pscore', 'Wscore', 'Gscore', 'Cscore')
colnames(Score_test) <- c('Pscore', 'Wscore', 'Gscore', 'Cscore')
Score_train$label <- 'Train'
Score_test$label <- 'Test'
Score <- rbind(Score_train, Score_test)
library(ggpubr)
ggscatter(Score, x = "Gscore", y = "Cscore", alpha = 0.5,
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "red")) +
  stat_cor(method = "spearman")
Score$id <- seq(1, nrow(Score), 1)
# 根据Score的值进行排序
Score <- Score %>%
  arrange(desc(Pscore))
l <- Score$id
Score_long <- Score %>%
  pivot_longer(cols = c('Pscore', 'Wscore', 'Gscore', 'Cscore'),
               names_to = 'Type', values_to = 'Score') %>%
  mutate(Type = factor(Type, levels = c('Pscore', 'Wscore', 'Gscore', 'Cscore')))
Score_long$id <- factor(Score_long$id, levels=l, ordered=T)
ggplot(Score_long, aes(x = id, y = Type)) +
  geom_raster(aes(fill = Score))  +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Type", fill = "Values") +
  ggtitle("Heatmap of Scores") + scale_fill_gradientn(colours=c("blue","green","yellow","red"))

data_s$Score_All2 <- ifelse(data_s$Score_All >= cutpoint$Score_All2, 1, 0)
data_s$psa <- log2(data_s$pros_dx_psa + 1)
reg <- coxph(Surv(OS.time, PCSM1) ~ Score_All + age + psa + Gleason_2 + cT,
             data = data_s)
sumreg <- summary(reg)
p <- round(sumreg$coef, 4) %>% as.data.frame()
p$`95% CI` <- paste0(round(sumreg$conf.int[, 3], 4), '-', round(sumreg$conf.int[, 4], 4))
p$n <- round(sumreg$n, 4)
p$p <- round(sumreg$coef[, 5], 4)
p <- p %>%
  select(n, `exp(coef)`, `95% CI`, p)
# 针对上述多变量回归中的相关临床变量，进行批量单变量cox回归
# 将多个变量的n/HR/95%ci/p整合到一个dataframe中
clin_var <- c('age', 'pros_dx_psa', 'Gleason', 'cT', 'Score_All2')
clin_var <- data.frame(clin_var)
clin_var <- clin_var %>%
  mutate(n = NA, HR = NA, `95% CI` = NA, p = NA)
for (i in 1:length(clin_var[, 1])){
  if (clin_var[i, 'clin_var'] == 'age'){
    reg <- coxph(Surv(OS.time, PCSM1) ~ age, data = data_s)
  } else if (clin_var[i, 'clin_var'] == 'pros_dx_psa'){
    reg <- coxph(Surv(OS.time, PCSM1) ~ log2(pros_dx_psa + 1), data = data_s)
  } else if (clin_var[i, 'clin_var'] == 'Gleason'){
    reg <- coxph(Surv(OS.time, PCSM1) ~ Gleason_2, data = data_s)
  } else if (clin_var[i, 'clin_var'] == 'cT'){
    reg <- coxph(Surv(OS.time, PCSM1) ~ cT, data = data_s)
  } else if (clin_var[i, 'clin_var'] == 'Score_All2'){
    reg <- coxph(Surv(OS.time, PCSM1) ~ Score_All, data = data_s)
  }
  sumreg <- summary(reg)
  clin_var[i, 'n'] <- round(sumreg$n, 4)
  clin_var[i, 'HR'] <- round(sumreg$coef[2], 4)
  clin_var[i, '95% CI'] <- paste0(round(sumreg$conf.int[3], 4), '-', round(sumreg$conf.int[4], 4))
  clin_var[i, 'p'] <- round(sumreg$sctest[3], 4)
}
clin_var

library(compareGroups)
data_s2 <- data_s %>% mutate(
  Gleason = as.character(Gleason),
  cT = case_when(
    cT == 1 ~ 'T1',
    cT == 2 ~ 'T2',
    cT == 3 ~ 'T3',
    cT == 4 ~ 'T4',
    TRUE ~ 'Unknown'
  ),
  OS = as.character(OS),
  PCSM = as.character(PCSM1),
  NCCN = case_when(
    NCCN == 1 ~ 'Low',
    NCCN == 2 ~ 'Intermediate',
    NCCN == 3 ~ 'High',
    TRUE ~ 'Unknown'
  )
)
p <- descrTable(~ age + pros_dx_psa + Gleason + cT + OS + PCSM + OS.time + NCCN, 
                data = data_s2,
                method = c(age = 2, pros_dx_psa = 2,
                           OS.time = 2))
export2xls(p, 'PLCO_Clin_Table.xlsx')

library(compareC)
# coxph(Surv(OS.time, PCSM) ~ Score_Clin2, data = data_s)
compareC(data_s$OS.time,
         data_s$PCSM1,
         data_s$Score_All,
         data_s$Score_Clin2)
compareC(data_s$OS.time,
         data_s$PCSM1,
         data_s$Score_All,
         data_s$Score_Gene)
compareC(data_s$OS.time,
         data_s$PCSM1,
         data_s$Score_All,
         data_s$Score)
compareC(data_s$OS.time,
         data_s$PCSM1,
         data_s$Score_All,
         data_s$NCCN)

#------------------------------------------------------
time <- function(data, name){
  library(timeROC)
  # p <- predict(cph, data)
  ROC <- timeROC(T = data$OS.time, #生存时间
                 delta = data$PCSM1,   #生存状态
                 marker = data[, name], #计算timeROC的变量
                 cause = 1,
                 weighting = "marginal",
                 times = c(5 * 365, 10 * 365, 15 * 365, 20 * 365,
                           25 * 365),
                 iid=TRUE)
  return(ROC)
}
auc_all <- time(data_s, 'Score_All')
auc_wsi <- time(data_s, 'Score')
auc_gene <- time(data_s, 'Score_Gene')
auc_clin <- time(data_s, 'Score_Clin2')
compare(auc_all, auc_clin)
compare(auc_all, auc_wsi)
compare(auc_all, auc_gene)
mult <- function(p_all, p_WSI, p_Gene, p_Clin, time){
  library(ggplot2)
  label <- case_when(
    time == 5 * 365 ~ 1,
    time == 10 * 365 ~ 2,
    time == 15 * 365 ~ 3,
    time == 20 * 365 ~ 4,
    time == 25 * 365 ~ 4
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
mult(auc_all, auc_wsi, auc_gene, auc_clin, 25 * 365)

#-----------------------------------
library(survminer)
draw_KM <- function(data_s, res.cut, name){
  # print(median(Score_train[, name]))
  data_s[, name] <- ifelse(data_s[, name] >= res.cut, 'High', 'Low')
  fit <- survfit(Surv(OS.time, PCSM1) ~ data_s[, name], data = data_s)
  fit
  p <- ggsurvplot(fit, data = data_s,
                  conf.int = TRUE, risk.table = TRUE,
                  legend.labs = c("High Score", "Low Score"))
  p$plot <- p$plot +
    ylab('Prostate Cancer Specific Survival')
  data.survdiff <- survdiff(Surv(OS.time, PCSM1) ~ data_s[, name], data = data_s)
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
cutpoint <- read.csv('cutpoint.csv', header = T, row.names = 1)
cutpoint$Score_All2 = 4.288867
data_s <- data_s %>%
  filter(Score_All <= 5.396036)
pic <- list()
for (i in 1:ncol(cutpoint)){
  if (i == 1) next
  if (i == 2) name = 'Score_All'
  if (i == 3) name = 'Score'
  if (i == 4) name = 'Score_Gene'
  if (i == 5) next
  if (i == 6) name = 'Score_Clin2'
  pic[[i]] <- draw_KM(data_s, cutpoint[, i], name)
}
library(cowplot)
plot_grid(pic[[2]]$plot, pic[[3]]$plot, pic[[4]]$plot, pic[[6]]$plot,
          labels = c("A", "B", "C",
                     'D', 'E', 'F'),
          ncol = 2, nrow = 2)

data_s2 <- data_s %>%
  filter(cT > 2)
summary(coxph(Surv(OS.time, PCSM1) ~ Score_All + Gleason,
               data = data_s2))
name = 'Score_All'
res.cut <- cutpoint$Score_All2
data_s2[, name] <- ifelse(data_s2[, name] >= res.cut, 'High', 'Low')
fit <- survfit(Surv(OS.time, PCSM1) ~ data_s2[, name], data = data_s2)
fit
p <- ggsurvplot(fit, data = data_s2,
                conf.int = TRUE, risk.table = TRUE,
                legend.labs = c("High Score", "Low Score"))
p$plot <- p$plot +
  ylab('Prostate Cancer Specific Survival')
data.survdiff <- survdiff(Surv(OS.time, PCSM1) ~ data_s2[, name], data = data_s2)
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

#-----------------------------------
data_s2 <- data_s %>% filter(NCCN == 3) %>%
  arrange(desc(Score_All))
data_s2$num <- seq(1, nrow(data_s2), 1)
ggplot(data_s2, aes(x = num, y = Score_All)) +
  geom_bar(stat = "identity", aes(fill = as.character(PCSM1))) +
  scale_fill_manual(values = c('0' = "grey", '1' = "red"),
                    labels = c('0' = "Not Die for prostate", '1' = "Die for prostate"),
                    name = "Status") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_vline(xintercept = data_s2[which.min(abs(rev(cummax(rev(data_s2$Score_All)))-cutpoint$Score_All2)), 'num'], 
             colour = "black", size = 1, linetype = 'dashed') +
  theme_bw() + xlab('') + ylab('Risk Score') +
  theme(axis.text.x = element_blank())

#-----------------------------------
data_s2 <- data_s %>%
  select(Score_All, OS.time, PCSM1, age, Gleason, cT, NCCN)
data_s2$Score_All <- ifelse(data_s2$Score_All >= cutpoint$Score_All2, 
                            'High risk', 'Low risk')
data_s2$age <- ifelse(data_s2$age >= 60, 'Age >= 60', 'Age < 60')
data_s2$NCCN <- ifelse(data_s2$NCCN == '1', 'Low risk', 
                       ifelse(data_s2$NCCN == '2', 'Intermediate risk', 'High risk'))
data_s2$Gleason <- ifelse(data_s2$Gleason >= 8, 'Gleason >= 8', 'Gleason <= 7')
data_s2 <- data_s2 %>%
  mutate(cT = case_when(
    cT <= 2 ~ 'T1-T2',
    cT >= 3 ~ 'T3-T4',
    TRUE ~ 'Unknown'
  )) %>%
  mutate(
    Score_All = factor(Score_All, levels = c('Low risk', 'High risk')),
  )
dfl <- data_s2 %>% 
  pivot_longer(cols = 4:ncol(.), names_to = "var",values_to = "value") %>% 
  arrange(var)
head(dfl)
ress <- dfl %>% 
  #group_by(var,value) %>% 
  group_nest(var,value) %>% 
  drop_na(value) %>% 
  mutate(
    model=map(data, ~ coxph(Surv(OS.time, PCSM1) ~ Score_All, data = .x)),
    res = map(model, broom::tidy, conf.int = T, exponentiate = T)
  ) %>% 
  dplyr::select(var,value,res)
ss <- dfl %>% 
  group_by(var, value, Score_All) %>% 
  drop_na(value) %>% 
  summarise(sample_size=n()) %>% 
  dplyr::select(var, value, Score_All, sample_size)
resss <- ress %>% 
  left_join(ss, b=c("var","value")) %>% 
  unnest(res, Score_All, sample_size) %>% 
  pivot_wider(names_from = "Score_All", values_from = "sample_size", names_prefix = "Pscore_") %>% 
  select(-c(term, std.error, statistic)) %>% 
  mutate(across(where(is.numeric), round,digits = 2)) %>% 
  mutate(`HR(95%CI)`=paste0(estimate,"(",conf.low,"-",conf.high,")"))

str(resss)
head(resss)

