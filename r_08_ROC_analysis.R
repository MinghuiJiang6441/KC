rm(list = ls())
gc()
setwd('/data/home/jiangminghui/KC/Project-0035/08_ROC_analysis/')

load('./Rdata_08_ROC_analysis.RData')
library(dplyr)
library(tidyr)
library(ggpubr)
library(pROC)

GSE77938_exp <-readRDS('../00_raw_data/combat_data.rds')

GSE151631_exp <-readRDS('../00_raw_data/GSE151631_exp.rds')

Candidate_biomarker <-readRDS('../07_ValidateExpressions/01_Candidate_biomarker.rds')

train_data<-list()

train_data$Group<- ifelse(grepl("^KC", colnames(GSE77938_exp)), 1,0)

plot_combined_roc <- function(data, biomarkers, exp_data) {
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity", main = " ")
  colors <- rainbow(length(biomarkers))
  
  for (i in seq_along(biomarkers)) {
    biomarker <- biomarkers[i]
    data$Biomarker_exp <- exp_data[biomarker, ] %>% as.numeric()
    
    roc_result <- roc(data$Group, data$Biomarker_exp)
    auc_value <- auc(roc_result)
    
    plot(roc_result, col = colors[i], add = TRUE)
    text(0.4, 0.4 - 0.05 * i, paste(biomarker, "AUC =", round(auc_value, 2)), col = colors[i])
  }
  
  legend("bottomleft", legend = biomarkers, col = colors, lwd = 2)
}

# 调用函数
pdf('./01_train_ROC_curve.pdf',width = 7,height = 7)
plot_combined_roc(train_data, Candidate_biomarker, GSE77938_exp)
dev.off()

png('./01_train_ROC_curve.png',width = 7,height = 7,units = 'in',res=300)
plot_combined_roc(train_data, Candidate_biomarker, GSE77938_exp)
dev.off()


valid_data <-list()
valid_data$Group <-ifelse(grepl("^DN|^LE", colnames(GSE151631_exp)), 0, 1)
valid_data$Biomarker_exp <- GSE151631_exp[Candidate_biomarker[2], ] %>% as.numeric()


pdf('./02_valid_ROC_curve.pdf',width = 7,height = 7)
plot_combined_roc(valid_data, Candidate_biomarker, GSE151631_exp)
dev.off()

png('./02_valid_ROC_curve.png',width = 7,height = 7,units = 'in',res=300)
plot_combined_roc(valid_data, Candidate_biomarker, GSE151631_exp)
dev.off()

saveRDS(train_data,'./01_train_data.rds')

saveRDS(valid_data,'./01_valid_data.rds')




# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'Rdata_',dir_name, ".RData"))