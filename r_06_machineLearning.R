rm(list = ls())
gc()

setwd('/data/home/jiangminghui/KC/Project-0035/06_machineLearning/')
#load('./rdata_06_machineLearning.RData')
library(e1071)
library(caret)
library(ggplot2)
library(cowplot)
library(glmnet)
library(ggvenn)
library(Boruta)
GSE77938 <-readRDS('../00_raw_data/combat_data.rds')
Overlapgenes <-readRDS('../03_Overlap_genes_Enrichment/01_Overlapgenes.rds')


gene_data <- data.frame(GSE77938[Overlapgenes,] %>% t())
gene_data$labels <- ifelse(grepl("^KC", rownames(gene_data)), "Disease", "Control")
gene_data$labels <- factor(gene_data$labels, levels = c('Disease', 'Control'))

set.seed(21) # 设置种子
control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)

num <- ncol(gene_data)-1
results <- rfe(x = gene_data[, 1:num], # 除去最后一列，其余列均为预测变量（也就是hubgene的表达量）
               y = gene_data$labels, # 分组信息
               sizes = c(1:num), 
               rfeControl = control,
               method = "svmRadial"
)

# 执行SVM-RFE分析
saveRDS(results,'./01_svmProfile.rds')
results<-readRDS('./01_svmProfile.rds')
svmrfe_result <- data.frame(symbol = predictors(results)) 
#NUAK1, GATA6, CCDC13, RBM24, CARD16
p1 <- plot(results, type = c("o"), xgap.axis = 1,lwd =6)

# 使用plot_grid函数组合图形
combined_plot <- plot_grid(p1, labels = "AUTO")

# 添加标题和主题元素
combined_plot <- combined_plot + 
  draw_label("SVM_RFE_analyse", x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25, fontface = "bold", color = "black") +
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 25),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = 8, y = 0.92, label = "n=8, Accuracy = 0.92", size = 5, hjust = 0.75, vjust = 0.75)

# 显示组合图形
print(combined_plot)

ggsave('./01_SVM_RFE_analyse.pdf',combined_plot,width = 7,height = 7)
ggsave('./01_SVM_RFE_analyse.png',combined_plot,width = 7,height = 7)
labels <-ifelse(grepl("^KC", rownames(gene_data)), "Disease", "Control")
# x <- as.matrix(gene_data[,-ncol(gene_data)])
# y <- ifelse(labels =='Disease',1,0)
# 
# num=111
# set.seed(num)
# fit <-glmnet(x,y,alpha = 1)
# 
# pdf('./02_Lassoplot_1.pdf',width = 6,height = 5)
# plot(fit,xvar="lambda",label=T)
# dev.off()
# 
# png('./02_Lassoplot_1.png',width = 6,height = 5,units = 'in',res = 300)
# plot(fit,xvar="lambda",label=T)
# dev.off()
# 
# # 执行LASSO回归分析
# set.seed(num)#不要改！！！！！！！
# cv_fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
# 
# 
# # 获取与最小偏似然偏差相对应的惩罚参数值（λ）
# best_lambda <- cv_fit$lambda.min
# set.seed(num)#不要改！！！！！！！
# 
# # 使用最佳λ值重新拟合模型
# pdf('./02_Lassoplot.pdf',width = 6,height = 5)
# plot(cv_fit, xvar = "lambda", label = TRUE)
# abline(v = log(best_lambda), col = " red", lty = 2) 
# dev.off()
# 
# png('./02_Lassoplot.png',width = 6,height = 5,units = 'in',res = 300)
# plot(cv_fit, xvar = "lambda", label = TRUE)
# abline(v = log(best_lambda), col = " red", lty = 2) 
# dev.off()
# 
# lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
# 
# # 获取β系数矩阵
# coef_matrix <- as.matrix(coef(lasso_model))
# 
# # 获取β系数不为0的特征基因
# selected_genes <- rownames(coef_matrix)[coef_matrix != 0]
# selected_genes <- selected_genes[selected_genes != "(Intercept)"]
# 
# # 打印选中的特征基因
# print(selected_genes)
# 



set.seed(1)

boruta <- Boruta(x=gene_data[,-ncol(gene_data)], y=gene_data$labels, pValue=0.01, mcAdj=T, 
                 maxRuns=300)

boruta
table(boruta$finalDecision)

boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(boruta)

gene_con <- subset(boruta.variable.imp,boruta.variable.imp$finalDecision == 'Confirmed')
boruta_signif <-unique(gene_con$Variable) %>% as.character()


feature_gene <-intersect(svmrfe_result$symbol,boruta_signif)
feature_gene

saveRDS(feature_gene,'./03_feature_gene.rds')


write.csv(feature_gene,file = './03_feature_gene.csv')

feature_gene <-readRDS('./03_feature_gene.rds')

venn <- list('SVM-RFE' =svmrfe_result$symbol,
             'boruta' = boruta_signif
)

ggvenn(venn)    

pdf('./03_Venn_Plot.pdf',width = 6,height = 6)
ggvenn(venn)    
dev.off()

png('./03_Venn_Plot.png',width = 6,height = 6,units = 'in',res = 300)
ggvenn(venn)    
dev.off()



# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'rdata_',dir_name, ".RData"))