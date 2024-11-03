rm(list = ls())
gc()
library(reshape2)
library(paletteer)
library(corrplot)
library(cowplot)
setwd('/data/home/jiangminghui/KC/Project-0035/04_Correlation_Analysis/')
load('./rdata_04_Correlation_Analysis.RData')
Overlapgenes <-readRDS('../03_Overlap_genes_Enrichment/01_Overlapgenes.rds')
GSE77938 <-readRDS('../00_raw_data/combat_data.rds')
library(paletteer)

library(psych)

gene_data <- data.frame(GSE77938[Overlapgenes,] %>%t())

# 计算相关性矩阵和p值矩阵
cor_matrix <- cor(gene_data, method = "spearman")
p_matrix <- cor.mtest(gene_data, method = "spearman")$p

# 设置阈值
threshold <- 0.6

# 过滤相关性矩阵
cor_matrix[abs(cor_matrix) <= threshold | p_matrix >= 0.05] <- 0

# 颜色设置
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10) # 整成10个色阶

pdf('./Heatmap.pdf',width = 6,height = 6)
# 绘制相关性圆圈图
corrplot(cor_matrix, type = "upper", 
         method = "pie",
         order = "hclust", 
         col = my_color,
         tl.col = "black", 
         tl.srt = 45)
dev.off()

png('./Heatmap.png',width = 6,height = 6,units = 'in',res = 300)
# 绘制相关性圆圈图
corrplot(cor_matrix, type = "upper", 
         method = "pie",
         order = "hclust", 
         col = my_color,
         tl.col = "black", 
         tl.srt = 45)
dev.off()

# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'rdata_',dir_name, ".RData"))