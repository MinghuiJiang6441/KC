rm(list =ls())
gc()
#.libPaths("/data/nas2/Software/miniconda3/envs/public_R/lib/R/library"  )
library(ggplot2)
library(reshape2)
library(RColorBrewer)
setwd("/data/home/jiangminghui/KC/Project-0035/11_Immunoinfiltration_Analysis/")
load('./rdata_11_Immunoinfiltration_Analysis.RData')
library(preprocessCore)

library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
#BiocManager::install("RcolorBrewer")
#install.packages("RColorBrewer")
library(RColorBrewer)
library(tidyHeatmap)
library(tidyverse)
library(RColorBrewer)
library(tidyr)
library(tibble)
library(ggsci)
library(IOBR)
library(psych)
library(reshape2)
library(ggsci)
library(psych)
library(Heatmap)
library(corrplot)
library(cowplot)　
library(paletteer)

 source('./CIBERSORT.R')
GSE77938 <-readRDS('../00_raw_data/combat_data.rds')

LM22.file<- read.table(file = "./LM22.txt", header=T, sep="\t", row.names=1,check.names=F)
GSE77938 <-as.data.frame(GSE77938)

GEO_cibersort.results <- CIBERSORT(LM22.file, GSE77938, perm = 50, QN = TRUE, absolute =F)

saveRDS(GEO_cibersort.results,'./01_GEO_cibersort.results.rds')

GEO_cibersort.results <-readRDS('./01_GEO_cibersort.results.rds')

GEO_cibersort.results <- as.data.frame(GEO_cibersort.results)

group <- ifelse(grepl("^KC", rownames(GEO_cibersort.results)), "Disease", "Control")

# GEO_cibersort.results <- GEO_cibersort.results %>% filter(`P-value` <= 0.05)
dim(GEO_cibersort.results)

cibersort_data <- as.data.frame(GEO_cibersort.results[,1:22])


diff_cells <- data.frame(Cell_Type = character(), p_value = numeric(), stringsAsFactors = FALSE)

# 对每种免疫细胞类型进行 Wilcoxon 检验
for (cell_type in colnames(cibersort_data)) {
  # 提取当前免疫细胞类型的数据
  cell_data <- cibersort_data[[cell_type]]
  
  # 进行 Wilcoxon 检验
  test_result <- wilcox.test(cell_data ~ group)
  
  if (!is.na(test_result$p.value) && test_result$p.value < 0.05) {
    diff_cells <- rbind(diff_cells, data.frame(Cell_Type = cell_type, p_value = test_result$p.value))
  }
}
diff_cells

diff_data <-cibersort_data[,diff_cells$Cell_Type]
diff_data_numeric <- as.data.frame(lapply(diff_data, as.numeric))

# 重新计算每个样本的细胞类型比例，使其合计为1
diff_data_normalized <- as.data.frame(t(apply(diff_data_numeric, 1, function(x) x / sum(x, na.rm = TRUE))))


diff_data_normalized$Sample <- rownames(diff_data)
diff_data_normalized$Group <- group
rownames(diff_data_normalized) <- NULL

# 将数据转换为长格式
long_data <- melt(diff_data_normalized, id.vars = c("Sample", "Group"))
colors <-( brewer.pal(n = length(unique(long_data$variable)), name = "Set3"))

p1 <-ggplot(long_data, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = colors) +
  labs(x = "Sample", y = "Proportion", fill = "Immune Cell Type") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), # 透明背景
    plot.background = element_rect(fill = "transparent", color = NA),  # 透明背景
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 14) 
  )

ggsave('./01_Training_set_immune_cell_infiltration.pdf',p1,width = 10,height = 4)

ggsave('./01_Training_set_immune_cell_infiltration.png',p1,width = 10,height = 4)




diff_data_numeric <-cibersort_data[,diff_cells$Cell_Type]
diff_data_numeric$Group <- group
diff_data_long <- melt(diff_data_numeric, id.vars = "Group", variable.name = "Cell_Type", value.name = "Infiltration_Rate")

# 标准化处理
diff_data_long$Infiltration_Rate <- log2(diff_data_long$Infiltration_Rate  +1)

# 绘制箱线图
p2 <-ggplot(diff_data_long, aes(x = Cell_Type, y = Infiltration_Rate, fill = Group)) +
  geom_boxplot() +
  labs(x = "Immune Cell Type", y = "Infiltration Rate (log2 + 1)", title = "Boxplot of Infiltration Rate by Group (Standardized)") +
  scale_fill_manual(values = c('#0339A6', '#A60A27')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = list(c("Disease", "Control")), 
              map_signif_level = TRUE, 
              test = "wilcox.test")

ggsave('./02_Training_set_immune_cell_infiltration.pdf',p2,width = 6,height = 6)

ggsave('./02_Training_set_immune_cell_infiltration.png',p2,width = 6,height = 6)

# 提取需要的数据
diff_data <- cibersort_data[, diff_cells$Cell_Type]
diff_data_numeric <- as.data.frame(lapply(diff_data, as.numeric))

cor_matrix <- cor(diff_data_numeric, method = "spearman")
p_matrix <- cor.mtest(diff_data_numeric, method = "spearman")$p

# 设置阈值
threshold <- 0.3

# 过滤相关性矩阵
cor_matrix[abs(cor_matrix) <= threshold | p_matrix >= 0.05] <- 0

# 颜色设置
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10) # 整成10个色阶

pdf('./03_Differential_immune_cell_correlation_heatmap.pdf',width =6,height = 5 )
corrplot(cor_matrix, type = "upper", 
         method = "pie",
         order = "hclust", 
         col = my_color,
         tl.col = "black", 
         tl.srt = 45)
dev.off()

png('./03_Differential_immune_cell_correlation_heatmap.png',width =6,height = 5 ,res = 300,units = 'in')
corrplot(cor_matrix, type = "upper", 
         method = "pie",
         order = "hclust", 
         col = my_color,
         tl.col = "black", 
         tl.srt = 45)
dev.off()
# 
# # 进行 Spearman 相关性分析
# cor_matrix <- corr.test(diff_data_numeric, method = "spearman")
# 
# # 提取相关性系数和 p 值
# cor_r <- cor_matrix$r
# cor_p <- cor_matrix$p
# 
# cor_r_filtered <- cor_r
# cor_r_filtered[abs(cor_r) <= 0.3 | cor_p >= 0.05] <- 0
# library(pheatmap)
# 
# # 假设cor_r_filtered是你的相关性矩阵
# # 创建颜色渐变
# color_palette <- colorRampPalette(c("#04198C", "white", "#D91E1E"))(50)
# 
# # 定义断点，使得0值显示为白色
# breaks <- seq(-1, 1, length.out = 51)
# # mid_point <- length(breaks) / 2
# # breaks[(mid_point - 1):(mid_point + 1)] <- c(breaks[mid_point - 1], 0, breaks[mid_point + 1])
# # sig_matrix <- ifelse(cor_p < 0.05, 1, NA)  # 显著的设置为1，不显著的设置为NA
# 
# # 创建显著性注释
# annotation_sig <- ifelse(sig_matrix == 1, "black", "transparent")
# # 创建热图
# pdf('./03_Differential_immune_cell_correlation_heatmap.pdf',width =6,height = 5 )
# pheatmap(cor_r_filtered, 
#          main = "Spearman Correlation Heatmap  \nof Differential Immune Cells", 
#          xlab = "Immune Cells", 
#          ylab = "Immune Cells", 
#          color = color_palette,
#          breaks = breaks)
# 
# dev.off()

# png('./03_Differential_immune_cell_correlation_heatmap.png',width =6,height = 5 ,res = 300,units = 'in')
# pheatmap(cor_r_filtered, 
#          main = "Spearman Correlation Heatmap  \nof Differential Immune Cells", 
#          xlab = "Immune Cells", 
#          ylab = "Immune Cells", 
#          color = color_palette,
#          breaks = breaks)
# 
# dev.off()

Biom <-readRDS('../07_ValidateExpressions/01_Candidate_biomarker.rds')

GSE77938<-as.matrix(GSE77938) 
gene_mat <-GSE77938[Biom,] %>%t()

cor_matrix <- corr.test(diff_data,gene_mat, method = "spearman")
# 提取相关性系数和 p 值
cor_r <- cor_matrix$r
cor_p <- cor_matrix$p

cor_r_filtered <- cor_r
cor_r_filtered[abs(cor_r) <= 0.3 | cor_p >= 0.05] <- 0
# 创建颜色渐变
color_palette <- colorRampPalette(c("#04198C", "white", "#D91E1E"))(50)

# 定义断点，使得0值显示为白色
breaks <- seq(-1, 1, length.out = 51)
# mid_point <- length(breaks) / 2
# breaks[(mid_point - 1):(mid_point + 1)] <- c(breaks[mid_point - 1], 0, breaks[mid_point + 1])
# sig_matrix <- ifelse(cor_p < 0.05, 1, NA)  # 显著的设置为1，不显著的设置为NA
significance_labels <- matrix(ifelse(cor_matrix$p < 0.0001, "****",
                                     ifelse(cor_matrix$p < 0.001, "***",
                                            ifelse(cor_matrix$p < 0.01, "**",
                                                   ifelse(cor_matrix$p < 0.05, "*", " ")))),
                              nrow = nrow(cor_matrix$p))


significance_labels[abs(cor_r) <= 0.3 ] <- " "
# 创建热图
pdf('./04_Heatmap_correlation_between_biomarkers_differential_immunecells.pdf',width =5,height = 5 )
pheatmap(cor_r_filtered, 
         main = "Heatmap of correlation between \nbiomarkers and differential immune cells", 
         # xlab = "Immune Cells", 
         # ylab = "Immune Cells", 
         color = color_palette,
         display_numbers = significance_labels, 
         breaks = breaks)

dev.off()

png('./04_Heatmap_correlation_between_biomarkers_differential_immunecells.png',width =5,height = 5,units = 'in',res = 300 )
pheatmap(cor_r_filtered, 
         main = "Heatmap of correlation between \nbiomarkers and differential immune cells", 
         # xlab = "Immune Cells", 
         # ylab = "Immune Cells", 
         color = color_palette,
         display_numbers = significance_labels, 
         breaks = breaks)

dev.off()


# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'rdata_',dir_name, ".RData"))