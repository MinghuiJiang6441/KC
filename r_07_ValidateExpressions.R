rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggpubr)

setwd('/data/home/jiangminghui/KC/Project-0035/07_ValidateExpressions/')

load('./Rdata_07_ValidateExpressions.RData')

GSE77938_exp <-readRDS('../00_raw_data/combat_data.rds')

GSE151631_exp <-readRDS('../00_raw_data/GSE151631_exp.rds')

feature_gene <-readRDS('../06_machineLearning/03_feature_gene.rds')
# 
# Candidate_biomarker <- setdiff(feature_gene, c('CARD16', 'CCDC13'))
# Candidate_biomarker <-feature_gene
# #[1] "GATA6" "RBM24" "STC2"  "HBEGF"

# 提取Overlapgenes的表达数据
expr_data <- GSE77938_exp[feature_gene, ]

# 转置表达数据以便每行代表一个样本
expr_data <- t(expr_data)

# 将表达数据转换为数据框
expr_data <- as.data.frame(expr_data)

sample_names <-rownames(expr_data)

# 创建分组向量
group <- ifelse(grepl("^KC", sample_names), "Disease", "Control")

# 添加样本信息
expr_data$SampleID <- rownames(expr_data)
expr_data$Group <- group

# 计算每个基因的Wilcoxon检验p值
p_values <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group) %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Expression ~ Group)$p.value)

# 查看结果
print(p_values)
# 合并表达数据和p值
expr_data_long <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group)

expr_data_long <- merge(expr_data_long, p_values, by = "Gene")

# 绘制基因表达柱状图并添加检验
p1 <- ggplot(expr_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  theme_minimal() +
  stat_compare_means(aes(group = Group, label = ..p.signif..),
                     method = "wilcox.test") +
  labs(title = "GSE77938 Gene Expression in KC vs Control Samples",
       x = "Group",
       y = "Expression") +
  scale_fill_manual(values = c("Control" = "#0468BF", "Disease" = "#BF0449")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# 显示图形
print(p1)

# 计算每个基因在Control和Disease组中的表达均值
mean_expression <- expr_data_long %>%
  group_by(Gene, Group) %>%
  summarise(mean_expression = mean(Expression, na.rm = TRUE))

# 合并均值数据和p值数据
combined_data <- merge(mean_expression, p_values, by = "Gene")

# 筛选Control组比Disease组表达值高且显著的基因
significant_genes <- combined_data %>%
  spread(Group, mean_expression) %>%
  filter(Control > Disease & p_value < 0.05)

# 查看结果
print(significant_genes)

# 提取Overlapgenes的表达数据
expr_data <- GSE151631_exp[feature_gene, ]

# 转置表达数据以便每行代表一个样本
expr_data <- t(expr_data)

# 将表达数据转换为数据框
expr_data <- as.data.frame(expr_data)

sample_names <-rownames(expr_data)
# 创建分组向量
group <- ifelse(grepl("^DN|^LE", sample_names), "Control", "Disease")


# 添加样本信息
expr_data$SampleID <- rownames(expr_data)
expr_data$Group <- group

# 计算每个基因的Wilcoxon检验p值
p_values <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group) %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Expression ~ Group)$p.value)

# 查看结果
print(p_values)
# 合并表达数据和p值
expr_data_long <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group)

expr_data_long <- merge(expr_data_long, p_values, by = "Gene")


p2 <- ggplot(expr_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  theme_minimal() +
  stat_compare_means(aes(group = Group, label = ..p.signif..),
                     method = "wilcox.test") +
  labs(title = "GSE77938 Gene Expression in KC vs Control Samples",
       x = "Group",
       y = "Expression") +
  scale_fill_manual(values = c("Control" = "#0468BF", "Disease" = "#BF0449")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

library(gridExtra)
grid.arrange(p1, p2, ncol = 1)

# Assuming p1 and p2 are your ggplot objects
combined_plot <- arrangeGrob(p1, p2, ncol = 1)


Candidate_biomarker <-c('GADD45B','GATA6','RBM24')
saveRDS(Candidate_biomarker,'./01_Candidate_biomarker.rds')

# 提取Overlapgenes的表达数据
expr_data <- GSE77938_exp[Candidate_biomarker, ]

# 转置表达数据以便每行代表一个样本
expr_data <- t(expr_data)

# 将表达数据转换为数据框
expr_data <- as.data.frame(expr_data)

sample_names <-rownames(expr_data)

# 创建分组向量
group <- ifelse(grepl("^KC", sample_names), "Disease", "Control")

# 添加样本信息
expr_data$SampleID <- rownames(expr_data)
expr_data$Group <- group

# 计算每个基因的Wilcoxon检验p值
p_values <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group) %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Expression ~ Group)$p.value)

# 查看结果
print(p_values)
# 合并表达数据和p值
expr_data_long <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group)

expr_data_long <- merge(expr_data_long, p_values, by = "Gene")

# 绘制基因表达柱状图并添加检验
p1 <- ggplot(expr_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  stat_compare_means(aes(group = Group, label = ..p.signif..),
                     method = "wilcox.test") +
  labs(title = "GSE77938 Gene Expression in KC vs Control Samples",
       x = "Group",
       y = "Expression") +
  scale_fill_manual(values = c("Control" = "#0468BF", "Disease" = "#BF0449")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# 显示图形
print(p1)

# 计算每个基因在Control和Disease组中的表达均值
mean_expression <- expr_data_long %>%
  group_by(Gene, Group) %>%
  summarise(mean_expression = mean(Expression, na.rm = TRUE))

# 合并均值数据和p值数据
combined_data <- merge(mean_expression, p_values, by = "Gene")

# 筛选Control组比Disease组表达值高且显著的基因
significant_genes <- combined_data %>%
  spread(Group, mean_expression) %>%
  filter(Control > Disease & p_value < 0.05)

# 查看结果
print(significant_genes)

# 提取Overlapgenes的表达数据
expr_data <- GSE151631_exp[Candidate_biomarker, ]

# 转置表达数据以便每行代表一个样本
expr_data <- t(expr_data)

# 将表达数据转换为数据框
expr_data <- as.data.frame(expr_data)

sample_names <-rownames(expr_data)
# 创建分组向量
group <- ifelse(grepl("^DN|^LE", sample_names), "Control", "Disease")


# 添加样本信息
expr_data$SampleID <- rownames(expr_data)
expr_data$Group <- group

# 计算每个基因的Wilcoxon检验p值
p_values <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group) %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Expression ~ Group)$p.value)

# 查看结果
print(p_values)
# 合并表达数据和p值
expr_data_long <- expr_data %>%
  gather(Gene, Expression, -SampleID, -Group)

expr_data_long <- merge(expr_data_long, p_values, by = "Gene")


p2 <- ggplot(expr_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  stat_compare_means(aes(group = Group, label = ..p.signif..),
                     method = "wilcox.test") +
  labs(title = "GSE77938 Gene Expression in KC vs Control Samples",
       x = "Group",
       y = "Expression") +
  scale_fill_manual(values = c("Control" = "#0468BF", "Disease" = "#BF0449")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

library(gridExtra)
grid.arrange(p1, p2, ncol = 1)

# Assuming p1 and p2 are your ggplot objects
combined_plot <- arrangeGrob(p1, p2, ncol = 1)

ggsave(filename = './01_combined_plot_Exp_plot.pdf',combined_plot,width = 6,height = 7)

ggsave(filename = './01_combined_plot_Exp_plot.png',combined_plot,width = 6,height = 7)


# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'Rdata_',dir_name, ".RData"))