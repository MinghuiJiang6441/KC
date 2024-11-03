rm(list = ls())
gc()
library(tidyverse)
library(GEOquery)
library(tinyarray)
library(biomaRt)
library(limma)
library(ggdendro)
library(ggplot2)
library(sva)

setwd('/data/home/jiangminghui/KC/Project-0035/00_raw_data/')
# 芯片数据
############################################################
#
#-------------------------GSE77938--------————————————————
#
############################################################
geo_number <- "GSE77938"

con <- gzfile("./GSE77938/GSE77938_discovery_gene_counts.txt.gz", "rt")
GSE77938_discovery <- read.table(con,header = T,comment.char = "!",row.names = 1)
# 关闭连接
close(con)

dim(GSE77938_discovery)
GSE77938_discovery[1:4,1:4]

con <- gzfile("./GSE77938/GSE77938_replication_transcript_counts.txt.gz", "rt")
GSE77938_replication <- read.table(con,header = T,comment.char = "!",row.names = 1)
# 关闭连接
close(con)

GSE77938_replication[1:4,1:4]
dim(GSE77938_replication)
colnames(GSE77938_replication)
range(GSE77938_replication)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取基因符号
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(GSE77938_discovery),
  mart = ensembl
)

# 检查并处理重复的基因符号
gene_symbols <- gene_symbols[gene_symbols$hgnc_symbol != "", ]  # 移除空的基因符号
gene_symbols <- gene_symbols[!duplicated(gene_symbols$ensembl_gene_id), ]  # 移除重复的Ensembl基因ID

# 创建一个映射表
gene_map <- setNames(gene_symbols$hgnc_symbol, gene_symbols$ensembl_gene_id)

# 替换矩阵的行名
new_rownames <- gene_map[rownames(GSE77938_discovery)]
new_rownames[is.na(new_rownames)] <- rownames(GSE77938_discovery)[is.na(new_rownames)]  # 保留没有匹配到的原始行名

# 确保行名唯一
new_rownames <- make.unique(new_rownames)

# 设置新的行名
rownames(GSE77938_discovery) <- new_rownames

GSE77938_replication[1:4,1:4]

# 获取基因符号
gene_symbols <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = rownames(GSE77938_replication),
  mart = ensembl
)

# 检查并处理重复的基因符号
gene_symbols <- gene_symbols[gene_symbols$hgnc_symbol != "", ]  # 移除空的基因符号
gene_symbols <- gene_symbols[!duplicated(gene_symbols$ensembl_transcript_id), ]  # 移除重复的Ensembl基因ID

# 创建一个映射表
gene_map <- setNames(gene_symbols$hgnc_symbol, gene_symbols$ensembl_transcript_id)

# 替换矩阵的行名
new_rownames <- gene_map[rownames(GSE77938_replication)]
new_rownames[is.na(new_rownames)] <- rownames(GSE77938_replication)[is.na(new_rownames)]  # 保留没有匹配到的原始行名

# 确保行名唯一
new_rownames <- make.unique(new_rownames)

# 设置新的行名
rownames(GSE77938_replication) <- new_rownames

GSE77938_replication[1:4,1:4]

saveRDS(GSE77938_replication,'./00_GSE77938_replication.rds')

saveRDS(GSE77938_discovery,'./00_GSE77938_discovery.rds')

GSE77938_replication <-readRDS('./00_GSE77938_replication.rds')
GSE77938_discovery <-readRDS('./00_GSE77938_discovery.rds')
# 获取共有基因
both_genes <- intersect(rownames(GSE77938_replication), rownames(GSE77938_discovery))

# 根据共有基因合并数据框
GSE77938_replication_common <- GSE77938_replication[both_genes, ]
GSE77938_discovery_common <- GSE77938_discovery[both_genes, ]

# 合并数据框
combined_data <- cbind(GSE77938_replication_common, GSE77938_discovery_common)


# 过滤表达过低的基因
threshold <- 0  # 设置表达量阈值
filtered_data <- combined_data[rowMeans(combined_data) > threshold, ]
dim(filtered_data)

saveRDS(filtered_data ,'./GSE77938_exp.rds') 

sample_dist <- dist(t(filtered_data))

# 进行层次聚类
sample_clustering <- hclust(sample_dist)

# 绘制聚类树
sample_groups <- rep(c("Batch1", "Batch2"), times = c(ncol(GSE77938_replication), ncol(GSE77938_discovery)))
group_colors <- c("Batch1" = "blue", "Batch2" = "red")

pdf('./Sup_Batch_effect_Clustering.pdf',width = 11,height = 6)
plot(sample_clustering,　 main = "Hierarchical Clustering of Samples", xlab = "", sub = "", cex = 0.9)

rect.hclust(sample_clustering, k = 2, border = c("blue", "red"))
dev.off()


png('./Sup_Batch_effect_Clustering.png',width = 11,height = 6,res = 300,units = 'in')
plot(sample_clustering,　 main = "Hierarchical Clustering of Samples", xlab = "", sub = "", cex = 0.9)

rect.hclust(sample_clustering, k = 2, border = c("blue", "red"))
dev.off()
group <- ifelse(grepl("^KC", colnames(combined_data)), "Disease", "Control")


　# 合并数据框
combined_data <- cbind(GSE77938_replication_common, GSE77938_discovery_common)

# 过滤表达过低的基因
threshold <- 0  # 设置表达量阈值
filtered_data <- combined_data[rowMeans(combined_data) > threshold, ]

# 创建批次信息
batch <- c(rep(1, ncol(GSE77938_replication_common)), rep(2, ncol(GSE77938_discovery_common)))

# 使用ComBat消除批次效应
filtered_data <-as.matrix(filtered_data)
combat_data_seq <- ComBat_seq(filtered_data, batch = batch,group = group)

saveRDS(combat_data_seq,'./combat_data.rds')

range(combat_data_seq)
sample_dist <- dist(t(combat_data_seq))

# 进行层次聚类
sample_clustering <- hclust(sample_dist)

pdf('./Sup_Batch_effect_Clustering_after.pdf',width = 11,height = 6)

plot(sample_clustering,　 main = "Hierarchical Clustering of Samples", xlab = "", sub = "", cex = 0.9)
rect.hclust(sample_clustering, k = 2, border = c("blue", "red"))
dev.off()

png('./Sup_Batch_effect_Clustering_after.png',width = 11,height = 6,res = 300,units = 'in')
plot(sample_clustering,　 main = "Hierarchical Clustering of Samples", xlab = "", sub = "", cex = 0.9)

rect.hclust(sample_clustering, k = 2, border = c("blue", "red"))
dev.off()
############################################################
#
#-------------------------GSE151631--------————————————————
#
############################################################

con <- gzfile("./GSE151631/GSE151631_Raw.counts.csv.gz", "rt")
GSE151631 <- read.csv(con,header = T,comment.char = "!")
# 关闭连接
close(con)

# 去重第一列
GSE151631_unique <- GSE151631[!duplicated(GSE151631[, 1]), ]

# 将第一列作为行名
rownames(GSE151631_unique) <- GSE151631_unique[, 1]
GSE151631_unique <- GSE151631_unique[, -1]

# 过滤表达过低的基因
threshold <- 0  # 设置表达量阈值
filtered_data <- GSE151631_unique[rowMeans(GSE151631_unique) > threshold, ]

saveRDS(filtered_data, './GSE151631_exp.rds')


