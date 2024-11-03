
rm(list = ls())
gc()
setwd('/data/home/jiangminghui/KC/Project-0035/01_KC-DEGs/')
library(DESeq2)
library(ggrepel)
library(ComplexHeatmap)
library(limma)

library(ggplot2)
library(circlize)
library(edgeR)

############################################################
#
#------------------DEGs -----------------————————————————
#
############################################################

GSE77938_exp <-readRDS('../00_raw_data/combat_data.rds')
# 获取列名
sample_names <- colnames(GSE77938_exp)

# 创建分组向量
group <- ifelse(grepl("^KC", sample_names), "Disease", "Control")

colData <- data.frame(
  row.names = sample_names,
  condition = group
)
# 查看分组向量
head(colData)

dds <- DESeqDataSetFromMatrix(countData = GSE77938_exp, colData = colData, design = ~ condition)

# 过滤低表达基因
dds <- dds[rowSums(counts(dds)) > 1, ]

# 差异表达分析
dds <- DESeq(dds)

# 提取结果
res <- results(dds)
# 查看结果
head(res)
res <-na.omit(res)
# design <- model.matrix(~ 0 + group)
# colnames(design) <- c("Control", "Disease")
# GSE77938_exp[GSE77938_exp < 0] <- 0
# 
# # 创建DGEList对象
# dge <- DGEList(counts = GSE77938_exp)
# 
# # 计算标准化因子
# dge <- calcNormFactors(dge)
# 
# # 进行voom转换
# v <- voom(dge, design)
# 
# # 拟合线性模型
# fit <- lmFit(v, design)
# 
# # 创建对比矩阵
# contrast.matrix <- makeContrasts(Disease - Control, levels = design)
# 
# # 应用对比矩阵
# fit2 <- contrasts.fit(fit, contrast.matrix)
# 
# # 进行贝叶斯检验
# fit2 <- eBayes(fit2)
# res <- topTable(fit2, adjust = "fdr", number = Inf)


saveRDS(res,'./01_DEGs_res.rds')
res <-readRDS('./01_DEGs_res.rds')
# 添加显著性标记
res$significant <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0.5, "Upregulated",
                          ifelse(res$padj < 0.05 & res$log2FoldChange < -0.5, "Downregulated", "Not Significant"))

table(res$significant)
# 将res对象转换为数据框
res_df <- as.data.frame(res)

# 筛选出上调基因中log2FoldChange最大的前20个基因
top_upregulated <- res_df %>%
  filter(significant == 'Upregulated') %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

# 筛选出下调基因中log2FoldChange最小的前20个基因
top_downregulated <- res_df %>%
  filter(significant == 'Downregulated') %>%
  arrange(log2FoldChange) %>%
  head(20)

# 合并上调和下调的基因
top_genes <- bind_rows(top_upregulated, top_downregulated)
top_genes$gene <- rownames(top_genes)

# 绘制火山图
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.8, size = 1.75) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "#A60A27", "Downregulated" = "#D9A60D")) +
  theme_minimal() +
  labs(title = "Volcano Plot of KC-DEGs", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +  # 添加log2FoldChange阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # 添加p值阈值线
  xlim(-20, 10) +  # 设置x轴范围
  ylim(0, 37) +  # 设置y轴范围
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    panel.background = element_blank(),  # 去除背景
    axis.line = element_line(color = "black")  # 添加轴线
  )+ 
  geom_text_repel(data = top_genes, aes(label = gene), 
                  size = 3, max.overlaps = Inf, 
                  box.padding = 0.5, point.padding = 0.5, 
                  segment.color = 'grey50', color = "black")  # 调整标签
# 显示火山图
print(volcano_plot)

 ggsave('./01_DEGs_volcano_plot.png',plot =volcano_plot,width = 9,height = 7 )
 ggsave('./01_DEGs_volcano_plot.pdf',plot =volcano_plot,width = 9,height = 7 )
 
 range(res_df$padj)
 top_genes_exp <- GSE77938_exp[rownames(GSE77938_exp) %in% top_genes$gene, ] %>%as.matrix()
 
saveRDS(top_genes_exp,'./02_top_genes_exp.rds')
############################################################
#
#------------------圆形热图 -----------------————————————————
#
############################################################

# 假设 GSE77938_exp 是一个矩阵，top_genes 是一个包含基因名的数据框
# 提取 top_genes 中的基因表达数据

top_genes_exp <-readRDS('./02_top_genes_exp.rds')
colnames(top_genes_exp)
sample_names_sorted <- c(sort(grep("^KC", colnames(top_genes_exp), value = TRUE)),
                         sort(grep("^KR", colnames(top_genes_exp), value = TRUE)))

# 按照排序后的列名重新排序表达矩阵
top_genes_exp_sorted <- top_genes_exp[, sample_names_sorted]
# # 标准化数据
# cir1 <- scale(top_genes_exp)
# range(cir1)
# 
# # 设置legend颜色
# mycol1 = colorRamp2(c(-1, 0, 5), c("#3068D9", "white", "#D93E30"))
# dev.new()
# # 绘制热图并保存为PDF
# pdf('./03_circos_Heatmap.pdf',width = 5,height = 5)
# circos.heatmap(cir1, col = mycol1, dend.side = "inside", # 聚类放在环形内侧
#                rownames.side = "outside")
# dev.off()
# 
# # 完全清除布局
# circos.clear()
# 
# png('./03_circos_Heatmap.png',width = 5,height = 5,units = 'in',res = 300)
# circos.heatmap(cir1, col = mycol1, dend.side = "inside", # 聚类放在环形内侧
#                rownames.side = "outside")
# dev.off()
split <- factor(c(rep("Group1", nrow(top_upregulated)), rep("Group2", nrow(top_downregulated))))

Heatmap(top_genes_exp, row_split = split)


mat1 = rbind(cbind(matrix(rnorm(50*5, mean = 1), nr = 50), 
                   matrix(rnorm(50*5, mean = -1), nr = 50)),
             cbind(matrix(rnorm(50*5, mean = -1), nr = 50), 
                   matrix(rnorm(50*5, mean = 1), nr = 50))
)

split <- sample(letters[1:5], 100, replace = TRUE)

split <- factor(split, levels = letters[1:5])

Heatmap(mat1, row_split = split)
col_fun <- colorRamp2(c(-2, 0, 2), c("#fc8d59", "#ffffbf", "#91bfdb"))
circos.heatmap(mat1, split = split, col = col_fun)
circos.clear()
circos.par(start.degree = 90, gap.degree = 10)
circos.heatmap(
  mat1, split = split, col = col_fun, 
  track.height = 0.4, bg.border = "red", 
  bg.lwd = 2, bg.lty = 2, show.sector.labels = TRUE
)
circos.clear()



dim(top_genes_exp)
col_fun <- colorRamp2(c(-2, 0, 2), c("#F87174", "#F2F5FA", "#5E8EC7"))

# 获取列名并排序
sorted_colnames <- colnames(top_genes_exp)[order(colnames(top_genes_exp))]

sorted_top_genes_exp <- top_genes_exp[, sorted_colnames] 

gene_groups <- c(rep("Upregulated", nrow(top_upregulated)), rep("Downregulated", nrow(top_downregulated)))

top_genes_exp_zscore <- as.data.frame(scale(top_genes_exp_sorted))

circos.clear()
#pdf('./03_circos_Heatmap.pdf',width = 7,height = 7)
png('./03_circos_Heatmap.png',width = 7,height = 7,res = 300,units = 'in')
circos.par(start.degree = 90, gap.degree = 10)
circos.heatmap(
  standardized_data[,1:20], split = gene_groups, col = col_fun, 
  track.height = 0.4, bg.border = "red", 
  bg.lwd = 2, bg.lty = 2, show.sector.labels = TRUE
)

circos.heatmap(
  standardized_data[,20:40], split = gene_groups, col = col_fun, 
  track.height = 0.4, bg.border = "red", 
  bg.lwd = 2, bg.lty = 2, show.sector.labels = TRUE
)
dev.off()



circos.clear()

range(standardized_data)
color <- colorRamp2(c(-0.5, 0, 7), c("blue", "white", "green"))
circos.par(gap.after = c(20))#间隔
circos.heatmap(standardized_data, #数据
               col = color,#颜色
               dend.side = "inside",#确定聚类结果放在圈内还是圈外
               rownames.side = "outside",#组名
               track.height = 0.4
               # clustering.method = "complete",#归一化处理
               # distance.method = "euclidean"#聚类方法，默认为欧氏距离
)
grid.draw(Legend(title = "Title", col_fun = color))

expr_scaled <- t(scale(t(expr[top_genes$gene,])))
group_colors <- list(Group = c("Control" = "#049DBF", "Disease" = "#8C030E"))

ha <- HeatmapAnnotation(df = data.frame(Group = group), col = group_colors)

# 绘制热图
heatmap <- Heatmap(standardized_data, 
                   name = "Expression", 
                   top_annotation = ha, 
                   show_row_names = T, 
                   show_column_names = FALSE,
                   cluster_rows = F, 
                   cluster_columns = F,
                   col = colorRamp2(c(min(standardized_data), 0, max(standardized_data)), c("#03178C", "white", "#A52226")))

# 显示热图
draw(heatmap)

pdf('./02_DEGs_heatmap.pdf', width = 14, height = 9)
draw(heatmap)
dev.off()

png('./02_DEGs_heatmap.png', width = 14, height = 9, units = "in", res = 300)
draw(heatmap)
dev.off()
