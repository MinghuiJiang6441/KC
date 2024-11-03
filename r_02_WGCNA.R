rm(list =ls())


setwd('/data/home/jiangminghui/KC/Project-0035/02_WGCNA/')
#.libPaths('/data/nas2/Software/miniconda3/envs/public_R/lib/R/library')
library(WGCNA)
library(reshape2)
library(stringr)
library(limma)
library(sva)
library(parallel)
library(foreach)
library(doParallel)


GSE77938_exp <-readRDS('../00_raw_data/combat_data.rds')
dim(GSE77938_exp)

load('./WGCNA.rdata')

type = "unsigned"
corType = "pearson"

corFnc = if (corType == "pearson") {
  corFnc = cor
} else {
  corFnc = bicor
}

maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

range(GSE77938_exp)
#expr1 <- log2(GSE77938_exp+1)
expr1 <-GSE77938_exp
range(expr1)

dim(expr1)

m.mad <- apply(expr1,1,mad,na.rm=T)

dataExprVar <- expr1[which(m.mad > 
                             max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]



## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
table(gsg$goodSamples)
gsg$allOK


if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
head(dataExpr)[,1:8]
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h = 3000000, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 3000000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.

keepSamples = (clust==1)


dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")

pdf('./01_sampleTree.pdf',width = 10,height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

png('./01_sampleTree.png',width = 10,height = 6,units = 'in',res = 300)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

#options(mc.cores = 30)

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)



par(mfrow = c(1,2))
cex1 = 0.9
power = sft$powerEstimate
power
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
pdf('./02_WGCNA_Soft_threshold_power.pdf',width = 8,height = 6)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")
dev.off()

png('./02_WGCNA_Soft_threshold_power.png',width = 8,height = 6,res = 300,units = 'in')
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")
dev.off()


pdf('./03_WGCNA_mean_Connectivity.pdf',width = 8,height = 6 )
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()


png('./03_WGCNA_mean_Connectivity.png',width = 8,height = 6 ,res = 300,units = 'in')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()

save.image('./WGCNA.rdata')


net <-blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                   TOMType = type, minModuleSize = 500,
                   reassignThreshold = 0, mergeCutHeight = 0.25,
                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                   saveTOMs = F, corType = corType, 
                   maxPOutliers = maxPOutliers, loadTOMs = F,
                   verbose = 3)

saveRDS(net, file = './04_WGCNA_net.rds')

save.image('./WGCNA.rdata')

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

pdf('./04_WGCNA_Cluster_tree.pdf',width = 8,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

png('./04_WGCNA_Cluster_tree.png',width = 8,height = 6,res = 300,units = 'in')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

##绘制模块之间的相关性图
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf('./Sup_Fig_WGCNA_Correlation_modules.pdf',height = 10,width = 8)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

sample_names <- rownames(dataExpr)

# 创建分组向量
group <- ifelse(grepl("^KC", sample_names), "Disease", "Control")

# 创建数据框
df <- data.frame(Disease = ifelse(group == "Disease", 1, 0),
                 Control = ifelse(group == "Control", 1, 0),
                 row.names = sample_names)

# 查看数据框
print(df)
df <-df[rownames(dataExpr),]

MEs_colpheno <- orderMEs(cbind(MEs_col, df))

#明确样本数和基因数
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(dataExpr), method = "average")


par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)


# 将 df$Disease 和 df$Control 转换为因子
group_factor <- factor(paste0("Disease_", df$Disease, "_Control_", df$Control))

# 生成颜色
sample_colors <- numbers2colors(as.numeric(group_factor), 
                                colors = c("#F20505", "white"), signed = FALSE)

pdf('./01_Sample_Cluster_Map.pdf',width = 10,height = 8)
#par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

dev.off()

png('./01_Sample_Cluster_Map.png',width = 10,height = 8,units = 'in',res = 300)
#par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

dev.off()


# 创建设计矩阵
design <- model.matrix(~ df$Disease + df$Control)
colnames(design) <- c("Intercept", "Disease", "Control")

moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(dataExpr, moduleColors)$eigengenes

# 保存模块特征值
# saveRDS(MEs0, file = './result/data/03_WGCNA_MEs0.rds')

MEs <- orderMEs(MEs0)  # 不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor <- cor(MEs, design[, -1], use = "p")  # 排除 Intercept 列
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 设置窗口大小
sizeGrWindow(10, 6)

# 显示相关性和p值
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
pdf('./05_WGCNA_Moduletrait_relationships.pdf',width = 8,height = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design)[-1],  # 排除 Intercept 列
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()

png('./05_WGCNA_Moduletrait_relationships.png',width = 8,height = 10,units = 'in',res = 300)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design)[-1],  # 排除 Intercept 列
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#module = c('magenta','blue','green')
module = c('yellow')
#module = "grey60"
pheno = "Disease"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(df))
# 获取模块内的基因
moduleGenes = moduleColors == module

table(moduleGenes)
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, df, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, df, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('./06_WGCNA_Scatterplot_GS_MM_yellow.pdf',width = 6,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#abline(h=0.6,v = 0.6,col="red")

dev.off()

# 与性状高度相关的基因，也是与性状相关的模型的关键基因
png('./06_WGCNA_Scatterplot_GS_MM_yellow.png',width = 6,height = 6,units = 'in',res = 300)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#abline(h=0.5,v = 0.5,col="red")

dev.off()


# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#module = c('magenta','blue','green')
module = c('green')
#module = "grey60"
pheno = "Control"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(df))
# 获取模块内的基因
moduleGenes = moduleColors == module

table(moduleGenes)
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, df, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, df, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('./06_WGCNA_Scatterplot_GS_MM_green.pdf',width = 6,height = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#abline(h=0.6,v = 0.6,col="red")

dev.off()

# 与性状高度相关的基因，也是与性状相关的模型的关键基因
png('./06_WGCNA_Scatterplot_GS_MM_green.png',width = 6,height = 6,units = 'in',res = 300)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in",  "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#abline(h=0.5,v = 0.5,col="red")

dev.off()


module = c('yellow','green')
#module = "grey60"
pheno = "Disease"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(df))
# 获取模块内的基因
moduleGenes = moduleColors == module
module_Genes <- rownames(geneModuleMembership[moduleGenes,])

length(module_Genes)

WGCNA_keyGene <- module_Genes

length(WGCNA_keyGene)

save(WGCNA_keyGene,file = './06_WGCNA_keyGene.RData')

save.image('./WGCNA.rdata')
