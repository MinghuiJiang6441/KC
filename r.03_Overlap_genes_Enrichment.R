rm(list = ls())
gc()

setwd('/data/home/jiangminghui/KC/Project-0035/03_Overlap_genes_Enrichment/')
library(openxlsx)
library(ggvenn)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
load('./03_Overlap_genes_Enrichment.RData')
load('../03_Overlap_genes_Enrichment/03_Overlap_genes_Enrichment.RData')
load('../02_WGCNA/06_WGCNA_keyGene.RData')
data <- read.xlsx("./tableS1.xlsx", sheet = 2)  # 读取第一个工作表
# 将数据框转换为矩阵
data_matrix <- as.matrix(data)

# 将矩阵转换为向量
genes_vector <- as.vector(data_matrix)

# 去除NA值
genes_vector <- genes_vector[!is.na(genes_vector)]

# 去除重复的基因名称
ISR_RGs <- unique(genes_vector)

# 查看提取的基因
print(ISR_RGs)

save(ISR_RGs,file = './00_ISR_RGs,rdata')

load('./00_ISR_RGs,rdata')

res <-readRDS('../01_KC-DEGs/01_DEGs_res.rds')

# 添加显著性标记
res$significant <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0.5, "Upregulated",
                          ifelse(res$padj < 0.05 & res$log2FoldChange < -0.5, "Downregulated", "Not Significant"))


res_df <- as.data.frame(res)

KC_DEGs <- subset(res_df, res_df$significant != "Not Significant") %>% rownames()
head(KC_DEGs)

venn <- list('KC_DEGs' =KC_DEGs,
             'KC_RGs' = WGCNA_keyGene,
             'ISR_RGs' =  ISR_RGs
)

str(venn)
ggvenn(venn)    

pdf('./01_Venn_Plot.pdf',width = 6,height = 6)
ggvenn(venn)    
dev.off()

png('./01_Venn_Plot.png',width = 6,height = 6,units = 'in',res = 300)
ggvenn(venn)    
dev.off()

Overlapgenes <- intersect(KC_DEGs,intersect(WGCNA_keyGene,ISR_RGs))
Overlapgenes

saveRDS(Overlapgenes,'./01_Overlapgenes.rds')


############################################################
#
#-------------------富集分析———————————————————————————
#
############################################################

Overlap_genes <-readRDS('./01_Overlapgenes.rds')

Overlap_genes_ENTREZID = bitr(Overlap_genes, #数据集
                              fromType="SYMBOL", #输入为SYMBOL格式
                              toType="ENTREZID",  # 转为ENTERZID格式
                              OrgDb="org.Hs.eg.db") #人类 数据库


head(Overlap_genes_ENTREZID)

ego_ALL <- enrichGO(gene = Overlap_genes_ENTREZID$ENTREZID,
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) 

saveRDS(ego_ALL,'./02_GO_res.rds')

write.csv(ego_ALL,'./02_GO_res.csv')


# 筛选p值小于0.05的通路
# 提取包含p值的信息
ego_df <- as.data.frame(ego_ALL)

# 筛选p值小于0.05的通路
ego_filtered <- ego_df[ego_df$p.adjust < 0.05, ]

# 确保ego_filtered仍然是一个enrichResult对象
ego_filtered <- ego_ALL
ego_filtered@result <- ego_filtered@result[ego_filtered@result$p.adjust < 0.05, ]

# 绘制点图
p <- dotplot(ego_filtered, showCategory = 10, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  theme(axis.text.y = element_text(size = 9))  # 调整 y 轴字体大小

# 显示图形
print(p)

pdf('./02_GOEnrichment_Plot.pdf',height =8,width = 6)
print(p)

dev.off()

png('./02_GOEnrichment_Plot.png',height =8,width = 6,units = 'in',res = 300)
print(p)

dev.off()

# 提取结果并按 p 值排序
ego_results <- as.data.frame(ego_ALL)

# 分别筛选 CC、BP 和 MF 中 p 值小于 0.05 的通路，并取前 10 个
top10_CC <- ego_results %>% filter(ONTOLOGY == "CC" & p.adjust < 0.05) %>% arrange(p.adjust) %>% head(10)
top10_BP <- ego_results %>% filter(ONTOLOGY == "BP" & p.adjust < 0.05) %>% arrange(p.adjust) %>% head(10)
top10_MF <- ego_results %>% filter(ONTOLOGY == "MF" & p.adjust < 0.05) %>% arrange(p.adjust) %>% head(10)

top10_BP$Description
# [1] "cellular response to decreased oxygen levels"      "cellular response to oxygen levels"               
# [3] "response to decreased oxygen levels"               "response to oxygen levels"                        
# [5] "cellular response to hypoxia"                      "response to hypoxia"                              
# [7] "endocardial cushion development"                   "positive regulation of mitochondrion organization" 

top10_CC$Description
# [1] "caveola"                      "plasma membrane raft"         "mitochondrial outer membrane"
# [4] "Bcl-2 family protein complex" "organelle outer membrane"     "outer membrane"              
# [7] "distal axon"                  "membrane raft"                "membrane microdomain"        
# [10] "calyx of Held"                    

top10_MF$Description
# [1] "phosphatidylinositol 3-kinase regulatory subunit binding"       
# [2] "type I transforming growth factor beta receptor binding"        
# [3] "activin binding"                                                
# [4] "transforming growth factor beta receptor activity"              
# [5] "pre-mRNA intronic binding"                                      
# [6] "sequence-specific mRNA binding"                                 
# [7] "potassium ion leak channel activity"                            
# [8] "transmembrane receptor protein serine/threonine kinase activity"
# [9] "leak channel activity"                                          
# [10] "narrow pore channel activity"   
cc_count <- sum(ego_ALL@result$ONTOLOGY == "CC")
bp_count <- sum(ego_ALL@result$ONTOLOGY == "BP")
mf_count <- sum(ego_ALL@result$ONTOLOGY == "MF")


genelist <- as.numeric(res_df[Overlap_genes,2]) 
names(genelist) <- row.names(res_df[Overlap_genes,])

cnetp2 <- cnetplot(ego_ALL, 
                  color.params = list(foldChange = genelist, edge = TRUE),
                  showCategory = 20,
                  node_label = 'gene',
                  circular = TRUE) + 
  theme(legend.text = element_text(size = 5),  # 调整图例字体大小
        plot.title = element_text(size = 5))  # 调整标题字体大小


png("02_GO_cnetplot.png", width = 5, height = 5,units = 'in',res = 300)

cnetp2
dev.off()

pdf("02_GO_cnetplot.pdf", width = 5, height = 5)

cnetp2
dev.off()

ego_filtered@result

############################################################
#
#--------------------KEGG富集分析————————————————————————————
#
############################################################

kk <-  enrichKEGG(
  gene          = Overlap_genes_ENTREZID$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 1
)

# #没有显著富集KEGG通路

 pdf('./03_KEGGEnrichment_Plot.pdf',height =5,width =7)
 dotplot(kk, showCategory = 15)
 dev.off()
#
 png('./03_KEGGEnrichment_Plot.png',height = 5,width = 7,res = 300 ,units = 'in')
 dotplot(kk, showCategory = 15)
 dev.off()

 gene_ids <- unlist(strsplit(kk$geneID, "/"))
 
 gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
 
 kk_df <- as.data.frame(kk)
 
 kk@result$geneSymbol <- sapply(strsplit(kk_df$geneID, "/"), function(ids) paste(gene_symbols[ids], collapse = "/"))
 
 gene_symbols <- mapIds(org.Hs.eg.db, keys =  kk@gene, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
   
 kk@gene <-gene_symbols
   
 netp_KEGG <- cnetplot(kk,
                   color.params = list(foldChange = genelist, edge = TRUE),
                   showCategory =10,
                   node_label = 'gene',
                   circular = TRUE)

 png("03_KEGG_cnetplot.png", width = 8, height = 8,units = 'in',res = 300)

 netp_KEGG
 dev.off()
#
 pdf("03_KEGG_cnetplot.pdf", width = 8, height = 8)

 netp_KEGG
 dev.off()


# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", dir_name, ".RData"))