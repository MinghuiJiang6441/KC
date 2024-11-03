rm(list = ls())
gc()
library(utils)
library(glmnet)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(psych)
library(Hmisc)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)
library(enrichplot)
library(limma)
library(ggplot2) 
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(msigdbr)
library(cowplot)
library(enrichplot)

setwd('/data/home/jiangminghui/KC/Project-0035/10_GSEA_enrichment_analysis/')
load('./Rdata_10_GSEA_enrichment_analysis.RData')
GSE77938 <-readRDS('../00_raw_data/combat_data.rds')
Candidate_biomarker <-readRDS('../07_ValidateExpressions/01_Candidate_biomarker.rds')



msigdb <- msigdbr(species = "Homo sapiens", category = "H")
msigdb_filtered <- msigdb %>% dplyr::select(gs_name, gene_symbol)

start_time <- Sys.time()
print(paste("开始时间:", start_time))
gsea_results_list <- list()

for (i in 1:length(Candidate_biomarker)) {
  # 提取当前预后基因数据
  prognostic_gene <- GSE77938[Candidate_biomarker[i], ]  %>% as.data.frame()
  
  # 提取所有基因数据
  all_genes <- GSE77938[!rownames(GSE77938) %in% Candidate_biomarker[i], ] %>% t()
  
  # 计算相关性
  cor_results <- corr.test(prognostic_gene, all_genes, method = "spearman")
  cor_coefficients <- cor_results$r
  
  # 准备基因列表
  geneList <- as.numeric(cor_coefficients[1, ])
  names(geneList) <- colnames(cor_coefficients)
  geneList <- sort(geneList, decreasing = TRUE)
  
  # 进行GSEA分析
  gsea_results <- GSEA(geneList, TERM2GENE = msigdb_filtered, pvalueCutoff = 0.05)
  
  # 将结果存储到列表中
  gsea_results_list[[Candidate_biomarker[i]]] <- gsea_results
  
  print(i)
}

end_time <- Sys.time()
print(paste("结束时间:", end_time))

names(gsea_results_list)


saveRDS(gsea_results_list, './01_gsea_results_list.rds')

gsea_results_list<-readRDS('./01_gsea_results_list.rds')


top5_plots <- list()

for (i in 1:length(gsea_results_list)) {

    significant_pathways <- gsea_results_list[[i]][gsea_results_list[[i]]$NES > 1, ]
  

      top5_pathways <- significant_pathways[order(significant_pathways$p.adjust), ][1:5, ]
  

    pathway <- top5_pathways[1:5, "ID"]
    plot <- gseaplot2(gsea_results_list[[i]], geneSetID = pathway, title = paste(Candidate_biomarker[i],'_', "Top 5 Pathways"))
    
    
    top5_plots[[paste("Candidate biomarker", i)]] <- plot
}



saveRDS(top5_plots,'./01_top5_plots.rds')

pdf("./01_GSEA_Top5_Plots.pdf", width = 8, height = 11)
for (i in 1:length(top5_plots)) {
  print(top5_plots[[i]])
}
dev.off()

for (i in 1:length(top5_plots)) {
  png(filename = paste0("01_GSEA_Top5_Plots_Candidate_biomarker_", Candidate_biomarker[i], ".png"), width = 8, height = 11,units = 'in',res = 300)
  print(top5_plots[[i]])
  dev.off()
}

pathway_ls<-list()
for (i in 1:length(gsea_results_list)) {
  
  significant_pathways <- gsea_results_list[[i]][gsea_results_list[[i]]$NES > 1, ]
  
  
  top5_pathways <- significant_pathways[order(significant_pathways$p.adjust), ][1:5, ]
  
  
  pathway <- top5_pathways[1:5, "ID"]
  
 

  pathway_ls[[i]] <-pathway
}

cat(pathway_ls[[1]],sep = ', ')



# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'Rdata_',dir_name, ".RData"))