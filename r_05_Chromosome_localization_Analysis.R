rm(list = ls())
gc()
library(biomaRt)
library(dplyr)
library(RCircos)

setwd('/data/home/jiangminghui/KC/Project-0035/05_Chromosome_localization_Analysis/')

Overlapgenes <-readRDS('../03_Overlap_genes_Enrichment/01_Overlapgenes.rds')


# Load the necessary packages

# Load the cytoband data (example data, replace with actual file if needed)
data(UCSC.HG19.Human.CytoBandIdeogram)
cytoBandIdeogram <- UCSC.HG19.Human.CytoBandIdeogram

# Initialize RCircos core components with the cytoband data
RCircos.Set.Core.Components(cyto.info = cytoBandIdeogram, chr.exclude = NULL, tracks.inside = 10, tracks.outside = 0)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取基因的位置
gene_positions <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                        filters = "hgnc_symbol",
                        values = Overlapgenes,
                        mart = ensembl)
gene_positions$chromosome_name <- paste0("chr", gene_positions$chromosome_name)

colnames(gene_positions) <- c("Gene", "Chromosome", "chromStart", "chromEnd")
gene_positions <-gene_positions[,c('Chromosome', 'chromStart', 'chromEnd'  , 'Gene')]

saveRDS(gene_positions,'./01_gene_positions.rds')

write.csv(gene_positions,file = './01_gene_positions.csv')
gene_positions<-readRDS('./01_gene_positions.rds')

# Plot the chromosome ideogram
par(mai=c(0.25, 0.25, 0.25, 0.25))
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot()

# Plot the gene locations
RCircos.Gene.Connector.Plot(gene_positions, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(gene_positions, name.col = 4, track.num = 2, side = "in")

# Save the plot to a file
pdf("./01_gene_locations_circos_plot.pdf",width =5,height = 5.5,onefile = F)
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot()

# Plot the gene locations
RCircos.Gene.Connector.Plot(gene_positions, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(gene_positions, name.col = 4, track.num = 2, side = "in")
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(gene_positions, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(gene_positions, name.col = 4, track.num = 2, side = "in")
dev.off()


png("./01_gene_locations_circos_plot.png",width = 5,height =5.5,units = 'in',res = 300)
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot()

# Plot the gene locations
RCircos.Gene.Connector.Plot(gene_positions, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(gene_positions, name.col = 4, track.num = 2, side = "in")
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(gene_positions, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(gene_positions, name.col = 4, track.num = 2, side = "in")
dev.off()
