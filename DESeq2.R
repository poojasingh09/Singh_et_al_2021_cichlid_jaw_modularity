#pooja.singh09@gmail.com
#Singh et al study on the transcriptional modularity of cichlid jaws

library("DESeq2")
library("BiocParallel")
sessionInfo()
packageVersion("DESeq2")

BPPARAM = MulticoreParam(workers=5)

setwd(".")

# reads files
directory <- ".."
sampleFiles <- grep(".count",list.files(directory),value=TRUE)

# import sample table
sampleTable <- read.table("../paper2_speciesinfo2.txt", header=T)

#start analysis (condition == species; module === jaw type)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition + module)
dds <- dds[ rowSums(counts(dds)) > 100, ]

#run analysis

dds <- DESeq(dds) 

res <- results(dds)
write.table(res, file = "DESeq_allsamples_fullresults.txt", sep = "\t", quote=F)

#save data
save(dds, file="./DESeq2_all.RData")

#plot dispersion plot
png("dispersion_plot.png")
plotDispEsts(dds)
dev.off()

#matrix of normalized counts

rld <- rlog(dds) 
vsd <- varianceStabilizingTransformation(dds)
matrix <- (assay(rld))
write.table(matrix, file = "DESeq_NormalizedCounts_rld.txt", sep = "\t", quote=F)

matrix2 <- (assay(vsd))
write.table(matrix2, file = "DESeq_NormalizedCounts_vsd.txt", sep = "\t", quote=F)

#########################################################################

#plot SD of transformed data using shifted logarithm transformation

library("vsn") 
notAllZero <- (rowSums(counts(dds))>5)


png("DESeq_meanSDplot_normalized_shiftedlog.png") 
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
dev.off()

png("DESeq_meanSDplot_normalized_regularizedlog.png")
meanSdPlot(assay(rld[notAllZero,]))
dev.off()

png("DESeq_meanSDplot_normalized_varianceStabilizingTransformation.png")
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

##########################################################################

#heatmap of normalized counts

library("pheatmap") 
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,] 
df <- as.data.frame(colData(dds)[,("condition")]) 

png("heatmap_20_highest_expr_genes_log2normalized.png")
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
dev.off()

png("heatmap_20_highest_expr_genes_rld.png")
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
dev.off()

png("heatmap_20_highest_expr_genes_vsd.png")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
dev.off()

#########################################################################

#make contrasts

res_OJPJ <- results(dds, contrast=c("module", "OJ", "LPJ"), alpha=0.05)
write.table(res_OJPJ, file = "OJ_LPJ_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)


res_OJPJ2 <- results(dds, contrast=c("module", "OJ", "LPJ"), alpha=0.01)
write.table(res_OJPJ, file = "OJ_LPJ_DESeq_fullresults_0.01.txt", sep = "\t", quote=F)


########################################################################

#summarise results


#summary of results

sum <- summary(res_OJPJ)
write.table(sum, file = "OJPJ_DESeq_summary_0.05.txt", sep = "\t", quote=F)

genesDE <- sum(res_OJPJ$padj < 0.05, na.rm=TRUE)
write.table(genesDE, file = "OJPJ_DESeq_genesDE_0.05.txt", sep = "\t", quote=F)

sum1 <- summary(res_OJPJ2)
write.table(sum1, file = "OJPJ_DESeq_summary_0.01.txt", sep = "\t", quote=F)

genesDE1 <- sum(res_OJPJ2$padj < 0.01, na.rm=TRUE)
write.table(genesDE1, file = "OJPJ_DESeq_genesDE_0.01.txt", sep = "\t", quote=F)

##########################################################################

#plot MA plot

png("DEseq_MAplot_0.05.png")
plotMA(res_LdTd, main="DESeq2_0.05", ylim=c(-2,2))
dev.off()

png("DEseq_MAplot_0.01.png")
plotMA(res_LdTd1, main="DESeq2_0.01", ylim=c(-2,2))
dev.off()

##########################################################################

#heatmaps of sample to samples distances: clustering

png("DEseq_clustering_heatmap_ofsamples_rld.png")
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer") 
sampleDistMatrix <- as.matrix(sampleDists) 
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
colnames(sampleDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=colors)
dev.off()


png("DEseq_clustering_heatmap_ofsamples_vsd.png")
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer") 
sampleDistMatrix <- as.matrix(sampleDists) 
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
colnames(sampleDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=colors)
dev.off()

#######################################################################


#principal component plot of samples

png("PCA_samples_rld.png")
plotPCA(rld, intgroup=c("condition", "module"))
dev.off()

png("PCA_samples_vsd.png")
plotPCA(vsd, intgroup=c("condition", "module"))
dev.off()

#netapp/home/tmp/RtmpAeSB6R/downloaded_packages

### plots for manuscript

library(ggplot2)
vsd <- varianceStabilizingTransformation(dds)
data <- plotPCA(vsd, intgroup=c("condition", "module"), returnData=TRUE) 
percentVar <- round(100 * attr(data, "percentVar")) 
pdf("deseq_modularity_PCA_2021.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=module)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_shape_manual(values=c(17, 16))
dev.off()

svg("deseq_modularity_PCA_2021.svg")
ggplot(data, aes(PC1, PC2, color=condition, shape=module)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_shape_manual(values=c(17, 16))
dev.off()

