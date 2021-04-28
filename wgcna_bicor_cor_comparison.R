#pooja,singh09@gmail.com
#running biweight 
#load libraries

library("WGCNA")
library("BiocParallel")

#important steps

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# read data

data <- read.csv("../DESeq_NormalizedCounts_vst_filt.txt", header=T, sep="\t")
dim(data)

# format data frame

datExpr0 = as.data.frame(t(data[, -c(1)]))
names(datExpr0) = data$gene
rownames(datExpr0) = names(data)[-c(1)]



#check if genes are good

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


# pick soft threshold

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, corFnc="cor")


#plot soft threshhold and mean connectivity as a function of soft threshold 

pdf(file = "1_soft_threshold_cor.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence using pearson correlation"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity using pearson correlation")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()


sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, corFnc="bicor")


#plot soft threshhold and mean connectivity as a function of soft threshold 

pdf(file = "1_soft_threshold_bicor.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence using biweight correlation"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity using biweight correlation")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()



