## WGCNA code used for global coexpression network analysis used in Singh et al 2021
## https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-021-01787-9
## pooja.singh09@gmail.com

#load libraries

library("WGCNA")
library("BiocParallel")
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# read data

data <- read.csv("DESeq_NormalizedCounts_vst_filt.txt", header=T, sep="\t")
dim(data)

# format data frame

datExpr0 = as.data.frame(t(data[, -c(1)]))
names(datExpr0) = data$gene
rownames(datExpr0) = names(data)[-c(1)]



# check if genes are good

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


# pick soft threshold

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)


#plot soft threshhold and mean connectivity as a function of soft threshold 

pdf(file = "1_soft_threshold.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()


# pick soft threshold with biweight

sft1 = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, corFnc = "bicor")

#plot soft threshhold and mean connectivity as a function of soft threshold 

pdf(file = "1_soft_threshold_bicor.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft1$fitIndices[,1], sft1$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity")) 
text(sft1$fitIndices[,1], sft1$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()


######################## step by step method ###########################

# calculate adjaceny

softPower = 6;
adjacency = adjacency(datExpr0, power = softPower, type = "signed hybrid");


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM


# Call the hierarchical clustering function

geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)


pdf(file = "6_clustering.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, method="hybrid");
modules <- table(dynamicMods)
write.table(modules, "7_modules.txt", sep="\t")

# Convert numeric lables into colors

dynamicColors = labels2colors(dynamicMods)
modulesC <- table(dynamicColors)
write.table(modulesC, "8_dynamicCmodules.txt", sep="\t")

# Plot the dendrogram and colors underneath

pdf(file = "9_network.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes


MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes

### calculate ME hist and median absolute distribution (2021)


pdf("MEs_6MAD_plots.pdf")

names<-names(MEs)
classes<-sapply(MEs,class)

for(name in names[classes == 'numeric'])
{
    #dev.new()
    mad <- mad(MEs[,name], center=median(MEs[,name]))
    med <- median(MEs[,name])
    hist(MEs[,name],xlim=c(-2,2), xlab=name, main="Histogram of MEs and mean absolute deviation") # subset with [] not $
    abline(v=med, lty=2)
    abline(v=med+(6*mad), lty=2, col="red")
    abline(v=med-(6*mad), lty=2, col="red")
    abline(v=med+(4*mad), lty=2, col="green")
    abline(v=med-(4*mad), lty=2, col="green")
}
dev.off()


###

# Calculate dissimilarity of module eigengenes

MEDiss = 1-cor(MEs);

# Cluster module eigengenes

METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result

pdf(file = "10_module_merge.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")

MEDissThres = 0.25

# Plot the cut line into the dendrogram

abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors

mergedColors = merge$colors;

# Eigengenes of the new merged modules:

mergedMEs = merge$newMEs;


#plot merged results

sizeGrWindow(12, 9)
pdf(file = "11_merged_modules.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

########### correlating to traits ########## make sure you convert categorical variables to binary contrasts https://peterlangfelder.com/2018/11/25/working-with-categorical-variables

#load trait data

traitData <- read.table("traits", header=T)
allTraits = traitData[, -c(2,4,8)];
dim(allTraits)
names(allTraits)


Samples = rownames(datExpr0);
traitRows = match(Samples, allTraits$sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
sizeGrWindow(12, 9)
svg(file = "12_sample_trait_clustering_2021.svg", wi = 9, he = 6)
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits),main = "Sample dendrogram and trait heatmap")
dev.off()


# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot

svg(file = "13_trait_modules_corr_2021.svg")
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = TRUE,cex.text = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

# Save module colors and labels for use in subsequent parts
save(datExpr0, datTraits, MEs, moduleLabels, moduleColors, geneTree, file = "multistep.RData")

#get genes with high correlation to trait

module_darkgreen <- names(datExpr0)[moduleColors=="darkgreen"]
write.table(module_darkgreen, file = "module_darkgreen.txt", sep="\t")

module_black <- names(datExpr0)[moduleColors=="black"]
write.table(module_black, file = "module_black.txt", sep="\t")

module_brown <- names(datExpr0)[moduleColors=="brown"]
write.table(module_brown, file = "module_brown.txt", sep="\t")

module_black <- names(datExpr0)[moduleColors=="black"]
write.table(module_black, file = "module_black.txt", sep="\t", quote=F, row.names=F)

module_purple <- names(datExpr0)[moduleColors=="purple"]
write.table(module_purple, file = "module_purple.txt", sep="\t", quote=F, row.names=F)

#Drerio annottate the modules



# Save module colors and labels for use in subsequent parts
save(datExpr0, datTraits, MEs, moduleLabels, moduleColors, geneTree, file = "all.RData")

#Gene relationship to trait and important modules: Gene Significance and conectivity


#We begin by calculating the intramodular connectivity for each gene. (In network literature, connectivity is often referred to as ”degree”.) The function intramodularConnectivity computes the whole network connectivity kTotal, the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal 

#Intramodular connectivity measures how connected, or co-expressed, a given gene is with respect to the genes of a particular module. The intramodular connectivity may be interpreted as a measure of module membership. 



GS1=as.numeric(cor(module,datExpr0, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, mergedColors, mean, na.rm=T)

pdf("15_GS.pdf",w=25,h=10)
sizeGrWindow(25,10)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,mergedColors)
dev.off()



ADJ1=abs(cor(datExpr0,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, mergedColors)
head(Alldegrees1)


pdf(file = "14_GS_Conn.pdf", w=18, h=12)
colorlevels=unique(mergedColors)
sizeGrWindow(18,12)
par(mfrow=c(4,5))
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (mergedColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
GeneSignificance[restrict1], col=mergedColors[restrict1],
main=whichmodule,
xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
} 
dev.off()

#7.c Generalizing intramodular connectivity for all genes on the array 
#The intramodular connectivity measure is only defined for the genes inside a given module. But in practice it can be very important to measure how connected a given genes is to biologically interesting modules. 
#we define a module eigengene-based connectivity measure for each gene as the correlation between a the gene expression and the module eigengene 
datKME=signedKME(datExpr0, MEs, outputColumnName="MM.")
head(datKME)

# nowWe have a module membership value for each gene in each module. 

# extract hub genes with high GS and MM  for modules of interest, the most signicant genes are also the most connected.
#Our previous analysis has shown that the brown module is an “interesting”module in that its module significance is high. Here we show how to find genes with high gene significance and high intramodular connectivity in the brown module. 

FilterGenes= abs(GS1)> .9 & abs(datKME$MM.blue) > .9
table(FilterGenes)
hubgenes <- dimnames(data.frame(datExpr0))[[2]][FilterGenes]
write.table(hubgenes, "module_blue_hub_0.9.txt", sep="\t")


FilterGenes= abs(GS1)> .95 & abs(datKME$MM.blue) > .95
table(FilterGenes)
hubgenes <- dimnames(data.frame(datExpr0))[[2]][FilterGenes]
write.table(hubgenes, "module_blue_hub_0.95.txt", sep="\t")

####################### plot connectivity ##########################

pdf("16_connectivity_MM.pdf")
sizeGrWindow(12,9)
par(mfrow=c(3,3))
which.color="blue";
restrictGenes=mergedColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                 (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                 col=which.color,
                 xlab="Intramodular Connectivity",
                 ylab="(Module Membership)^6")
which.color="cyan";
restrictGenes=mergedColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                 (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                 col=which.color,
                 xlab="Intramodular Connectivity",
                 ylab="(Module Membership)^6")
which.color="plum2";
restrictGenes=mergedColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                 (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                 col=which.color,
                 xlab="Intramodular Connectivity",
                 ylab="(Module Membership)^6")

which.color="floralwhite";
restrictGenes=mergedColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                 (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                 col=which.color,
                 xlab="Intramodular Connectivity",
                 ylab="(Module Membership)^6")

which.color="saddlebrown";
restrictGenes=mergedColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                 (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                 col=which.color,
                 xlab="Intramodular Connectivity",
                 ylab="(Module Membership)^6")


which.color="grey";
restrictGenes=mergedColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                 (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                 col=which.color,
                 xlab="Intramodular Connectivity",
                 ylab="(Module Membership)^6")
dev.off()






#####################################


#visualization of inter-module relationship

dissTOM = 1-TOMsimilarityFromExpr(datExpr0, power = 6);
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


