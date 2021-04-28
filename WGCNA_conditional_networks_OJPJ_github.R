## WGCNA code used for conditional OJ and PJ coexpression network analysis used in Singh et al 2021
## https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-021-01787-9
## pooja.singh09@gmail.com
## you can choose whether to use the 1 step method or the step-by step method for making coexpression networks


#load libraries

library("WGCNA")
library("BiocParallel")

#important steps

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

#read data

ojdata <- read.csv("DESeq_NormalizedCounts_vst_OJ_filt.txt", header=T, sep="\t")
pjdata <- read.csv("DESeq_NormalizedCounts_vst_LPJ_filt.txt", header=T, sep="\t")

dim(ojdata)
dim(pjdata)

names(ojdata)
names(pjdata)

datExprOJ = as.data.frame(t(ojdata[, -c(1)]))
colnames(datExprOJ) = ojdata$gene
rownames(datExprOJ) = names(ojdata)[-c(1)]

datExprPJ = as.data.frame(t(pjdata[, -c(1)]))
colnames(datExprPJ) = pjdata$gene
rownames(datExprPJ) = names(pjdata)[-c(1)]


#make data comparable

OJ2PJ = match(colnames(datExprOJ), colnames(datExprPJ));
all.equal(colnames(datExprOJ), colnames(datExprPJ))


####################### OJ step by step method ###########################

# calculate adjaceny

softPower = 18;
adjacencyOJ = adjacency(datExprOJ, power = softPower, type = "signed hybrid");


# Turn adjacency into topological overlap
TOMOJ = TOMsimilarity(adjacencyOJ, TOMType = "signed");
dissTOMOJ = 1-TOMOJ


# Call the hierarchical clustering function

geneTreeOJ = hclust(as.dist(dissTOMOJ), method = "average");

# Plot the resulting clustering tree (dendrogram)


pdf(file = "6_clustering_OJ.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plot(geneTreeOJ, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;

# Module identification using dynamic tree cut

dynamicModsOJ = cutreeDynamic(dendro = geneTreeOJ, distM = dissTOMOJ,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, method="hybrid");
modulesOJ <- table(dynamicModsOJ)
write.table(modulesOJ, "7_modules_OJ.txt", sep="\t")

# Convert numeric lables into colors

dynamicColorsOJ = labels2colors(dynamicModsOJ)
modulesCOJ <- table(dynamicColorsOJ)
write.table(modulesCOJ, "8_dynamicCmodules_OJ.txt", sep="\t")

# Plot the dendrogram and colors underneath

pdf(file = "9_network_OJ.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plotDendroAndColors(geneTreeOJ, dynamicColorsOJ, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes


MEListOJ = moduleEigengenes(datExprOJ, colors = dynamicColorsOJ)
MEsOJ = MEListOJ$eigengenes

# Calculate dissimilarity of module eigengenes

MEDissOJ = 1-cor(MEsOJ);

# Cluster module eigengenes

METreeOJ = hclust(as.dist(MEDissOJ), method = "average");

# Plot the result

pdf(file = "10_module_merge_OJ.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
plot(METreeOJ, main = "Clustering of module eigengenes",xlab = "", sub = "")

MEDissThresOJ = 0.45

# Plot the cut line into the dendrogram

abline(h=MEDissThresOJ, col = "red")
dev.off()

# Call an automatic merging function

mergeOJ = mergeCloseModules(datExprOJ, dynamicColorsOJ, cutHeight = MEDissThresOJ, verbose = 3)

# The merged module colors

mergedColorsOJ = mergeOJ$colors;

# Eigengenes of the new merged modules:

mergedMEsOJ = mergeOJ$newMEs;

#plot merged results


sizeGrWindow(12, 9)
pdf(file = "11_merged_modules_OJ.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeOJ, cbind(dynamicColorsOJ, mergedColorsOJ),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColorsOJ = mergedColorsOJ
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabelsOJ = match(moduleColorsOJ, colorOrder)-1;
MEsOJ = mergedMEsOJ;
#########

####################### PJ step by step method ###########################

# calculate adjaceny

softPower = 18;
adjacencyPJ = adjacency(datExprPJ, power = softPower, type = "signed hybrid");


# Turn adjacency into topological overlap
TOMPJ = TOMsimilarity(adjacencyPJ, TOMType = "signed");
dissTOMPJ = 1-TOMPJ


# Call the hierarchical clustering function

geneTreePJ = hclust(as.dist(dissTOMPJ), method = "average");

# Plot the resulting clustering tree (dendrogram)


pdf(file = "6_clustering_PJ.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plot(geneTreePJ, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;

# Module identification using dynamic tree cut

dynamicModsPJ = cutreeDynamic(dendro = geneTreePJ, distM = dissTOMPJ,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, method="hybrid");
modulesPJ <- table(dynamicModsPJ)
write.table(modulesPJ, "7_modules_PJ.txt", sep="\t")

# Convert numeric lables into colors

dynamicColorsPJ = labels2colors(dynamicModsPJ)
modulesCPJ <- table(dynamicColorsPJ)
write.table(modulesCPJ, "8_dynamicCmodules_PJ.txt", sep="\t")

# Plot the dendrogram and colors underneath

pdf(file = "9_network_PJ.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plotDendroAndColors(geneTreePJ, dynamicColorsPJ, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes


MEListPJ = moduleEigengenes(datExprPJ, colors = dynamicColorsPJ)
MEsPJ = MEListPJ$eigengenes

# Calculate dissimilarity of module eigengenes

MEDissPJ = 1-cor(MEsPJ);

# Cluster module eigengenes

METreePJ = hclust(as.dist(MEDissPJ), method = "average");

# Plot the result

pdf(file = "10_module_merge_PJ.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
plot(METreePJ, main = "Clustering of module eigengenes",xlab = "", sub = "")

MEDissThresPJ = 0.45

# Plot the cut line into the dendrogram

abline(h=MEDissThresPJ, col = "red")
dev.off()

# Call an automatic merging function

mergePJ = mergeCloseModules(datExprPJ, dynamicColorsPJ, cutHeight = MEDissThresPJ, verbose = 3)

# The merged module colors

mergedColorsPJ = mergePJ$colors;

# Eigengenes of the new merged modules:

mergedMEsPJ = mergePJ$newMEs;

#plot merged results


sizeGrWindow(12, 9)
pdf(file = "11_merged_modules_PJ.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreePJ, cbind(dynamicColorsPJ, mergedColorsPJ),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColorsPJ = mergedColorsPJ
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabelsPJ = match(moduleColorsPJ, colorOrder)-1;
MEsPJ = mergedMEsPJ;


#########1 step preservation ########

#oj in pj

pdf("OJ_PJ_16_signedhybrid_preservation.pdf",w=11,h=7)
sizeGrWindow(11, 7)
# Set up appropriate screen sectioning
layout(matrix(c(1:4), 4, 1), heights = rep(c(0.8, 0.2), 2));
# Plot the female dendrogram
plotDendroAndColors(geneTreeOJ, moduleColorsOJ,"OJ modules", main = "OJ gene dendrogram and OJ module colors",dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4,cex.main = 2, cex.lab = 1.4, cex.axis = 1.4,addGuide = TRUE);
# Plot the male dendrogram with female module colors
plotDendroAndColors(geneTreePJ, moduleColorsOJ,"OJ modules", main = "LPJ gene dendrogram and OJmodule colors",dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4, cex.main = 2, cex.lab = 1.4, cex.axis = 1.4,addGuide = TRUE);
dev.off()

#pj in oj
pdf("PJ_OJ_16_signedhybrid_preservation.pdf",w=11,h=7)
sizeGrWindow(11, 7)
# Set up appropriate screen sectioning
layout(matrix(c(1:4), 4, 1), heights = rep(c(0.8, 0.2), 2));
# Plot the female dendrogram
plotDendroAndColors(geneTreePJ, moduleColorsPJ,"LPJ modules", main = "LPJ gene dendrogram and LPJ module colors",dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4,cex.main = 2, cex.lab = 1.4, cex.axis = 1.4,addGuide = TRUE);
# Plot the male dendrogram with female module colors
plotDendroAndColors(geneTreeOJ, moduleColorsPJ,"LPJ modules", main = "OJ gene dendrogram and LPJmodule colors",dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4, cex.main = 2, cex.lab = 1.4, cex.axis = 1.4,addGuide = TRUE);
dev.off()

###module OJ presevation in PJ plots


setLabels = c("OJ", "PJ");
multiExpr = list(OJ = list(data = datExprOJ), PJ = list(data = datExprPJ));
multiColor = list(OJ = moduleColorsOJ, PJ = moduleColorsPJ);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
nPermutations = 200,
randomSeed = 1,
verbose = 3)
} );
# Save the results
save(mp, file = "OJ_PJ_preservation.RData");

load(file = "OJ_PJ_preservation.RData");
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 7);
pdf("OJ_PJ_preservation.pdf", wi=10, h=7)
par(mfrow = c(1,2))

for (p in 1:2)
{
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
if (min > -max/10) min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
} else
ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
main = mains[p],
cex = 2.2,
ylab = mains[p], xlab = "Module size", log = "x",
ylim = ylim,
xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.2, cex.axis=0.8)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1.2, offs = 0.08);
# For Zsummary, add threshold lines
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}
# If plotting into a file, close it
dev.off();

##########
###module PJ presevation in OJ plots


setLabels = c("PJ", "OJ");
multiExpr = list(PJ = list(data = datExprPJ), OJ = list(data = datExprOJ));
multiColor = list(PJ = moduleColorsPJ, OJ = moduleColorsOJ);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
nPermutations = 200,
randomSeed = 1,
verbose = 3)
} );
# Save the results
save(mp, file = "PJ_OJ_preservation.RData");


load(file = "PJ_OJ_preservation.RData");
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results

modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 7);
pdf("PJ_OJ_preservation.pdf", wi=10, h=7)
par(mfrow = c(1,2))

for (p in 1:2)
{
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
if (min > -max/10) min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
} else
ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
main = mains[p],
cex = 2.2,
ylab = mains[p], xlab = "Module size", log = "x",
ylim = ylim,
xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.2,cex.axis=0.8)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
# For Zsummary, add threshold lines
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}
# If plotting into a file, close it
dev.off();
########

