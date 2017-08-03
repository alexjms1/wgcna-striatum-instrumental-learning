library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()
pdf(file = "strMoreStrains%.03d.pdf", width = 15, height = 15, onefile=FALSE, compress=FALSE)
StriatData=read.csv("striatum_expression1.csv")
datExpr0=as.data.frame(t(StriatData[,-c(1)]))
names(datExpr0)=StriatData$Probe_Id
rownames(datExpr0)=names(StriatData)[-c(1)]
gsg = goodSamplesGenes(datExpr0,verbose=3)
gsg$allOK
sampleTree=flashClust(dist(datExpr0),method="average")
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=16,col="red")
clust=cutreeStatic(sampleTree,cutHeight=16,minSize=10)
table(clust)
keepSamples=(clust==1)
datExpr=datExpr0[keepSamples,]
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

#load phenotype data
PhenoData=read.csv("Rates_BC_and_non.csv")

Striat=rownames(datExpr)
traitRows=match(Striat,PhenoData$Strains)
datTraits = PhenoData[traitRows,-1]
rownames(datTraits)=PhenoData[traitRows,1]
collectGarbage()

sampleTree2 = flashClust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits,file="Striatum-1step-dataInput.RData")

#Step 2: Network Construction

lnames=load(file="Striatum-1step-dataInput.RData")
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

#soft-thresholding

powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#power of 12 looks good.  
#network construction
#note single block here - this takes a while and a lot of RAM - consider virtual memory (especially if using an SSD)

bwnet = blockwiseModules(datExpr, maxBlockSize=100000,
power = 12, TOMType = "unsigned", minModulesize=30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE,
saveTOMs = TRUE,
saveTOMFileBase = "strTOM-block",
verbose = 3)

#dendrogram

table(bwnet$colors)
mergedColors=labels2colors(bwnet$colors)
moduleLabels=bwnet$colors
bwmoduleColors=labels2colors(bwnet$colors)
bwLabels=matchLabels(bwnet$colors,moduleLabels)
bwModuleColors=labels2colors(bwLabels)

sizeGrWindow(6,6)
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
"Module colors", main = "Gene dendrogram and module colors in single block",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

MEs=bwnet$MEs
geneTree=bwnet$dendrograms[[1]]
save(MEs, moduleLabels, bwmoduleColors,geneTree,file="str-1step-networkConstruction-auto.RData")