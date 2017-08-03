
rm(list=ls())

memory.size(max=FALSE)
memory.limit(40000)

setwd("C:/Users/ajame/Google Drive/FastLMM/WGCNA/One-step-str")


library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

pdf(file = "strOneStepInput%.03d.pdf", onefile=FALSE, compress=FALSE)

StriatData=read.csv("striatum_expression1.csv")
datExpr0=as.data.frame(t(StriatData[,-c(1)]))
names(datExpr0)=StriatData$Probe_Id
rownames(datExpr0)=names(StriatData)[-c(1)]
gsg = goodSamplesGenes(datExpr0,verbose=3)
gsg$allOK
sampleTree=flashClust(dist(datExpr0),method="average")
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
dim(PhenoData)
names(PhenoData)

Striat=rownames(datExpr)
traitRows=match(Striat,PhenoData$Strains)
datTraits = PhenoData[traitRows, -1]
rownames(datTraits)=PhenoData[traitRows, 1]
collectGarbage()

sampleTree2 = flashClust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits,file="Striatum-1step-dataInput.RData")
dev.off()
collectGarbage()
rm(list=ls())
