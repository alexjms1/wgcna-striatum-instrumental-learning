rm(list=ls())
memory.size(max=FALSE)
memory.limit(40000)


setwd("C:/Users/Ajame/Google Drive/FastLMM/WGCNA/One-step-str")


library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()
pdf(file = "heatMap%.03d.pdf", onefile=FALSE, compress=FALSE)
lnames=load(file="Striatum-1step-dataInput.RData")
lnames=load(file="str-1step-networkConstruction-auto.Rdata")

#quantifying module-trait associations

nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

#recalculating MEs

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitCor,file="module_trait_cor.csv")
write.csv(moduleTraitPvalue,file="module_trait_cor_Pvalue.csv")
#correlation of ME with trait in heatmap
#sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

#correlations in heatmap

labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = substr(names(MEs),3,nchar(names(MEs))),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.lab = 0.5,
cex.text = 0.34,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()
#gene significance (GS) and module membership (MM)

#MM
modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
write.csv(geneModuleMembership, file="gene_module_membership.csv")
write.csv(MMPvalue, file="module_membership_p_values.csv")

#GS (must run for each trait of interest)
ActiveBC = as.data.frame(datTraits$ActiveBC) 
names(ActiveBC) = "ActiveBC"
geneTraitSignificance = as.data.frame(cor(datExpr, ActiveBC, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(ActiveBC), sep="");
names(GSPvalue) = paste("p.GS.", names(ActiveBC), sep="");


#correlation between GS and MM


pdf(file = "morefigures%.03d.pdf", onefile=FALSE, compress=FALSE)
intModules = c("brown", "salmon", "turquoise", "green", "greenyellow", "yellow","red", "cyan", "blue")
for (module in intModules)
{
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for total lever presses",
main = paste("Module membership vs. gene significance\n"),
par(c(cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, cex=1.2)) , col = module)
}

dev.off()

pdf(file = "morefiguresTraitME%.03d.pdf", onefile=FALSE, compress=FALSE)
#correlation between trait and module eigengene
intModules = c("brown", "salmon", "turquoise", "green", "greenyellow", "yellow","red", "cyan", "blue")
for (module in intModules)
{
eigen=MEs[, paste("ME",module, sep="")]
par(mfrow = c(1,1));
verboseScatterplot(eigen,datTraits$ActiveBC,
xlab = paste(toupper(substring(module,1,1)), substring(module,2), " eigengene", sep=""),
ylab = "Total lever presses",
main = paste("Module eigengene vs. total lever presses\n"), 
abline=TRUE,
par(c(cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, cex=1.2)) , col = module)
}
dev.off()
#probe-gene annotation (must do for all traits)
annot=read.csv(file="annot.csv")
dim(annot)
probes=names(datExpr)
probes2annot=match(probes,annot[,1])
# The following is the number of probes without annotation:
sum(is.na(probes2annot))

geneInfo0 = data.frame(Probe = probes,
geneSymbol = annot$Gene[probes2annot],
chromosome = annot$chromosome[probes2annot],
chrStart = annot$chrStart[probes2annot],
moduleColor = moduleColors,
geneTraitSignificance,
GSPvalue)

#order modules by significance for trait (must run for each trait of interest)
modOrder = order(-abs(cor(MEs, ActiveBC)));

#add module membership info in chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#order genes in geneInfo variable first by module color, then by gene-trait significance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.ActiveBC));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_ActiveBC.csv")


#Step 4: interfacing network analysis with functional annotation and gene ontology

#output gene lists for use with online software and services
probes = names(datExpr)
probes2annot=match(probes,annot[,1])
allLLIDs = annot$Entrez_Gene_ID[probes2annot];



intModules = c("brown", "salmon", "turquoise", "green", "greenyellow", "yellow","red", "cyan", "purple", "blue")
for (module in intModules)
{
modGenes = (moduleColors==module)
modLLIDs = allLLIDs[modGenes];
fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
write.table(as.data.frame(modLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)
}

#all probes for analysis
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)

#enrichment analysis directly within R
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
row.names = FALSE, col.names = FALSE)

#Step 5: visualizing network of eigengenes

#recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
ActiveBC = as.data.frame(datTraits$ActiveBC);
names(ActiveBC)="ActiveBC"
# Add the trait to existing module eigengenes - considering it as an eigengene
MET = orderMEs(cbind(MEs, ActiveBC))

pdf(file = "dendroHeatmap%.03d.pdf", onefile=FALSE, compress=FALSE)
#plot relationship among eigengenes and trait
#sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.75, xLabelsAngle
= 90)

#split the two graphs
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
collectGarbage()



par(cex = 1.0)
pdf(file = "bandedPlots%.03d.pdf", onefile=FALSE, compress=FALSE)
intModules = c("brown", "salmon", "turquoise", "green", "greenyellow", "yellow","red", "cyan", "purple", "blue")
for (module in intModules)
{
par(cex.main=1.2)
which.module=module
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=1.8)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="Eigengene expression",xlab="Strains")

}

dev.off()

#pdf(file = "pairwisePlots%.03d.pdf", onefile=FALSE, compress=FALSE)
#par(mfrow=c(1,1), cex.main=1.0)
#par(cex=0.8)
#par(mar=c(3.1, 2.1,  2.1, 0.1))
#pair-wise plot
#plotMEpairs(MEs,y=t(as.vector(ActiveBC)))
#plotMEpairs(MEs,y=t(as.vector(ActiveBC)), labels = names(c(ActiveBC,substr(names(MEs),3,nchar(names(MEs)))))
#plotMEpairs(MEs,y=t(as.vector(ActiveBC)), yaxt='n', xaxt='n')

#plotMEpairs(MEs,y =t(as.vector(ActiveBC)), par(cex.labels = .3, cex.cor),yaxt='n', xaxt='n',)

#plotMEpairs(MEs,y=t(as.vector(ActiveBC)), labels = names(c(ActiveBC,substr(names(MEs),3,nchar(names(MEs)))), par(cex.labels = .3),yaxt='n', xaxt='n')
#plotMEpairs(MEs,y=t(as.vector(ActiveBC)), par(cex.labels = .3))
#plotMEpairs(MEs,y =t(as.vector(ActiveBC)), par(cex.labels = .3, cex.cor),yaxt='n', xaxt='n',)
#plotMEpairs(MEs,y=t(as.vector(ActiveBC)), labels = c("ActiveBC",substr(names(MEs),3,nchar(names(MEs))), par(cex.labels = .3))
#plotMEpairs(MEs,y=t(as.vector(ActiveBC)), labels = c("ActiveBC",substr(names(MEs),3,nchar(names(MEs))), par(cex.labels = .3),yaxt='n', xaxt='n')

#par(mfrow=c(1,1), cex=1.0)
#Part 6: Exporting gene network to external visualization software

#Exporting to VisANT

#TOM = TOMsimilarityFromExpr(datExpr, power = 12);
#plotTOM=(1-TOM)^7
#diag(plotTOM) = NA;
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
lnames=load(file="visant-TOM.RData")
annot = read.csv(file = "annot.csv"); 
intModules = c("green", "greenyellow", "yellow","red", "cyan", "brown", "purple")
for (module in intModules)
{
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];

modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0,
probeToGene = data.frame(annot$Probe, annot$Gene))
}

#Step 7: MM, intramodular connectivity, and screening for intramodular hub genes

#measure of MS as average GS
pdf(file = "meansig%.03d.pdf", onefile=FALSE, compress=FALSE)
GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
GeneSignificance=abs(GS1)
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColors)
signif(cor(datTraits,MEs,use="p"),2)
corME=signif(cor(datTraits,MEs,use="p"),2)
dev.off()
par(cex=1)

write.csv(corME,"correlation_ME_traits.csv")

datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
write.csv(datKME,file="MM_values_all_probes.csv")

#relationship between GS and intramodular connectivity
#colorlevels=unique(moduleColors)
#par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
#par(mar = c(4,5,3,1))
#for (i in c(1:length(colorlevels)))
#{
#whichmodule=colorlevels[[i]];
#restrict1 = (moduleColors==whichmodule);
#verboseScatterplot(Alldegrees1$kWithin[restrict1],
#GeneSignificance[restrict1], col=moduleColors[restrict1],
#main=whichmodule,
#xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
#}


#modules of interest --> c("brown", "salmon", "turquoise", "green", "greenyellow", "yellow","red", "cyan", "purple", "blue")


GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.salmon)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.salmon[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_salmon_ActiveBC.csv")


GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.green)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.green[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_green_ActiveBC.csv")


GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.brown)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.brown[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_brown_ActiveBC.csv")

GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.turquoise)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.turquoise[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_turquoise_ActiveBC.csv")

GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.red)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.red[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_red_ActiveBC.csv")

GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.yellow)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.yellow[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_yellow_ActiveBC.csv")

GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.greenyellow)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.greenyellow[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_greenyellow_ActiveBC.csv")

GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.cyan)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.cyan[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_cyan_ActiveBC.csv")

GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.purple)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.purple[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_purple_ActiveBC.csv")


GS1=as.numeric(cor(ActiveBC,datExpr, use="p"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.blue)>.8
table(FilterGenes)
hub_genes=dimnames(data.frame(datExpr))[[2]][FilterGenes]
hub_genes_output=data.frame(hub_genes,GS1[FilterGenes],datKME$MM.blue[FilterGenes],geneInfo0$geneSymbol[FilterGenes])
write.csv(hub_genes_output,file="hub_genes_blue_ActiveBC.csv")
collectGarbage()
dev.off()
