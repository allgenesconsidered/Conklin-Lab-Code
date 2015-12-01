source('plotRNASeqFun.R')
rpkm <- read.csv('./data/2015_11_19_rpm_RNAseq1.csv')
graphTitles <- c("iPS WTC", "CRISPRn No Guide", "CRISPRi No Guide", "CRISPRi GcAMP", "CRISPRi BAG3", "CRISPRi aActinin")


# library('biomaRt')
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- rpkm$X
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),values=genes,mart= ensembl)
# rpkm <- merge(rpkm,G_list,by = 1 ,all = TRUE)

for(i in c(7)){
  xIndex1 = i
  yIndex1 = i + 6
  graph <- pairwiseRNASeqPlot(rpkm,
                                   xIndex = xIndex1,
                                   yIndex =  yIndex1,
                                   plotTitle = graphTitles[i - 1],
                                   xAxisTitle = "- Dox (RPM)",
                                   yAxisTitle = "+ Dox (RPM)",
                                   genesOfInterest = c('aActinin', 'Oct4', "Sox2", 'Nanog'), 
                                   geneLabels = c('ENSG00000151929', 'ENSG00000204531','ENSG00000181449','ENSG00000111704')),
                                   #geneColors = c('#C98A1F', 'steelblue', 'firebrick','maroon2'))
  #saveName <- as.character(paste('./figures/',graphTitles[i-1],'_RPM_Black.png', sep = ''))
  #ggsave(filename = saveName, width=9, height=9)
}

CRn <- read.csv('./data/20151116_crisprNfeats.txt', sep = '\t', header = FALSE)
CRn <- data.frame(lapply(CRn, as.character), stringsAsFactors=FALSE)
colnames(CRn) <- CRn[2,]
CRn <- CRn[3:5,]

CRi <- read.csv('./data/20151116_crisprifeats.txt', sep = '\t', header = FALSE)
CRi <- data.frame(lapply(CRi, as.character), stringsAsFactors=FALSE)
colnames(CRi) <- CRi[2,]
CRi <- CRi[3:5,]


merged = merge(CRi,CRn, all = TRUE)
readData <- read.csv('./data/20151111_combined_featcounts.txt', sep = '\t', header = TRUE)

names <- c(paste(shQuote(readNames$Sample.Name..Technical.Rep.A., type="cmd"), collapse=", "))





###################################################
##################################################
##################################################


#PCA and Heatmap in DESeq

library('DESeq')
readCounts = read.csv("./data/20151111_combined_featcounts.txt", sep = '\t', header = TRUE)
row.names(readCounts) <- readCounts[,1]
readCounts <- readCounts[,7:30]
names <- read.csv("data/2015_11_13_rpkm_RNAseq1.csv", header = TRUE)
colnames(readCounts) <- colnames(names)[2:25]
readDesign <- data.frame(row.names = colnames(readCounts),
                         condition = c("- Dox","- Dox","- Dox","- Dox","- Dox","- Dox",
                                       "+ Dox","+ Dox","+ Dox","+ Dox","+ Dox","+ Dox",
                                       "- Dox","- Dox","- Dox","- Dox","- Dox","- Dox",
                                       "+ Dox","+ Dox","+ Dox","+ Dox","+ Dox","+ Dox"),
                         cellLine = c('WTC','CRISPRn','CRISPRi','CRISPRi GcAMP', 'CRISPRi BAG3', 'TetO aActinin',
                                      'WTC','CRISPRn','CRISPRi','CRISPRi GcAMP', 'CRISPRi BAG3', 'TetO aActinin',
                                      'WTC','CRISPRn','CRISPRi','CRISPRi GcAMP', 'CRISPRi BAG3', 'TetO aActinin',
                                      'WTC','CRISPRn','CRISPRi','CRISPRi GcAMP', 'CRISPRi BAG3', 'TetO aActinin'))
cdsFull <- newCountDataSet(readCounts, readDesign)
cdsFull = estimateSizeFactors(cdsFull)
cdsFullBlind = estimateDispersions(cdsFull, method = 'blind')
vsdFull = varianceStabilizingTransformation( cdsFullBlind)

library(RColorBrewer)
library(gplots)
select = order(rowMeans(counts(cdsFull)), decreasing = TRUE)[1:60]
# genes <- c('ENSG00000151929', 'GCaMP6f', 'ACTN2mKate')
# genes <-as.numeric(match(genes,row.names(readCounts)))
# for(i in 1:length(genes)){select[[length(select)+1]] = genes[i]}
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(exprs(vsdFull)[select,], col = hmcol,
          scale="row", key=T, keysize=1.2,margins = c(12,4),
          density.info="none", trace="none",cexCol=0.9, labRow=NA)

