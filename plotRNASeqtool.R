source('plotRNASeqFun.R')
rpkm <- read.csv("./data/2015_11_30_rpm_RNAseq.csv", stringsAsFactors = FALSE)
graphTitles <- c("iPS WTC", "CRISPRn No Guide", "CRISPRi No Guide", "CRISPRi GcAMP", "CRISPRi BAG3", "CRISPRi aActinin")


# library('biomaRt')
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- rpkm$X
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),values=genes,mart= ensembl)
# rpkm <- merge(rpkm,G_list,by = 1 ,all = TRUE)

for(i in 2:7){
  xIndex1 = i
  yIndex1 = i + 6
  graph <- pairwiseRNASeqPlot(rpkm,
                                   xIndex = xIndex1,
                                   yIndex =  yIndex1,
                                   plotTitle = graphTitles[i - 1],
                                   xAxisTitle = "- Dox (RPM)",
                                   yAxisTitle = "+ Dox (RPM)",
                                   genesOfInterest = c('Oct4', "Sox2", 'Nanog', 'HERG','ROCK1','BAG3','GSK3B','dCas9-KRABS'), 
                                   geneLabels = c('ENSG00000204531','ENSG00000181449',
                                                  'ENSG00000111704','ENSG00000055118',
                                                  'ENSG00000067900','ENSG00000151929',
                                                  'ENSG00000100744','KRAB-dCAS9-P2A-mCherry'),
                                   geneColors = c('#d13632', '#1d829e', '#ec883a','#96bf33','#8a2aa7',
                                                  '#479e1b', '#503fa9','orange'))
  saveName <- as.character(paste('./figures/',graphTitles[i-1],'_RPM_Cas9.png', sep = ''))
  ggsave(filename = saveName, width=9, height=9)
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

readNames <- read.csv('./data/RNA-seq Samples.csv', header = TRUE)
names <- c("Geneid", "Chr", "Start", "End", "Strand","Length", "WTC - No Dox - Sample A", "CRISPRn no guide- No Dox - Sample A", "CRISPRi no guide- No Dox - Sample A", "CRISPRi GCaMP - No Dox - Sample A", "CRISPRi BAG3 - No Dox - Sample A", "TetO-AlphaActinin - No Dox - Sample A", "WTC - Plus Dox - Sample A", "CRISPRn no guide- Plus Dox - Sample A", "CRISPRi no guide- Plus Dox - Sample A", "CRISPRi GCaMP - Plus Dox - Sample A", "CRISPRi BAG3 - Plus Dox - Sample A", "TetO-AlphaActinin - Plus Dox - Sample A", "WTC - No Dox - Sample B", "CRISPRn no guide- No Dox - Sample B", "CRISPRi no guide- No Dox - Sample B", "CRISPRi GCaMP - No Dox - Sample B", "CRISPRi BAG3 - No Dox - Sample B", "TetO-AlphaActinin - No Dox - Sample B", "WTC - Plus Dox - Sample B", "CRISPRn no guide- Plus Dox - Sample B", "CRISPRi no guide- Plus Dox - Sample B", "CRISPRi GCaMP - Plus Dox - Sample B", "CRISPRi BAG3 - Plus Dox - Sample B", "TetO-AlphaActinin - Plus Dox - Sample B")
colnames(readData) <- names

names1 <- c("Geneid", "Chr", "Start", "End", "Strand","Length","CRISPRi no guide- No Dox - Sample A",
            "CRISPRi no guide- No Dox - Sample B", "CRISPRi no guide- Plus Dox - Sample A",
            "CRISPRi no guide- Plus Dox - Sample B","CRISPRi GCaMP - No Dox - Sample A",
            "CRISPRi GCaMP - No Dox - Sample B","CRISPRi GCaMP - Plus Dox - Sample A",
            "CRISPRi GCaMP - Plus Dox - Sample B","CRISPRi BAG3 - No Dox - Sample A",
            "CRISPRi BAG3 - No Dox - Sample B","CRISPRi BAG3 - Plus Dox - Sample A",
            "CRISPRi BAG3 - Plus Dox - Sample B","CRISPRn no guide- No Dox - Sample A",
            "CRISPRn no guide- No Dox - Sample B","CRISPRn no guide- Plus Dox - Sample A",
            "CRISPRn no guide- Plus Dox - Sample B")
colnames(merged) <- names1
library(plyr)
readCounts <- rbind.fill(readData,merged)
for( i in 7:30){readCounts[,i] <- as.numeric(readCounts[,i])}
readCounts <- readCounts[-64705,]

row.names(readCounts) <- readCounts$Geneid

#Sum replicates (A in col 7-18, B in col 19-30)
for( i in 7:18){
  readCounts[,i] <- rowSums(readCounts[,c(i,i+12)])
}

readCounts <- readCounts[,1:18]

#Generate a table of only reads, and variables for RPKM calculation
readsTable <- readCounts[,7:18]
readsTable[is.na(readsTable)] <- 0
#Normalize data to count (first part of FPKM) and then calculate the rest of the FPKM.
normToCount <- sweep(readsTable,2,colSums(readsTable), '/')
rpm <- (normToCount * 10^9)

#Read the names of the samples, chance the column names, and filter out '0' values.
readNames <- read.csv('./data/RNA-seq Samples.csv', header = TRUE)
colnames(rpkm) <- readNames[1:12,1] 
rpkm[rpkm == 0] <- NA 
#Remove rows with all 'NA' values to tidy the dataset. 
rpm <- rpm[rowSums(is.na(rpm)) != ncol(rpm),]
if(!file.exists("./data/2015_11_30_rpkm_RNAseq.csv")){
  write.csv(rpm, file = "./data/2015_11_30_rpkm_RNAseq.csv")
}

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

