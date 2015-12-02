####################################
####################################
###   ______  ___  _____  _____  ###
###  | ___ \/ _ \|  __ \|____ |  ###
###  | |_/ / /_\ \ |  \/    / /  ###
###  | ___ \  _  | | __     \ \  ###
##   | |_/ / | | | |_\ \.___/ /  ###
###  \____/\_| |_/\____/\____/   ###
####################################
####################################

library('DESeq')
readCounts = read.csv("./data/20151111_combined_featcounts.txt", sep = '\t', header = TRUE)
row.names(readCounts) <- readCounts[,1]
readCounts <- readCounts[,7:30]
names <- read.csv("data/2015_11_13_rpkm_RNAseq1.csv", header = TRUE)
colnames(readCounts) <- colnames(names)[2:25]


#DESeq Analysis

#Step 1, set up data to be analyzed
readCountsBag3 = readCounts[,c(5,17,11,23)] #Only Bag3
datDesign <- factor(c("untreated","untreated","treated","treated")) #Add metadata (treated = Dox)
cds <- newCountDataSet( readCountsBag3, datDesign)
#Step 2, Normalize to Size Factor and estimate variance
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
#Step 3, Compare 2 experemental conditions using a normal binomial test
res = nbinomTest(cds, "untreated", "treated")


#Plot Volcano
#Add a -log10(p-value) column
res$negLog10p = -(log10(res$pval))
#Remove NA values
res = res[complete.cases(res),]
res = res[abs(res$log2FoldChange) < 20,]

library('biomaRt')
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- res$id
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),values=genes,mart= ensembl)
res <- merge(res,G_list,by = 1 ,all = TRUE)


for(i in 1:length(res$external_gene_name)){
  if(!is.na(res[i,10])){
    res[i,1] = res[i,10]
  }
}



write.csv(res, file = './data/2015_11_19_stats_RNAseq1_BAG3.csv')
resBAG3 <- read.csv('./data/2015_11_19_stats_RNAseq1_BAG3.csv')
resBAG3 <- resBAG3[,2:length(resBAG3)]

#change the colr of genes where pval <= 0.005
resBAG3$sig = ifelse(resBAG3$pval <= 0.05 & abs(resBAG3$log2FoldChange) >= 2.0, 'Y','N')
mycolours <- c("Y" = "firebrick", "N" = "white")
sigsub <- resBAG3[resBAG3$sig == 'Y',]

#Plot!
library(ggplot2)
require(scales)
require(grid)

ggplot() +
  geom_point(data = resBAG3, alpha = 0.6, aes(x = log2FoldChange, y = negLog10p, text = id, colour = sig)) +
  scale_color_manual(values = mycolours) +
  geom_text(data = sigsub, size = 3, colour = 'white', aes(x = log2FoldChange, y = negLog10p, label=id, vjust = 2)) +
  scale_x_continuous('log2(Fold Change)', limit = c(-9,9), breaks = c(-8,-6,-4,-2,2,4,6,8)) +
  scale_y_log10("-log10(p-value)",limits = c(0.1,150)) +
  annotation_logticks(sides = 'l', colour = 'white') +
  ggtitle('CRISPRi BAG3')+ 
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(colour = 'white'),
        rect = element_rect(fill = 'black', colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=16),
        axis.title.y = element_text(face="bold",size=16),
        plot.title = element_text(size = 20, face="bold"),
        axis.ticks.length=unit(-0.20, "cm"),
        axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = 'white'),
        panel.border = element_rect(size = 1.2, colour = 'white'),
        plot.background = element_rect(fill = 'black',colour = 'black'),
        panel.background = element_rect(fill = 'black'))

ggplot() + 
  geom_bar(data = sigBAG3, aes(x = log2FoldChange, y = id)) 
library('biomaRt')
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- res$id
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),values=genes,mart= ensembl)
res <- merge(res,G_list,by = 1 ,all = TRUE)

ggsave(filename = 'BAG3_Volcano_black.png', width=7.3, height=6, dpi = 720)