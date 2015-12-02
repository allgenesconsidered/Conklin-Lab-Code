rpkm <- read.csv('./data/2015_11_19_rpm_RNAseq1.csv', stringsAsFactors = FALSE)
graphTitles <- c("iPS WTC", "CRISPRn No Guide", "CRISPRi No Guide", "CRISPRi GcAMP", "CRISPRi BAG3", "CRISPRi aActinin")


pairwiseRNASeqPlot <- function(RNAdata, xIndex, yIndex,
                               plotTitle = 'Graph Title',
                               xAxisTitle = 'Condition A (RPKM)',
                               yAxisTitle = 'Condition B (RPKM)',
                               genesOfInterest = c('Oct4', "Sox2", 'Nanog', 'HERG','ROCK1','BAG3','GSK3B'), 
                               geneLabels = c('ENSG00000204531','ENSG00000181449',
                                              'ENSG00000111704','ENSG00000055118',
                                              'ENSG00000067900','ENSG00000151929','ENSG00000100744'),
                               geneColors = c('#d13632', '#1d829e', '#ec883a','#96bf33','#8a2aa7','#479e1b', '503fa9')){
  require(ggplot2)
  require(scales)
  require(grid)
  geneIndex = RNAdata[0,] 
  for( i in 1:length(genesOfInterest)) {
    geneIndex <-rbind(geneIndex,RNAdata[RNAdata[,1] == geneLabels[i],])
  }
  geneIndex <- cbind(geneIndex, colors = geneColors)
  geneIndex$X <- genesOfInterest
  pairwisePlot <- ggplot(data = RNAdata, aes(x = RNAdata[,xIndex], y = RNAdata[,yIndex]) ,
                         environment = environment()) +
    geom_abline(intercept = 0, slope = 1, size = 1, colour = 'grey0') +
    geom_point(size = 2, colour = 'black') +
    geom_point(data = geneIndex, size = 4, colour = geneIndex$colors, aes(x = geneIndex[,xIndex],
                                                                          y = geneIndex[,yIndex])) +  
    geom_text(data = geneIndex, colour = geneIndex$colors, aes(x = geneIndex[,xIndex],
                                                               y = geneIndex[,yIndex] / 1.5, label = geneIndex$X)) +
    scale_x_log10(name = xAxisTitle, breaks = c(10^-2,10^-1,10^0,10^1,10^2,10^3,10^4),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name = yAxisTitle, breaks = c(10^-2,10^-1,10^0,10^1,10^2,10^3,10^4),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() +
    ggtitle(plotTitle) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_text(face="bold",size=14), axis.text.y = element_text(face="bold",size=14),
          axis.title.x = element_text(face="bold",size=16), axis.title.y = element_text(face="bold",size=16),
          plot.title = element_text(size = 20, face="bold"), axis.ticks.length=unit(-0.20, "cm"), 
          axis.ticks.margin=unit(0.5, "cm"),panel.border = element_rect(size = 1.2, colour = 'black'))
  
  print(pairwisePlot)
  return(pairwisePlot)
}

pairwiseRNASeqPlotBlack <- function(RNAdata, xIndex, yIndex,
                               plotTitle = 'Graph Title',
                               xAxisTitle = 'Condition A (RPKM)',
                               yAxisTitle = 'Condition B (RPKM)',
                               genesOfInterest = c('GCaMP6f', 'Oct4', "Sox2", 'Nanog'), 
                               geneLabels = c('GCaMP6f', 'ENSG00000204531','ENSG00000181449','ENSG00000111704'),
                               geneColors = c('forestgreen', 'steelblue', 'firebrick','maroon2') ){
  require(ggplot2)
  require(scales)
  require(grid)
  geneIndex = RNAdata[0,] 
  for( i in 1:length(genesOfInterest)) {
    geneIndex <-rbind(geneIndex,RNAdata[RNAdata[,1] == geneLabels[i],])
  }
  geneIndex <- cbind(geneIndex, colors = geneColors)
  geneIndex$X <- genesOfInterest
  
  pairwisePlot <- ggplot(data = RNAdata, aes(x = RNAdata[,xIndex], y = RNAdata[,yIndex]) ,
                         environment = environment()) +
    geom_abline(intercept = 0, slope = 1, size = 1, colour = 'grey50') +
    geom_point(size = 2, colour = 'white') +
    geom_point(data = geneIndex, size = 4, colour = geneIndex$colors, aes(x = geneIndex[,xIndex],
                                                                          y = geneIndex[,yIndex])) +  
    geom_text(data = geneIndex, colour = geneIndex$colors, aes(x = geneIndex[,xIndex],
                                                               y = geneIndex[,yIndex] / 1.5, label = geneIndex$X)) +
    scale_x_log10(name = xAxisTitle, breaks = c(10^-2,10^-1,10^0,10^1,10^2,10^3,10^4),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name = yAxisTitle, breaks = c(10^-2,10^-1,10^0,10^1,10^2,10^3,10^4),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() +
    ggtitle(plotTitle) +
    theme_bw() +
    theme(rect = element_rect(fill = 'black', colour = 'black'),
          text = element_text(colour = 'white'),
          line = element_line(colour = 'white'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.key = element_rect(linetype = 'blank'),
          legend.title = element_text(face="bold",size=14),
          legend.text = element_text(size=14),
          axis.text.x = element_text(face="bold",size=14), 
          axis.text.y = element_text(face="bold",size=14),
          axis.title.x = element_text(face="bold",size=16), 
          axis.title.y = element_text(face="bold",size=16),
          axis.ticks.length = unit(-0.20, "cm"),
          axis.ticks.margin = unit(0.5, "cm"),
          axis.ticks = element_line(colour = 'white'),
          panel.border = element_rect(size = 1.2, colour = 'white'), 
          plot.background = element_rect(fill = 'black',colour = 'black'),
          panel.background = element_rect(fill = 'black'),
          plot.title = element_text(size = 20, face="bold"))
  
  print(pairwisePlot)
  return(pairwisePlot)
}