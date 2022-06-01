#!/usr/bin/env Rscript
##################################################################################################
##################################################################################################
##################################################################################################

#Load required packages for script
library("DiffBind")
library("tidyverse")
library("ggplot2")
library("ggfortify")
library("DESeq2")
library("ComplexHeatmap")
library("circlize")
library("gplots")
library("dplyr")

#Set up command-line argument input
args = commandArgs(trailingOnly=TRUE)

#Set static Variables
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1234)
RUNNER <- args[1]
INPUT_DIR <- args[2]
OUTPUT_DIR <- args[3]
FIGURE_DIR <- args[4]

#Function for reading in motif files for RUN2
ReadMotifFile = function(motiffile, col.prefix) {
  
  df = read.table(file=motiffile, head = T, sep = ",")
  #df = df[, c(1,3,5,7,9, 10:14)]
  #colnames(df) = c(paste0(col.prefix, ".motif"), paste0(col.prefix, ".pval"), paste0(col.prefix, ".fdr"), paste0(col.prefix, ".fg"), paste0(col.prefix, ".bg"), paste0(col.prefix, ".FC"), paste0(col.prefix, ".logFC"), paste0(col.prefix, ".Escore"), paste0(col.prefix, ".clusterId"), "finalGene")
  df = df[, c(3,5,7,9, 10:14)]
  colnames(df) = c(paste0(col.prefix, ".pval"), paste0(col.prefix, ".fdr"), paste0(col.prefix, ".fg"), paste0(col.prefix, ".bg"), paste0(col.prefix, ".FC"), paste0(col.prefix, ".logFC"), paste0(col.prefix, ".Escore"), paste0(col.prefix, ".clusterId"), "finalGene")
  
  df =  df %>% 
    group_by(finalGene) %>% 
    summarise(across(everything(), list(min)))
  
  return(df)
}

##################################################################################################
##################################################################################################
##################################################################################################
if (RUNNER == 'RUN1'){
  #Read in peak and reads sample sheet formatted following DiffBind manual
  sample_sheet = paste(INPUT_DIR, 'sample_sheet.csv', sep = '/')
  samples <- read.csv(sample_sheet)


  #Generate DBA object from peak files and calculate the consensus peak and differential peak sets
  DBdata <- dba(sampleSheet=samples)
  DBdata <- dba.count(DBdata)
  DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, minMembers = 2)
  DBdata <- dba.analyze(DBdata, method = DBA_DESEQ2) #Can also be changed to EDGER


  #Extract the concensus peakset from the DBA object and generate an output file
  consensusPeaks = dba.peakset(DBdata, peak.format = "bed", consensus = T, bRetrieve = T)
  consensusPeaksDF = data.frame(chromosome=seqnames(consensusPeaks), starts=start(consensusPeaks), ends=end(consensusPeaks),  names=paste0("peak_", names(consensusPeaks)), scores=c(rep("0", length(consensusPeaks))), strand=strand(consensusPeaks), consensusPeaks)
  write.table(consensusPeaksDF, file = paste0(OUTPUT_DIR, "consensusPeaks_rpkm.txt"), sep = "\t", quote = F, row.names = F)


  #Extract the differential peakset from the DBA object and generate output files
  filename = paste0(output_dir, DEcounts[i,1], "-", DEcounts[i,3]) 
  report = dba.report(DBdata, contrast = i, method=DBA_DESEQ2, )
  df1 <- data.frame(seqnames=seqnames(report), starts=start(report)-1, ends=end(report), names=paste0("peak_", names(report)), scores=c(rep("0", length(report))), strand=c(rep(".", length(report))), fold=elementMetadata(report)$Fold)


  #Write pairwise differential peaksets from DBA object and generate bed files for each
  write.table(df1[df1$fold > 0,1:6], file = paste0(filename, "_DE-up_DEseq2.bed"), sep = "\t", quote = F, row.names = F)
  write.table(df1[df1$fold < 0,1:6], file = paste0(filename, "_DE-down_DEseq2.bed"), sep = "\t", quote = F, row.names = F)
  write.table(df1[,1:6], file = paste0(filename, "_DE-all_Deseq2.bed"), sep = "\t", quote = F, row.names = F)


  #Generate master differential peak data table and write it to file
  df2 = data.frame(report)
  df2$start = ((df2$start)-1)
  finaldf = data.frame(df1[,1:6], df2[,c(4,6:11)])
  write.table(finaldf, file = paste0(filename, "_DE-all_Deseq2.csv"), sep = ",", quote = F, row.names = F) 


  #Load in consensus peak files
  Input <- read.table(paste(OUTPUT_DIR,'consensusPeaks_rpkm.txt', sep = '/'), header = T)


  #Make all entries integers and omit rows containing NAs
  for(j in 1:ncol(Input)) { 
    Input[,j] <- as.integer(Input[,j]) 
    }
  Nfur <- na.omit(Input)


  #Normalize peak data and generate principle components
  Nfur <- vst(as.matrix(Nfur))
  Nfur <- t(Nfur)
  test1 <- prcomp(Nfur, scale.=T) 


  #Reattach metadata to peak set
  sNfur <- c('Pre-Dia', 'Pre-Dia', 'Dia(1M)', 'Dia(1M)')
  dNfur <- c('Dev', 'Dev', 'Dia', 'Dia')
  Nfur <- cbind(sNfur, dNfur, Nfur)


  #Plot PCA using peaks from samples 
  pdf(paste(FIGURE_DIR,'Figure2B.pdf', sep = '/'))
  autoplot(test1, x=1, y=2, data=Nfur, size=10, colour='dNfur', shape='nNfur') +
    scale_colour_manual(values=c('orange', 'light blue')) +
    scale_shape_manual(values=c(19)) +
    labs(title= 'N. furzeri PCA')
  dev.off()
}

##################################################################################################
##################################################################################################
##################################################################################################

if (RUNNER == 'RUN2'){

  
}



                  