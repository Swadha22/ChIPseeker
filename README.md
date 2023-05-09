
# ChIPseeker to annotate the peak with WB genes 
      ## installed ChIPseeker, verion 1.26.2


### 1. Make sure to install all the packages from Bioconductor (see the very bottom of the code)
      rm(list=ls())
      library(ChIPseeker)
      library(TxDb.Celegans.UCSC.ce11.refGene)
      library(TxDb.Celegans.UCSC.ce11.ensGene)
      library(clusterProfiler)
      library(GenomicRanges)
      library(magrittr)
      library(tidyverse)
      library(ReactomePA)


  rm(list=ls()) #removes all objects
  setwd("/Users/swadh/Desktop/chipseeker/to_graduate/chipseeker_V2/")
  dir()
  
 
### 2 read the output Excel files (csv) from MACS2
    ### note: make sure to modify the column names in Excel (csv) files

###  R does not like space 

  S39.HTAS1.C <-read.csv("HTAS-1.C.S39.csv")
  S39.HTZ1 <-read.csv("HTZ-1.S39.csv")
  S41.HTAS1.C <-read.csv("HTAS-1.C.S41.csv")
  S41.HTAS1.N <-read.csv("HTAS-1.N.S41.csv")
  S41.HTZ1 <-read.csv("HTZ-1.S41.csv")
  All_HTZ1 <- read.csv("All_HTZ1.csv")
  All_HTAS1 <- read.csv("All_HTAS.1.csv")
 
  head( S39.HTAS1.C)

  
 ### analysis of each column. Let's check Peak Length
 
  summary(S39.HTAS1.C)
  hist( S39.HTAS1.C$length, border="black", col="lightblue",
       xlim=c(0,3000), breaks=100)
  
  summary(S39.HTZ1)
  hist( S39.HTZ1$length, border="black", col="lightblue",
        xlim=c(0,3000), breaks=100)
  
  summary(S41.HTZ1)
  hist( S41.HTZ1$length, border="black", col="lightblue",
        xlim=c(0,3000), breaks=100)
  
  summary(S41.HTAS1.C)
  hist( S41.HTAS1.C$length, border="black", col="lightblue",
        xlim=c(0,3000), breaks=100)
  
  summary(S41.HTAS1.N)
  hist( S41.HTAS1.N$length, border="black", col="lightblue",
        xlim=c(0,3000), breaks=100)
  
  summary(All_HTAS1)
  hist( All_HTAS1$length, border="black", col="lightblue",
        xlim=c(0,3000), breaks=100)
  
  summary(All_HTZ1)
  hist( All_HTZ1$length, border="black", col="lightblue",
        xlim=c(0,3000), breaks=100)
  
  
### analysis of each column. Let's check Lo10.qvalue
  summary(S39.HTAS1.C$LOG10.qvalue)
  hist(S39.HTAS1.C$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  summary(S41.HTAS1.C$LOG10.qvalue)
  hist(S41.HTAS1.C$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  summary(S41.HTAS1.N$LOG10.qvalue)
  hist(S41.HTAS1.N$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  summary(S41.HTZ1$LOG10.qvalue)
  hist(S41.HTZ1$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  summary(S39.HTZ1$LOG10.qvalue)
  hist(S39.HTZ1$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  
  summary(All_HTAS1$LOG10.qvalue)
  hist(All_HTAS1$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  
  summary(All_HTZ1$LOG10.qvalue)
  hist(All_HTZ1$qvalue, border="black", col="lightblue",
       xlim=c(0,20), breaks=100)
  
  #analysis of each column. Let's check "fold_enrichment
  summary(S39.HTAS1.C$fold_enrichment)
  hist(S39.HTAS1.C$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  summary(S39.HTZ1$fold_enrichment)
  hist(S39.HTZ1$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  summary(S41.HTAS1.C$fold_enrichment)
  hist(S41.HTAS1.C$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  summary(S41.HTAS1.N$fold_enrichment)
  hist(S41.HTAS1.N$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  summary(S41.HTZ1$fold_enrichment)
  hist(S41.HTZ1$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  summary(All_HTAS1$fold_enrichment)
  hist(All_HTAS1$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  summary(All_HTZ1$fold_enrichment)
  hist(All_HTZ1$fold_enrichment, border="black", col="lightblue",
       xlim=c(0,6), breaks=10)
  
  
### making gr range object for auto
  S39.HTAS1.C <- with(S39.HTAS1.C, 
                    GRanges(chr,
                            IRanges(start,end),
                    length=length,
                    log10q = qvalue
                    ))

  S39.HTZ1 <- with(S39.HTZ1, 
                           GRanges(chr,
                                   IRanges(start,end, names=name),
                                   pileup=pileup,
                                   FoldEnrichment= fold_enrichment,
                                   log10q = qvalue
                           ))
 
  
  S41.HTAS1.C <- with(S41.HTAS1.C, 
                           GRanges(chr,
                                   IRanges(start,end, names=name),
                                   pileup=pileup,
                                   FoldEnrichment= fold_enrichment,
                                   log10q = qvalue
                           ))

  
  S41.HTAS1.N <- with(S41.HTAS1.N, 
                           GRanges(chr,
                                   IRanges(start,end, names=name),
                                   pileup=pileup,
                                   FoldEnrichment= fold_enrichment,
                                   log10q = qvalue
                           ))
  
  S41.HTZ1 <- with(S41.HTZ1, 
                        GRanges(chr,
                                IRanges(start,end, names=name),
                                pileup=pileup,
                                FoldEnrichment= fold_enrichment,
                                log10q = qvalue
                        ))
  
  
  All_HTAS1 <- with(All_HTAS1, 
                   GRanges(chr,
                           IRanges(start,end, names=name),
                           pileup=pileup,
                           FoldEnrichment= fold_enrichment,
                           log10q = qvalue
                   ))
  
  All_HTZ1 <- with(All_HTZ1, 
                    GRanges(chr,
                            IRanges(start,end, names=name),
                            pileup=pileup,
                            FoldEnrichment= fold_enrichment,
                            log10q = qvalue
                    ))
  
  
  
  

 
### make a list called "peak, using two gr range objects defined above 

  peak1 <-list(S39.HTAS1.C, S39.HTZ1)
  peak2 <-list(S39.HTAS1.C, S41.HTAS1.C, S41.HTAS1.N)
  peak3 <-list(S39.HTZ1,S41.HTZ1)
  peak4 <-list(S41.HTAS1.C, S41.HTAS1.N,S41.HTZ1)
  peak5 <- list(All_HTAS1, All_HTZ1)
  
  
### profile ChIP peaks biding to TSS region. UCSC ref Gene; also there is a choice for UCSC.ce11.ensGene
  
  #txdb <- TxDb.Celegans.UCSC.ce11.refGene
  txdb <- TxDb.Celegans.UCSC.ce11.ensGene
  
  
  txdb
  txidinfo <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
  head(txidinfo)
  
  promoter <- getPromoters(TxDb = txdb, upstream=3000, downstream = 3000)
  
  
###  Average Profile of ChIP peaks binding to TSS region

 p1 <- plotAvgProf2(peak1, TxDb=txdb, upstream=3000, downstream=3000,
               xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf = 0.95,  resample = 1000)
  
 p1+scale_color_hue(labels = (c("HTAS-1.C","HTZ-1")))
  
 
  
  p2 <- plotAvgProf2(peak2[], TxDb=txdb, upstream=3000, downstream=3000,
               xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)
  
  p2+scale_color_hue(labels = (c("Rep1.HTAS1.C", "Rep2.HTAS1.C", "Rep2.HTAS1.N")))
  
  
  p3 <- plotAvgProf2(peak3[], TxDb=txdb, upstream=3000, downstream=3000,
               xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)
  
  p3+scale_color_hue(labels = (c("Rep1.HTZ-1", "Rep2.HTZ-1")))
  
  
  
  p4 <- plotAvgProf2(peak4[], TxDb=txdb, upstream=3000, downstream=3000,
               xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)
  
  p4+scale_color_hue(labels = (c("HTAS1.C", "HTAS1.N", "HTZ1")))
  
  
  
  p5 <- plotAvgProf2(peak5, TxDb=txdb, upstream=3000, downstream=3000,
                     xlab="Genomic Region (5'->3')", conf = 0.95, ylab = "Read Count Frequency",  resample = 1000)
  
  p5+scale_color_hue(labels = (c("HTAS-1","HTZ-1")))
  

###   Average Profile of ChIP peaks binding to TSS and TTS region

  ##### peak[[1]] for dup01, and peak[[2]] for auto 
  p1 <- plotPeakProf2(peak5, upstream = 3000, downstream = 3000,
                by = "gene", type = "body", nbin = 100,
                TxDb = txdb,ignore_strand = F)
  
  
  p1+scale_color_hue(labels = (c("HTAS-1","HTZ-1")))
  
  
  p2 <- plotPeakProf2(peak4, upstream = rel(.2), downstream = rel(.3),
                      by = "gene", type = "body", nbin = 100,
                      TxDb = txdb,ignore_strand = F)
  
  
  p2+scale_color_hue(labels = (c("HTAS-1.C", "HTAS-1.N", "HTZ-1")))
  
  
###### peak annotation ############ here I define tssRegions as -1000bp upstream and 300bp downstream

  peakAnno.auto <- annotatePeak(All_HTZ1[], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Ce.eg.db")

  peakAnno.df.auto <-as.data.frame(peakAnno.auto)
  write.csv(peakAnno.df.auto, "peakAnno.HTZ1.S41.csv") 
  
  #do the similar codes for dup01 by using peak[[1]]
  
  plotAnnoPie(peakAnno.auto) #or peakAnno.auto

  plotAnnoBar(peakAnno.auto) #or peakAnno.auto
  
  vennpie(peakAnno.auto, r=0.2, cex=0.8) #or peakAnno.auto
  
  
 ##### UPSET PLOT
 
  library(UpSetR)
  library(ggupset)
  library(ggplotify)
  library(ggimage)
  
  peakAnno.auto <- annotatePeak(S39.HTAS1.C[], tssRegion=c(-1000, 3000), TxDb=txdb, annoDb="org.Ce.eg.db")
  upsetplot(peakAnno.auto)
  write.csv(peakAnno.auto, "test.csv") 
 
  peakAnno.auto <- annotatePeak(S39.HTZ1[], tssRegion=c(-1000, 3000), TxDb=txdb, annoDb="org.Ce.eg.db")
  upsetplot(peakAnno.auto)

  
  peakAnno.auto <- annotatePeak(S41.HTAS1.C[], tssRegion=c(-1000, 3000), TxDb=txdb, annoDb="org.Ce.eg.db")
  upsetplot(peakAnno.auto)

  
  peakAnno.auto <- annotatePeak(S41.HTAS1.N[], tssRegion=c(-1000, 3000), TxDb=txdb, annoDb="org.Ce.eg.db")
  upsetplot(peakAnno.auto)

  
  peakAnno.auto <- annotatePeak(S41.HTZ1[], tssRegion=c(-1000, 3000), TxDb=txdb, annoDb="org.Ce.eg.db")
  upsetplot(peakAnno.auto)

   
  plotDistToTSS(peakAnno.auto)
  
###### LIST ######

  
  
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  tagMatrixList <- lapply(peak1, getTagMatrix, windows=promoter)
  
  plotAvgProf(tagMatrixList, xlim=c(-3000,3000))
  
  plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
  
  tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

  
  
 ###### START HERE: comparison of annotated genes #####
  rm(list=ls()) #removes all objects
  setwd("~/Documents/wd/Chu/ChIPseq.revisit.May2021")
  dir()
  
  #note: modified column to with Excel with Left function 
      # 3' UTR is now "3"
      # 5' UTR is now "5"
      # Distal Intergenic is now "D"
      # Downstream (<1kb) is now "D"
      # Downstream (1-2kb) is now "D"
      # Downstream (2-3kb) is now "D"
      # Exon (bunch) is now "E"
      # Intron (bunch) is now "I"
      # Promoter is now "P"
  
  auto <-read.csv("3.peakAnno.auto.csv")
  dup01<- read.csv("3.peakAnno.dup01.csv")
    
    #select only the peaks associated with genes 
    auto%<>%filter(annotation2%in%c("3","5","E","I","P"))
    dup01%<>%filter(annotation2%in%c("3","5","E","I","P"))
    
    write.csv(auto,"4.auto.7882.gene-associated.peaks.csv")
    write.csv(dup01,"4.dup01.3030.gene-associated.peaks.csv")
    head(auto)
  
  
  ##### remove duplicates based on geneId columns
  #which one does R choose to keep? The first one that appears?
  auto.unique <-auto %>% distinct(geneId, .keep_all = TRUE)
  dup01.unique <-dup01 %>% distinct(geneId, .keep_all = TRUE)
  
  overlap <- auto.unique$geneId%in%dup01.unique$geneId
  sum(overlap)
  
##### venn
  library(grid)
  library(futile.logger)
  library(VennDiagram)
  library(RColorBrewer)

  
  a <- auto.unique$geneId
  d <- dup01.unique$geneId
  
  venn.diagram(
    x=list(a,d), 
    category.names = c("auto","dup01"),
    filename = "annotation.geneId.auto.vs.dup01_5.2.2021.png",
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    cex = .6,
    cat.cex = 0.6)
    
############ START FROM HERE ##########################################
  
  rm(list=ls()) #removes all objects
  setwd("~/Documents/wd/Chu/ChIPseq.revisit.May2021")
  dir()
  
  auto<-read.csv("4.auto.7882.gene-associated.peaks.csv")
  dup <-read.csv("4.dup01.3030.gene-associated.peaks.csv")
  auto$wbgene <- substr(auto$geneId,1,nchar(auto$geneId)-2)
  dup$wbgene <- substr(dup$geneId,1,nchar(dup$geneId)-2)
    head(auto)
    head(dup)

  a <-unique(auto$wbgene) #5495
  d <-unique(dup$wbgene) #2352
  
  #transfer the csv frile from Chu/1.14.2021.HTSeq.data/wChIPseq_ce11_v1
  DE<- read.csv("1.DEseq2.HTAS.LFC.output.1.28.21.csv")

  data <- DE %>% mutate(ChIPOverlap.auto=ifelse(DE$Wbgene%in%a,1,0),
                        ChIPOverlap.dup =ifelse(DE$Wbgene%in%d,1,0))
  data%<>%mutate(ChIP.Overlap=case_when(
                            ChIPOverlap.auto==1&ChIPOverlap.dup==1 ~ "common",
                            ChIPOverlap.auto==0&ChIPOverlap.dup==0 ~ "none",
                            ChIPOverlap.auto==1&ChIPOverlap.dup==0 ~ "auto.only",
                            ChIPOverlap.auto==0&ChIPOverlap.dup==1 ~ "dup01.only"))
  head(data,20)
  sum(data$ChIP.Overlap=="common")
  sum(data$ChIP.Overlap=="auto.only")
  sum(data$ChIP.Overlap=="dup01.only")
  
  #volcano plot
  #data <-res
  x <-data$HTAS.log2FoldChange
  y<-(-log10(data$HTAS.padj))
  sig<-data$HTAS.padj<0.05
  sig %>%as.numeric() %>% sum(na.rm = T) 
  common <-data$ChIP.Overlap=="common"
  auto.only <-data$ChIP.Overlap=="auto.only"
  dup01.only <-data$ChIP.Overlap=="dup01.only"
  
  ###########need to work on how to move the dots that ###########
  #for (i in volcano$X.log10.mes4.q) if(i >=3.5) volcano$X.log10.mes4.q <- "3.5"
  #volcano$X.log10.mes4.q[volcano$X.log10.mes4.q >= 3.5] <- 3.5
  
  temp <-cbind(x, y)
  plot(temp,ylim=c(0,20), xlim=c(-6,6),
       col= densCols(temp,colramp = colorRampPalette(c("gray","black"))),pch=20, ann=FALSE)
  abline (v=0, h=0)
  abline (v=-log2(2))
  abline (v=log2(2), h=-log10(0.01))
  #points(temp[sig,], col=densCols(temp[sig,],colramp = colorRampPalette(c("pink","deeppink4"))), pch=20)
  points(temp[common,], col=densCols(temp[common,],colramp = colorRampPalette(c("red","red4"))), pch=20)
  points(temp[auto.only,], col=densCols(temp[auto.only,],colramp = colorRampPalette(c("blue","blue4"))), pch=20)
  points(temp[dup01.only,], col=densCols(temp[dup01.only,],colramp = colorRampPalette(c("green","green4"))), pch=20)
  
  
  ######## MA plot  #####
  y <-data$HTAS.log2FoldChange
  x<-data$HTAS.baseMean
  x<-log2(x+1)

  
  temp <-cbind(x, y)
  plot(temp, ylim=c(-5,5),
       col= densCols(temp,colramp = colorRampPalette(c("gray","black"))),pch=20, ann=FALSE)
  abline (v=0, h=0)
  points(temp[common,], col=densCols(temp[common,],colramp = colorRampPalette(c("red","red4"))), pch=20)
  points(temp[auto.only,], col=densCols(temp[auto.only,],colramp = colorRampPalette(c("blue","blue4"))), pch=20)
  points(temp[dup01.only,], col=densCols(temp[dup01.only,],colramp = colorRampPalette(c("green","green4"))), pch=20)
  
  
  
  #venn
  library(grid)
  library(futile.logger)
  library(VennDiagram)
  library(RColorBrewer)
  
  common.genes <-data %>% filter(ChIP.Overlap=="common")
  common.genes<-common.genes$Wbgene
  
  all.auto.genes <-data %>% filter(ChIPOverlap.auto==1)
  all.auto.genes<-all.auto.genes$Wbgene
