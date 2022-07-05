# ChIPseeker:To annotate and analyze ChIP peaks. 

## installed ChIPseeker, verion 1.28.3

      library(ChIPseeker)
      library(TxDb.Celegans.UCSC.ce11.refGene)
      library(TxDb.Celegans.UCSC.ce11.ensGene)
      library(clusterProfiler)
      library(GenomicRanges)
      library(magrittr)
      library(tidyverse)
      library(ggplot)

        rm(list=ls()) #removes all objects
        setwd("/Users/swadh/Desktop/chipseeker/to_graduate/distribution/")
        dir()
        
        
       #2 read the output Excel files (csv) from MACS2
       #note: make sure to modify the column names in Excel (csv) files
        S39.HTAS1.downstream <-read.csv("S39.HTAS-1.downstream.csv")
        S39.HTAS1.overlap <-read.csv("S39.HTAS-1.overlap.csv")
        S39.HTAS1.upstream <-read.csv("S39.HTAS-1.upstream.csv")
        S39.HTZ1.downstream <-read.csv("S39.HTZ-1.downstream.csv")
        S39.HTZ1.overlap <-read.csv("S39.HTZ-1.overlap.csv")
        S39.HTZ1.upstream <-read.csv("S39.HTZ-1.upstream.csv")
        HTAS1_S39 <-read.csv("HTAS-1_S39.csv")
        HTZ1_S39.csv <-read.csv("HTZ-1_S39.csv") 


        head(whittle_HTZ1_emb.peaks1)
        head(ChIP.HTZ1_S41)
        
        #analysis of each column. Let's check Peak Length
        summary(whittle_HTZ1_emb.peaks1)
        hist(whittle_HTZ1_emb.peaks1$length, border="black", col="lightblue",
             xlim=c(0,3000), breaks=100)

        #analysis of each column. Let's check Lo10.qvalue
        summary(whittle_HTZ1_emb.peaks1$LOG10.qvalue)
        hist(whittle_HTZ1_emb.peaks1$qvalue, border="black", col="lightblue",
             xlim=c(0,20), breaks=100)

        #analysis of each column. Let's check "fold_enrichment
        summary(ChIP.HTAS1.C_S39$fold_enrichment)
        hist(ChIP.HTAS1.C_S39$fold_enrichment, border="black", col="lightblue",
             xlim=c(0,6), breaks=10)



          # making gr range object for auto
            S39.HTAS1.downstream <- with(S39.HTAS1.downstream, 
                              GRanges(chr,
                                      IRanges(start,end),
                              length=length,
                              log10q = log10q
                              ))
            S39.HTAS1.downstream


            S39.HTAS1.overlap <- with(S39.HTAS1.overlap, 
                                     GRanges(chr,
                                             IRanges(start,end, names=name),
                                             pileup=pileup,
                                             FoldEnrichment= fold_enrichment,
                                             log10q = log10q
                                     ))
            S39.HTAS1.overlap


            S39.HTAS1.upstream <- with(S39.HTAS1.upstream, 
                                     GRanges(chr,
                                             IRanges(start,end, names=name),
                                             pileup=pileup,
                                             FoldEnrichment= fold_enrichment,
                                             log10q = log10q
                                     ))
            S39.HTAS1.upstream

            S39.HTZ1.downstream <- with(S39.HTZ1.downstream, 
                                     GRanges(chr,
                                             IRanges(start,end, names=name),
                                             pileup=pileup,
                                             FoldEnrichment= fold_enrichment,
                                             log10q = log10q
                                     ))
            S39.HTZ1.downstream

            S39.HTZ1.overlap <- with(S39.HTZ1.overlap, 
                                  GRanges(chr,
                                          IRanges(start,end, names=name),
                                          pileup=pileup,
                                          FoldEnrichment= fold_enrichment,
                                          log10q = log10q
                                  ))
            S39.HTZ1.overlap

            S39.HTZ1.upstream <- with(S39.HTZ1.upstream, 
                                     GRanges(chr,
                                             IRanges(start,end, names=name),
                                             pileup=pileup,
                                             FoldEnrichment= fold_enrichment,
                                             log10q = log10q
                                     ))
            S39.HTZ1.upstream

            HTAS1_S39 <- with(HTAS1_S39, 
                                         GRanges(chr,
                                                 IRanges(start,end),
                                                 length=length,
                                                 log10q = log10q
                                         ))
            HTAS1_S39

            HTZ1_S39.csv <- with(HTZ1_S39.csv, 
                              GRanges(chr,
                                      IRanges(start,end),
                                      length=length,
                                      log10q = log10q
                              ))
            HTZ1_S39.csv

  
 
            #make a list called "peak, using two gr range objects defined above 
            peak1 <-list(S39.HTZ1.upstream,S39.HTZ1.overlap,S39.HTZ1.downstream)
            peak2 <-list(S39.HTAS1.upstream,S39.HTAS1.overlap,S39.HTAS1.downstream)
            peak3 <-list(S39.HTZ1.upstream,S39.HTZ1.overlap,S39.HTZ1.downstream,S39.HTAS1.upstream,S39.HTAS1.overlap,S39.HTAS1.downstream)
            peak4 <-list(HTZ1_S39.csv,HTAS1_S39)

            peak4[]

            #plot cov plot 
            # peak is a list with 1=dup01, and 2=auto.
            # here, I am calling 2 = auto in the list called "peak" to plot it in covPlot
            covplot(peak[[1]], weightCol = "qvalue")


          # profile ChIP peaks biding to TSS regions
            # UCSC ref Gene; also there is a choice for UCSC.ce11.ensGene

            #txdb <- TxDb.Celegans.UCSC.ce11.refGene
            txdb <- TxDb.Celegans.UCSC.ce11.ensGene
            txdb
            txidinfo <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
            head(txidinfo)

            promoter <- getPromoters(TxDb = txdb, upstream=5000, downstream = 5000)
            #tagMatrix1 <- getTagMatrix(peak[[1]], windows=promoter) # this does not work!
            #tagMatrix2 <- getTagMatrix(peak[[2]], windows=promoter) # this does not work!

            #tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red") # FAIL
            peakHeatmap(peak1[[1]], TxDb=txdb, upstream=3000, downstream=3000, color="red")


          #Average Profile of ChIP peaks binding to TSS region
          # peak[[1]] for dup01, and peak[[2]] for auto  
            plotAvgProf2(peak1[], TxDb=txdb, upstream=3000, downstream=3000,
                         xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)


            plotAvgProf2(peak2[], TxDb=txdb, upstream=3000, downstream=3000,
                         xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)

            plotAvgProf2(peak3[], TxDb=txdb, upstream=3000, downstream=3000,
                         xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)

            plotAvgProf2(peak4[], TxDb=txdb, upstream=10000, downstream=10000,
                         xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", resample = 1000)

          #### peak annotation ############
          # default: -3kb and +3kb TSS regions 
          # here I define tssRegions as -1000bp upstream and 300bp downstream
            files <- HTZ1_S39.csv
            print(HTZ1_S39.csv)
            peakAnno.auto <- annotatePeak(peak4[], tssRegion=c(-1000, 5000), TxDb=txdb, annoDb="org.Ce.eg.db")

            peakAnno.df.auto <-as.data.frame(peakAnno.auto)
            write.csv(peakAnno.df.auto, "peakAnno.HTAS1.csv") 

            #do the similar codes for dup01 by using peak[[1]]

            plotAnnoPie(peakAnno.auto) #or peakAnno.auto

            plotAnnoBar(peakAnno.auto) #or peakAnno.auto

            vennpie(peakAnno.auto, r=0.2, cex=0.8) #or peakAnno.auto

            library(UpSetR)
            library(ggupset)
            library(ggplotify)
            library(ggimage)

            upsetplot(peakAnno)
            upsetplot(peakAnno.dup01, vennpie=TRUE)

            plotDistToTSS(peakAnno)

          #### LIST ######



            promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
            tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter)

            plotAvgProf(tagMatrixList, xlim=c(-3000,3000))

            plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

            tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

  
  
          #### START HERE: comparison of annotated genes #####
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


          #remove duplicates based on geneId columns
            #which one does R choose to keep? The first one that appears?
            auto.unique <-auto %>% distinct(geneId, .keep_all = TRUE)
            dup01.unique <-dup01 %>% distinct(geneId, .keep_all = TRUE)

            overlap <- auto.unique$geneId%in%dup01.unique$geneId
            sum(overlap)

            #venn
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
    
                    ########## START FROM HERE ##########################################

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

                      # need to work on how to move the dots that 
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


                      # MA plot 
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

                      up.genes <-data %>% filter(HTAS.padj<0.01&HTAS.log2FoldChange>0)
                      up.genes<-up.genes$Wbgene #191

                      down.genes <-data %>% filter(HTAS.padj<0.01&HTAS.log2FoldChange<0)
                      down.genes<-down.genes$Wbgene #343

                      myList <-list(common.genes,up.genes,down.genes)
                      Reduce(intersect,list(common.genes,up.genes,down.genes))
                      intersect(myList[[1]],myList[[3]])

                      venn.diagram(
                        x=list(common.genes,up.genes,down.genes), 
                        category.names = c("bound.common.genes (1640)","up.genes (191)","down.genes(343)"),
                        filename = "S39.HTAS1C_6.11.2021.png",
                        # Output features
                        imagetype="png" ,
                        height = 480 , 
                        width = 480 , 
                        resolution = 300,
                        compression = "lzw",
                        cex = .6,
                        cat.cex = 0.6)

                      venn.diagram(
                        x=list(all.auto.genes,up.genes,down.genes), 
                        category.names = c("all.auto.genes (5495)","up.genes (191)","down.genes(343)"),
                        filename = "test.png",
                        # Output features
                        imagetype="png" ,
                        height = 480 , 
                        width = 480 , 
                        resolution = 300,
                        compression = "lzw",
                        cex = .6,
                        cat.cex = 0.6)





                      venn.diagram(
                        x=list(common.genes,up.genes), 
                        category.names = c("bound.common.genes (1640)","up.genes (191)"),
                        filename = "test.png",
                        # Output features
                        imagetype="png" ,
                        height = 480 , 
                        width = 480 , 
                        resolution = 300,
                        compression = "lzw",
                        cex = .6,
                        cat.cex = 0.6)


                      venn.diagram(
                        x=list(common.genes,down.genes), 
                        category.names = c("bound.common.genes (1640)","down.genes(343)"),
                        filename = "test2.png",
                        # Output features
                        imagetype="png" ,
                        height = 480 , 
                        width = 480 , 
                        resolution = 300,
                        compression = "lzw",
                        cex = .6,
                        cat.cex = 0.6)



                          ######################################################################  
                          # ways to install special packages from bioconductor

                          #install TxDb.Celegans.UCSC.ce11.refGene
                            if (!requireNamespace("BiocManager", quietly = TRUE))
                              install.packages("BiocManager")
                            BiocManager::install("TxDb.Celegans.UCSC.ce11.refGene")

                          #install TxDb.Celegans.UCSC.ce11.ensGene
                            if (!requireNamespace("BiocManager", quietly = TRUE))
                              install.packages("BiocManager")
                            BiocManager::install("TxDb.Celegans.UCSC.ce11.ensGene")

                            if (!requireNamespace("BiocManager", quietly = TRUE))
                              install.packages("BiocManager")

                          #install clusterProfiler
                            if (!requireNamespace("BiocManager", quietly = TRUE))
                              install.packages("BiocManager")
                            BiocManager::install("clusterProfiler")

                          #install genome wide annotation for worm
                          if (!requireNamespace("BiocManager", quietly = TRUE))
                              install.packages("BiocManager")
                            BiocManager::install("org.Ce.eg.db")

                          # reactomePA
                          if (!requireNamespace("BiocManager", quietly = TRUE))
                              install.packages("BiocManager")
                            BiocManager::install("ReactomePA")

