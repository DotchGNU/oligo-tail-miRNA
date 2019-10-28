setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
armmirbase <- unique(miRBase[,c(1,16)])
armmirbase <- armmirbase[which(armmirbase$STRAND=="5P"|armmirbase$STRAND=="3P"),]

setwd("~/Desktop/R_scripts/F2L3 AGO2IP comparison/")
# sequences <- read.table(file="comparison3_F2L3_WT_AGO2IP.isomir.sequence_info.tsv", header = TRUE, fill = TRUE, sep =  "\t")
# sequences$percent <- gsub("%", "", sequences$RATIO)
# sequences$percent <- as.numeric(sequences$percent)
# 
# sequences$SAMPLE <- gsub("WT-WT","WT",sequences$SAMPLE)
# sequences$SAMPLE <- gsub("WT-F2L3","F2L3",sequences$SAMPLE)
# sequences$SAMPLE <- gsub("WT-PAZ","PAZ",sequences$SAMPLE)
# 
# sequences$SAMPLE <- gsub("-endo.fastq_ready","",sequences$SAMPLE)
# 
# sequences$templated <- NA
# i <- 1
# #i <- 13028
# for (i in 1:nrow(sequences)) {
#   sequence_tested <- as.character(sequences$SEQUENCE[i])
#   found <- length(grep(sequence_tested, miRBase$EXTENDED.SEQUENCE))
#   if (found==0) {
#     sequences$templated[i] <- FALSE
#   }
#   if (found!=0) {
#     sequences$templated[i] <- TRUE
#   }
#   print(i/nrow(sequences)*100)
# }
# rm(found, i, sequence_tested)
# write.table(sequences, "comparison3_F2L3_WT_AGO2IP.isomir.sequence_info_templ_untempl.tsv", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)
sequences <- read.table(file="comparison3_F2L3_WT_AGO2IP.isomir.sequence_info_templ_untempl.tsv", header = TRUE, fill = TRUE, sep =  "\t")

summary <- read.table(file="comparison3_F2L3_WT_AGO2IP.isomir.tsv", header = TRUE, fill = TRUE, sep =  "\t")
summary$CPM <- summary$TOTAL_READS/summary$TOTAL_READS_IN_SAMPLE*1000000

summary$SAMPLE <- gsub("WT-WT","WT",summary$SAMPLE)
summary$SAMPLE <- gsub("WT-F2L3","F2L3",summary$SAMPLE)
summary$SAMPLE <- gsub("WT-PAZ","PAZ",summary$SAMPLE)
summary$SAMPLE <- gsub("-endo.fastq_ready","",summary$SAMPLE)

summary_WT <- summary[which(summary$SAMPLE=="WT"),c(2,18)]
summary_WT <- summary_WT[order(-summary_WT$CPM),] 

sizes_sequences <- aggregate(sequences$percent, by=list(sequences$SAMPLE, sequences$MIRNA, sequences$LEN_READ), FUN=sum)
sizes_WT <- sizes_sequences[sizes_sequences$Group.1=="WT",c(2:4)]
colnames(sizes_WT) <- c("MIRNA", "nt","percent_WT")
sizes_F2L3 <- sizes_sequences[sizes_sequences$Group.1=="F2L3",c(2:4)]
colnames(sizes_F2L3) <- c("MIRNA", "nt","percent_F2L3")
sizes_WT_F2L3 <- merge(sizes_WT, sizes_F2L3, by=c("MIRNA","nt"))
rm(sizes_sequences, sizes_WT, sizes_F2L3)







summary_WT <- summary[which(summary$SAMPLE=="WT"),c(2,18)]
summary_F2L3 <- summary[which(summary$SAMPLE=="F2L3"),c(2,18)]
summary_PAZ <- summary[which(summary$SAMPLE=="PAZ"),c(2,18)]
summary_WT_F2L3 <- merge(summary_WT, summary_F2L3, by="MIRNA")
summary_WT_F2L3 <- merge(summary_WT_F2L3, summary_PAZ, by="MIRNA")
colnames(summary_WT_F2L3) <- c("MIRNA", "CPM_WT","CPM_F2L3","CPM_PAZ")
summary_WT_F2L3 <-  merge(summary_WT_F2L3, armmirbase, by="MIRNA")
summary_WT_F2L3$FC_F2L3_293T <- log2(summary_WT_F2L3$CPM_F2L3/summary_WT_F2L3$CPM_WT)
summary_WT_F2L3$FC_PAZ_293T <- log2(summary_WT_F2L3$CPM_PAZ/summary_WT_F2L3$CPM_WT)

#summary_WT_F2L3 <- summary_WT_F2L3[which(summary_WT_F2L3$CPM_WT>=median(summary_WT_F2L3$CPM_WT)),]
summary_strand <- aggregate(summary_WT_F2L3$FC, by=list(summary_WT_F2L3$STRAND), FUN=mean)
p15 <- ggplot(summary_WT_F2L3, aes(x = log10(CPM_WT), y = log10(CPM_F2L3), fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
p15
write.table(summary_WT_F2L3, "CPM_F2L3_WT_AGO2IP.tsv", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)
test3 <- ggplot(summary_WT_F2L3, aes(x = log10(CPM_WT), y = FC_F2L3_293T, fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
test3


p152 <- ggplot(summary_WT_F2L3, aes(x = log10(CPM_WT), y = log10(CPM_PAZ), fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
p152
rm(summary_WT, summary_F2L3, summary_PAZ)

summary_WT_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-WT"),c(2,18)]
summary_F2L3_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-F2L3"),c(2,18)]
summary_PAZ_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-PAZ"),c(2,18)]
summary_WT_F2L3_DIS3L2 <- merge(summary_WT_DIS3L2, summary_F2L3_DIS3L2, by="MIRNA")
summary_WT_F2L3_DIS3L2 <- merge(summary_WT_F2L3_DIS3L2, summary_PAZ_DIS3L2, by="MIRNA")
colnames(summary_WT_F2L3_DIS3L2) <- c("MIRNA", "CPM_WT","CPM_F2L3","CPM_PAZ")
summary_WT_F2L3_DIS3L2 <-  merge(summary_WT_F2L3_DIS3L2, armmirbase, by="MIRNA")
summary_WT_F2L3_DIS3L2$FC_F2L3_DIS3L2 <- log2(summary_WT_F2L3_DIS3L2$CPM_F2L3/summary_WT_F2L3_DIS3L2$CPM_WT)
summary_WT_F2L3_DIS3L2$FC_PAZ_DIS3L2 <- log2(summary_WT_F2L3_DIS3L2$CPM_PAZ/summary_WT_F2L3_DIS3L2$CPM_WT)
p15bis <- ggplot(summary_WT_F2L3_DIS3L2, aes(x = log10(CPM_WT), y = log10(CPM_F2L3), fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
p15bis
p15bis2 <- ggplot(summary_WT_F2L3_DIS3L2, aes(x = log10(CPM_WT), y = log10(CPM_PAZ), fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
p15bis2
rm(summary_WT_DIS3L2, summary_F2L3_DIS3L2, summary_PAZ_DIS3L2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(p15,p152, p15bis,p15bis2, cols=2)

summary_WT_F2L3_DIS3L2 <- summary_WT_F2L3_DIS3L2[,c(1,6,7)]
summary_WT_F2L3 <- merge(summary_WT_F2L3, summary_WT_F2L3_DIS3L2, by="MIRNA")
summary_WT_F2L3 <- unique(summary_WT_F2L3)

summary_WT_F2L3 <- summary_WT_F2L3[order(-summary_WT_F2L3$CPM_WT),] 
write.table(summary_WT_F2L3, "FC_F2L3_PAZ_293T_DIS3L2KO.tsv", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)



selected <- summary_WT_F2L3[which(summary_WT_F2L3$CPM_WT>=median(summary_WT_F2L3$CPM_WT)),]
#selected <- summary_WT_F2L3[which(summary_WT_F2L3$CPM_WT>=150),]

breaks <-  seq(-10, 10, by=0.001)

datLog2FC <-  as.numeric(selected$FC_F2L3_293T) #FC
cumfreq0 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(selected))) 
plot(breaks, cumfreq0, type="n", xlab = "fold-change Log2", ylab = "Cumulative fraction", xlim = c(-5,2) )
lines(breaks, cumfreq0)

datLog2FC <-  as.numeric(selected$FC_PAZ_293T) #FC
cumfreq1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(selected)))
lines(breaks, cumfreq1, col="red")

datLog2FC <-  as.numeric(selected$FC_F2L3_DIS3L2) #FC
cumfreq2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(selected))) 
lines(breaks, cumfreq2, col="blue")

datLog2FC <-  as.numeric(selected$FC_PAZ_DIS3L2) #FC
cumfreq3 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(selected))) 
lines(breaks, cumfreq3, col="green")

wilcox.test(selected$FC_F2L3_293T, selected$FC_F2L3_DIS3L2)
wilcox.test(selected$FC_F2L3_293T, selected$FC_PAZ_293T)
t.test(selected$FC_F2L3_293T,selected$FC_F2L3_DIS3L2,paired=TRUE)
t.test(selected$FC_PAZ_293T,selected$FC_PAZ_DIS3L2,paired=TRUE)
t.test(selected$FC_PAZ_293T,mu=0)

