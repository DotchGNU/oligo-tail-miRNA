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
summary$SAMPLE <- gsub("-endo.fastq_ready","",summary$SAMPLE)

summary_WT <- summary[which(summary$SAMPLE=="WT"),c(2,18)]
summary_F2L3 <- summary[which(summary$SAMPLE=="F2L3"),c(2,18)]
summary_WT_F2L3 <- merge(summary_WT, summary_F2L3, by="MIRNA")
colnames(summary_WT_F2L3) <- c("MIRNA", "CPM_WT","CPM_F2L3")
plot(log10(summary_WT_F2L3$CPM_WT),log10(summary_WT_F2L3$CPM_F2L3))
summary_WT_F2L3 <-  merge(summary_WT_F2L3, armmirbase, by="MIRNA")
summary_WT_F2L3$FC <- summary_WT_F2L3$CPM_F2L3/summary_WT_F2L3$CPM_WT
#summary_WT_F2L3 <- summary_WT_F2L3[which(summary_WT_F2L3$CPM_WT>=median(summary_WT_F2L3$CPM_WT)),]
summary_strand <- aggregate(summary_WT_F2L3$FC, by=list(summary_WT_F2L3$STRAND), FUN=mean)
p15 <- ggplot(summary_WT_F2L3, aes(x = log10(CPM_WT), y = log10(CPM_F2L3), fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
p15

summary_WT_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-WT"),c(2,18)]
summary_F2L3_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-F2L3"),c(2,18)]
summary_WT_F2L3_DIS3L2 <- merge(summary_WT_DIS3L2, summary_F2L3_DIS3L2, by="MIRNA")
colnames(summary_WT_F2L3_DIS3L2) <- c("MIRNA", "CPM_WT","CPM_F2L3")
summary_WT_F2L3_DIS3L2 <-  merge(summary_WT_F2L3_DIS3L2, armmirbase, by="MIRNA")
summary_WT_F2L3_DIS3L2$FC_DIS3L2 <- summary_WT_F2L3_DIS3L2$CPM_F2L3/summary_WT_F2L3_DIS3L2$CPM_WT
p15bis <- ggplot(summary_WT_F2L3_DIS3L2, aes(x = log10(CPM_WT), y = log10(CPM_F2L3), fill = STRAND)) + geom_point(aes(color=STRAND)) + geom_abline(intercept = 0)
p15bis


total_mir <- aggregate(summary$CPM, by=list(summary$MIRNA), FUN=sum)
colnames(total_mir) <- c("MIRNA","CPM_avg")
total_mir$CPM_avg <- total_mir$CPM_avg/2
total_mir <- total_mir[order(-total_mir$CPM_avg),]
total_mir_selected <- total_mir[c(1:200),]

sequences$SEQUENCE <- as.character(sequences$SEQUENCE)
sequences$hasN <- nchar(sequences$SEQUENCE)-nchar(gsub("N","",sequences$SEQUENCE))
sequences <- sequences[which(sequences$hasN==0),]

#Canonicals definition 1 (based on miRBase)
# canonicals <- sequences[which(sequences$LEN_TRIM==0 & sequences$LEN_TAIL==0 & sequences$VAR_5P==0),]
# canonicals_WT <- sequences[which(sequences$LEN_TRIM==0 & sequences$LEN_TAIL==0 & sequences$VAR_5P==0 & sequences$SAMPLE=="WT"),c(2,3,12,15)]
# canonicals_F2L3 <- sequences[which(sequences$LEN_TRIM==0 & sequences$LEN_TAIL==0 & sequences$VAR_5P==0 & sequences$SAMPLE=="F2L3"),c(3,12,15)]

#Canonicals definition 2 (based on most abundant templated read)
canonicals <- sequences[which(sequences$templated==TRUE & sequences$SAMPLE=="WT" & sequences$READS>=10),]
#canonicals <- canonicals[ (canonicals$MIRNA %in% total_mir_selected$MIRNA), ]
canonicals_mir <- aggregate(canonicals$CPM, by=list(canonicals$MIRNA), FUN=max)
colnames(canonicals_mir) <- c("MIRNA", "CPM")
canonicals_mir <- merge(canonicals, canonicals_mir, by=c("MIRNA","CPM"))
canonicals_mir <- unique(canonicals_mir)
canonicals_mir <- canonicals_mir[which(canonicals_mir$SEQUENCE!="CTTGGCACCTAGCAAGCACTCA"),]

canonicals <- sequences[ (sequences$SEQUENCE %in% canonicals_mir$SEQUENCE), ]
canonicals <- merge(canonicals, armmirbase, by="MIRNA")
canonicals$percentIsomiR <-  100-canonicals$percent
canonicals <- canonicals[ (canonicals$MIRNA %in% total_mir_selected$MIRNA), ]
canonicals_WT <- canonicals[which(canonicals$SAMPLE=="WT"),c(1,3,12,15)]
canonicals_F2L3 <- canonicals[which(canonicals$SAMPLE=="F2L3"),c(3,12,15)]

total_mir_selected <- merge(total_mir_selected, canonicals_WT, by="MIRNA")
total_mir_selected <- merge(total_mir_selected, canonicals_F2L3, by="SEQUENCE")
total_mir_selected <- unique(total_mir_selected)
total_mir_selected$decrease <- total_mir_selected$percent.x > total_mir_selected$percent.y
t.test(total_mir_selected$percent.x,total_mir_selected$percent.y,paired=TRUE)
summary(total_mir_selected)

total_mir_selected <- merge(total_mir_selected, armmirbase, by="MIRNA")
total_mir_selected$percentIsomirWT <- 100-total_mir_selected$percent.x
total_mir_selected$percentIsomirF2L3 <- 100-total_mir_selected$percent.y
total_mir_selected <- total_mir_selected[order(total_mir_selected$percentIsomirWT),]
write.table(total_mir_selected, "canonicals_isolate_miRNA2.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)

amr5p <- total_mir_selected[which(total_mir_selected$STRAND=="5P"),]
amr3p <- total_mir_selected[which(total_mir_selected$STRAND=="3P"),]
t.test(amr5p$percentIsomirWT,amr5p$percentIsomirF2L3,paired=TRUE)
t.test(amr3p$percentIsomirWT,amr3p$percentIsomirF2L3,paired=TRUE)
write.table(amr5p, "amr5p_isomir_percent.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)
write.table(amr3p, "amr3p_isomir_percent.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)


plot(total_mir_selected$percentIsomirWT, total_mir_selected$percentIsomirF2L3, xlim = c(0,100), ylim = c(0,100))
lines(x = c(0,100), y = c(0,100))

library(ggplot2)
#5p vs 3p comparison
p10 <- ggplot(canonicals, aes(x = SAMPLE, y = percentIsomiR, fill = STRAND)) +
  geom_boxplot() + theme_bw()
p10

# detect sequence preference on generation of isomiRs
trimmed <- sequences[which(sequences$LEN_TRIM>0),]
trimmed_summary <- aggregate(trimmed$percent, by=list(trimmed$MIRNA, trimmed$SAMPLE), FUN=sum)
colnames(trimmed_summary) <- c("MIRNA", "SAMPLE","Trimmed_percentage")
trimmed_summary <- merge(trimmed_summary, armmirbase, by="MIRNA")
trimmed_summary <- trimmed_summary[ (trimmed_summary$MIRNA %in% total_mir_selected$MIRNA), ]
p11 <- ggplot(trimmed_summary, aes(x = SAMPLE, y = Trimmed_percentage, fill = STRAND)) +
  geom_boxplot() + theme_bw()
p11

tailed <- sequences[which(sequences$LEN_TAIL>0),]
tailed_summary <- aggregate(tailed$percent, by=list(tailed$MIRNA, tailed$SAMPLE), FUN=sum)
colnames(tailed_summary) <- c("MIRNA", "SAMPLE","Trailed_percentage")
tailed_summary <- merge(tailed_summary, armmirbase, by="MIRNA")
tailed_summary <- tailed_summary[ (tailed_summary$MIRNA %in% total_mir_selected$MIRNA), ]
p12 <- ggplot(tailed_summary, aes(x = SAMPLE, y = Trailed_percentage, fill = STRAND)) +
  geom_boxplot() + theme_bw()
p12

len_dist <-  aggregate(sequences$CPM, by=list(sequences$LEN_READ, sequences$SAMPLE), FUN=sum)
len_WT <- len_dist[which(len_dist$Group.2=="WT"),]
len_WT$x <- len_WT$x/sum(len_WT$x)*100
len_F2L3 <- len_dist[which(len_dist$Group.2=="F2L3"),]
len_F2L3$x <- len_F2L3$x/sum(len_F2L3$x)*100
len_dist <- merge(len_WT,len_F2L3,by=c("Group.1"))
colnames(len_dist) <- c("LEN_READ","SAMPLE_WT","REL_WT","SAMPLE_F2L3","REL_F2L3")
write.table(len_dist, "len_dist_F2L3.tsv", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)
len_dist2 <- rbind(len_WT, len_F2L3)
colnames(len_dist2) <- c("LEN_READ","SAMPLE","REL")
p13 <- ggplot(len_dist2, aes(x = LEN_READ, y = REL, fill = SAMPLE)) + geom_line(aes(color=SAMPLE))+
                geom_point(aes(color=SAMPLE))+xlim(10,40)
p13


tailed_trimmed <- sequences[which(sequences$LEN_TAIL>0 & sequences$LEN_TRIM>0),]
tailed_trimmed_summary <- aggregate(tailed_trimmed$percent, by=list(tailed_trimmed$MIRNA, tailed_trimmed$SAMPLE), FUN=sum)
colnames(tailed_trimmed_summary) <- c("MIRNA", "SAMPLE","Tailed_trimmed_percentage")
tailed_trimmed_summary <- merge(tailed_trimmed_summary, armmirbase, by="MIRNA")
tailed_trimmed_summary <- tailed_trimmed_summary[ (tailed_trimmed_summary$MIRNA %in% total_mir_selected$MIRNA), ]
p14 <- ggplot(tailed_trimmed_summary, aes(x = SAMPLE, y = Tailed_trimmed_percentage, fill = STRAND)) +
  geom_boxplot() + theme_bw()
p14

WT <- sequences[which(sequences$SAMPLE=="WT"),c(2:4,7:11,14,16,12,15)]
F2L3 <- sequences[which(sequences$SAMPLE=="F2L3"),c(3,12,15)]
#ALL_percent <- merge(WT, F2L3, by="SEQUENCE", all = TRUE)
ALL_percent <- merge(WT, F2L3, by="SEQUENCE", all = FALSE)
ALL_percent$delta <- ALL_percent$percent.y-ALL_percent$percent.x

ALL_percent <- ALL_percent[which(ALL_percent$CPM.x>=100),]
ALL_percent$lastNT <- substr(ALL_percent$SEQUENCE, nchar(ALL_percent$SEQUENCE), nchar(ALL_percent$SEQUENCE))
ALL_percent$firstNT <- substr(ALL_percent$SEQUENCE, 1, 1)

len_delta <- aggregate(ALL_percent$delta, by=list(ALL_percent$LEN_READ), FUN=mean)
NT_delta <- aggregate(ALL_percent$delta, by=list(ALL_percent$lastNT), FUN=mean)
NT1_delta <- aggregate(ALL_percent$delta, by=list(ALL_percent$firstNT), FUN=mean)

plot(log10(ALL_percent$CPM.x), log10(ALL_percent$CPM.y))
lines(x = c(0,100), y = c(0,100))



trimmed_F2L3 <- trimmed_summary[which(trimmed_summary$SAMPLE=="F2L3"),c(1,3)]
trimmed_WT <- trimmed_summary[which(trimmed_summary$SAMPLE=="WT"),c(1,3)]
trimmed_WT_F2L3 <- merge(trimmed_WT, trimmed_F2L3, by="MIRNA")
summary_WT_F2L3 <- merge(summary_WT_F2L3, trimmed_WT_F2L3, by="MIRNA")
summary_WT_F2L3 <- merge(summary_WT_F2L3, armmirbase, by="MIRNA")
summary_WT_F2L3$delta <- summary_WT_F2L3$Trimmed_percentage.y-summary_WT_F2L3$Trimmed_percentage.x
summary_WT_F2L3$delta2 <- summary_WT_F2L3$CPM_F2L3-summary_WT_F2L3$CPM_WT
p16 <- ggplot(summary_WT_F2L3, aes(x = log10(CPM_WT), y = log10(CPM_F2L3), fill = delta)) + geom_point(aes(color=delta)) + geom_abline(intercept = 0)
p16

p17 <- ggplot(summary_WT_F2L3, aes(x = delta, y = delta2, fill = STRAND)) + geom_point(aes(color=STRAND))
p17
