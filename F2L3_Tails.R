setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
miRBase <- miRBase[which(miRBase$high_conf==T),]
#armmirbase <- unique(miRBase[,c(1,16)])
#armmirbase <- armmirbase[which(armmirbase$STRAND=="5P"|armmirbase$STRAND=="3P"),]

setwd("~/Desktop/P6 Tiejuan/02_F2L3 AGO2IP comparison/")
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
sequences$percent <- gsub("%", "", sequences$RATIO)
sequences$percent <- as.numeric(sequences$percent)
sequences <- unique(sequences)
sequences <- sequences[which(sequences$MATCH==""),]
sequences$SEQ_TAIL <- as.character(sequences$SEQ_TAIL)
sequences$Tscore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("T","",sequences$SEQ_TAIL))
sequences$Ascore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("A","",sequences$SEQ_TAIL))
sequences$Gscore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("G","",sequences$SEQ_TAIL))
sequences$Cscore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("C","",sequences$SEQ_TAIL))

summary <- read.table(file="comparison3_F2L3_WT_AGO2IP.isomir.tsv", header = TRUE, fill = TRUE, sep =  "\t")
summary$CPM <- summary$TOTAL_READS/summary$TOTAL_READS_IN_SAMPLE*1000000
summary$SAMPLE <- gsub("-endo.fastq_ready","",summary$SAMPLE)
summary$SAMPLE <- gsub("WT-","",summary$SAMPLE)

sequences <- sequences[ (sequences$MIRNA %in% miRBase$MIRNA),]
summary <- summary[ (summary$MIRNA %in% miRBase$MIRNA),]

global2 <- summary[,c(1:2)]
i <- 1
global2$CPM_total <- NA
global2$CPM_canonical <- NA
global2$CPM_tailed <- NA
global2$CPM_trimmed <- NA
global2$CPM_bothTT <- NA

for (i in 1:nrow(global2)) {
  miRNA_test <- as.character(global2$MIRNA[i])
  sample_test <- as.character(global2$SAMPLE[i])
  
  sequences_test <- sequences[which(sequences$MIRNA==miRNA_test & sequences$SAMPLE==sample_test),]
  CPM_total <- sum(sequences_test$CPM)
  #CPM_canonical <- max(sequences_test[which(sequences_test$templated>0),12])
  CPM_canonical <- sum(sequences_test[which(sequences_test$LEN_TRIM==0 & sequences_test$LEN_TAIL==0),12])
  CPM_templed_tail <- sum(sequences_test[which(sequences_test$LEN_TRIM==0 & sequences_test$LEN_TAIL>0 & sequences_test$templated>0),12])
  CPM_tailed <- sum(sequences_test[which(sequences_test$LEN_TRIM==0 & sequences_test$LEN_TAIL>0 & sequences_test$templated==0),12])
  CPM_trimmed <- sum(sequences_test[which(sequences_test$LEN_TRIM>0 & sequences_test$LEN_TAIL==0),12])
  CPM_bothTT <- sum(sequences_test[which(sequences_test$LEN_TRIM>0 & sequences_test$LEN_TAIL>0),12])
  
  global2$CPM_total[i] <- CPM_total
  global2$CPM_canonical[i] <- CPM_canonical+CPM_templed_tail
  global2$CPM_tailed[i] <- CPM_tailed
  global2$CPM_trimmed[i] <- CPM_trimmed
  global2$CPM_bothTT[i] <- CPM_bothTT
  print(i)
  rm(sequences_test,miRNA_test,sample_test,CPM_total, CPM_canonical, CPM_tailed, CPM_trimmed, CPM_bothTT,CPM_templed_tail)
}

total_mir <- aggregate(global2$CPM_canonical, by=list(global2$MIRNA), FUN=sum)
colnames(total_mir) <- c("MIRNA","CPM_avg")
total_mir$CPM_avg <- total_mir$CPM_avg/6
total_mir <- total_mir[order(-total_mir$CPM_avg),]
total_mir <- unique(total_mir[c(1:100),])

sequences <- sequences[ (sequences$MIRNA %in% total_mir$MIRNA),]
summary <- summary[ (summary$MIRNA %in% total_mir$MIRNA),]
global2 <- global2[ (global2$MIRNA %in% total_mir$MIRNA),]

plot(log10(global2$CPM_canonical),log10(global2$CPM_tailed))
plot(log10(global2$CPM_canonical),log10(global2$CPM_trimmed))
plot(log10(global2$CPM_canonical),log10(global2$CPM_bothTT))
plot(log10(global2$CPM_canonical),log10(global2$CPM_total))
plot(log10(global2$CPM_total),log10(global2$CPM_canonical))
plot(log10(global2$CPM_tailed),log10(global2$CPM_trimmed))
abline(a=0, b=1)

global2$percent_tailNT <- global2$CPM_tailed/global2$CPM_total*100
global2$percent_trimmed <- global2$CPM_trimmed/global2$CPM_total*100

WT <- global2[which(global2$SAMPLE=="WT"),c(2:9)]
F2L3 <- global2[which(global2$SAMPLE=="F2L3"),c(2:9)]
merged_Ctrl <- merge(WT, F2L3, by="MIRNA")
wilcox.test(merged_Ctrl$percent_trimmed.x, merged_Ctrl$percent_trimmed.y, paired = T, alternative = "less")
wilcox.test(merged_Ctrl$percent_tailNT.x, merged_Ctrl$percent_tailNT.y, paired = T, alternative = "greater")
#write.table(merged_Ctrl, "CPM_by_type_F2L3_mutant.tsv", sep="\t", append = FALSE, row.names = F)

WT_DIS3L2 <- global2[which(global2$SAMPLE=="DISKO-WT"),c(2:9)]
F2L3_DIS3L2 <- global2[which(global2$SAMPLE=="DISKO-F2L3"),c(2:9)]
merged_DIS3L2 <- merge(WT_DIS3L2, F2L3_DIS3L2, by="MIRNA")
#write.table(merged_DIS3L2, "CPM_by_type_DIS3L2_F2L3_mutant.tsv", sep="\t", append = FALSE, row.names = F)

merged_F2L3 <- merge(F2L3, F2L3_DIS3L2, by="MIRNA")
#write.table(merged_F2L3, "CPM_by_type_F2L3_DIS3L2effect_mutant.tsv", sep="\t", append = FALSE, row.names = F)


#miRNA_selected <- "hsa-miR-20a-5p"
#miRNA_selected <- "hsa-miR-92a-3p-1-2"
#miRNA_selected <- "hsa-miR-7-5p-1-2-3"
#miRNA_selected <- "hsa-miR-103a-3p-1-2"
#miRNA_selected <- "hsa-miR-221-3p"
miRNA_selected <- "hsa-miR-148a-3p"
miRNA_selected_group <- sequences[which(sequences$MIRNA==miRNA_selected & sequences$LEN_TAIL>=2),]
miRNA_selected_T_scores <- aggregate(miRNA_selected_group$CPM, by=list(SAMPLE=miRNA_selected_group$SAMPLE,Tscore=miRNA_selected_group$Tscore), FUN=sum)
miRNA_selected_T_scores2 <- aggregate(miRNA_selected_group$CPM, by=list(SAMPLE=miRNA_selected_group$SAMPLE), FUN=sum)
miRNA_selected_T_scores <- merge(miRNA_selected_T_scores, miRNA_selected_T_scores2, by="SAMPLE")
miRNA_selected_T_scores$percent <- miRNA_selected_T_scores$x.x/miRNA_selected_T_scores$x.y*100
a <- as.character(unique(miRNA_selected_T_scores$SAMPLE))
b <- c(0:max(miRNA_selected_T_scores$Tscore))
empty_table <- expand.grid(a,b)
colnames(empty_table) <- c("SAMPLE","Tscore")
miRNA_selected_T_scores <- merge(miRNA_selected_T_scores, empty_table, by=c("SAMPLE","Tscore"), all = T)
rm(miRNA_selected_T_scores2,miRNA_selected_group,empty_table,a,b)
miRNA_selected_T_scores[is.na(miRNA_selected_T_scores)] <- 0
miRNA_selected_T_scores <- miRNA_selected_T_scores[,c(1,2,5)]
miRNA_selected_T_scores <- miRNA_selected_T_scores[order(miRNA_selected_T_scores$SAMPLE, miRNA_selected_T_scores$Tscore),]
setwd("~/Desktop/P6 Tiejuan/02_F2L3 AGO2IP comparison/F2L3_Tscore_individial_miRNA/")
write.table(miRNA_selected_T_scores, paste0("T_score_",miRNA_selected,".tsv"), sep="\t", append = FALSE, row.names = F)


WT <- summary[which(summary$SAMPLE=="WT"),c(2,12)]
F2L3 <- summary[which(summary$SAMPLE=="F2L3"),c(2,12)]
WT_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-WT"),c(2,12)]
F2L3_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-F2L3"),c(2,12)]
TRIM_only <- merge(WT,F2L3, by="MIRNA")
TRIM_only <- merge(TRIM_only,WT_DIS3L2, by="MIRNA")
TRIM_only <- merge(TRIM_only,F2L3_DIS3L2, by="MIRNA")
colnames(TRIM_only) <- c("MIRNA","WT","F2L3","WT_DIS3L2KO","F2L3_DIS3L2KO")
write.table(TRIM_only, "Percentage_TRIM_only_F2L3_DIS3L2KO.tsv", sep="\t", append = FALSE, row.names = F)
wilcox.test(TRIM_only$F2L3, TRIM_only$F2L3_DIS3L2KO, paired = T)
wilcox.test(TRIM_only$WT, TRIM_only$WT_DIS3L2KO, paired = T)

WT <- summary[which(summary$SAMPLE=="WT"),c(2,14)]
F2L3 <- summary[which(summary$SAMPLE=="F2L3"),c(2,14)]
WT_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-WT"),c(2,14)]
F2L3_DIS3L2 <- summary[which(summary$SAMPLE=="DISKO-F2L3"),c(2,14)]
TAIL_only <- merge(WT,F2L3, by="MIRNA")
TAIL_only <- merge(TAIL_only,WT_DIS3L2, by="MIRNA")
TAIL_only <- merge(TAIL_only,F2L3_DIS3L2, by="MIRNA")
TAIL_only$dif_DIS3L2 <- TAIL_only$F2L3_DIS3L2KO-TAIL_only$F2L3
colnames(TAIL_only) <- c("MIRNA","WT","F2L3","WT_DIS3L2KO","F2L3_DIS3L2KO")
wilcox.test(TAIL_only$F2L3, TAIL_only$F2L3_DIS3L2KO, paired = T)
wilcox.test(TAIL_only$WT, TAIL_only$WT_DIS3L2KO, paired = T, alternative = "greater")
write.table(TAIL_only, "Percentage_TAIL_only_F2L3_DIS3L2KO.tsv", sep="\t", append = FALSE, row.names = F)
