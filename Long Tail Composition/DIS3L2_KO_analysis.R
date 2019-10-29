setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
miRBase <- miRBase[which(miRBase$high_conf==T),]

#### COMPOSITION miRNA
setwd("~/Desktop/P6 Tiejuan/00_DIS3L2_KO_AG2017/new analysis DIS3L2_KO/")

global <- read.table(file="DIS3L2_KO_AGO2IP_1_99_Jun20.isomir.tsv", header = TRUE, fill = TRUE, sep =  "\t")
global <- global[ (global$MIRNA %in% miRBase$MIRNA),]
global$CPM <- global$TOTAL_READS/global$TOTAL_READS_IN_SAMPLE*1000000
global$SET <- substr(global$SAMPLE, 18, 19)
#global <- global[which(global$SET=="WT" | global$SET=="DI"),]
global$SAMPLE <- gsub("Endo-AGO-IP-293T-","",global$SAMPLE)
global$SAMPLE <- gsub("_L001_R1_001.fastq_ready","",global$SAMPLE)
global$SAMPLE <- gsub("_S","",global$SAMPLE)

# sequences <- read.table(file="DIS3L2_KO_AGO2IP_1_99_Jun20.isomir.sequence_info.tsv", header = TRUE, fill = TRUE, sep =  "\t")
# sequences$templated <- NA
# i <- 1
# for (i in 1:nrow(sequences)) {
#   seq_test <- as.character(sequences$SEQUENCE[i])
#   miRBase_test <- miRBase
#   miRBase_test <- grep(seq_test, miRBase_test$PRI.SEQUENCE)
#   miRBase_test <- length(miRBase_test)
# 
#   sequences$templated[i] <- miRBase_test
#   print(i)
#   rm(seq_test, miRBase_test)
# }
# rm(found, i, sequence_tested)
# write.table(sequences, "DIS3L2_KO_AGO2IP_1_99_Jun20.isomir.sequence_info_v2.tsv", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)
sequences <- read.table(file="DIS3L2_KO_AGO2IP_1_99_Jun20.isomir.sequence_info_v2.tsv", header = TRUE, fill = TRUE, sep =  "\t")
sequences <- sequences[ (sequences$MIRNA %in% miRBase$MIRNA),]

sequences <- sequences[which(sequences$MATCH==""),]
sequences$SET <- substr(sequences$SAMPLE, 18, 19)
#sequences <- sequences[which(sequences$SET=="WT" | sequences$SET=="DI"),]
sequences$SAMPLE <- gsub("Endo-AGO-IP-293T-","",sequences$SAMPLE)
sequences$SAMPLE <- gsub("_L001_R1_001.fastq_ready","",sequences$SAMPLE)
sequences$SAMPLE <- gsub("_S","",sequences$SAMPLE)

sequences$SEQ_TAIL <- as.character(sequences$SEQ_TAIL)
sequences$RATIOL <- as.numeric(sequences$RATIO)
sequences$Tscore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("T","",sequences$SEQ_TAIL))
sequences$Ascore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("A","",sequences$SEQ_TAIL))
sequences$Gscore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("G","",sequences$SEQ_TAIL))
sequences$Cscore <- nchar(sequences$SEQ_TAIL)-nchar(gsub("C","",sequences$SEQ_TAIL))


global2 <- global[,c(1:2)]
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
total_mir$CPM_avg <- total_mir$CPM_avg/9
total_mir <- total_mir[order(-total_mir$CPM_avg),]
total_mir <- unique(total_mir[c(1:200),])

sequences <- sequences[ (sequences$MIRNA %in% total_mir$MIRNA),]
global <- global[ (global$MIRNA %in% total_mir$MIRNA),]
global2 <- global2[ (global2$MIRNA %in% total_mir$MIRNA),]
global2$percent_tailNT <- global2$CPM_tailed/global2$CPM_total*100
global2$percent_trimmed <- global2$CPM_trimmed/global2$CPM_total*100

#miRNA_selected <- "hsa-miR-92a-3p-1-2"
#miRNA_selected <- "hsa-miR-20a-5p"
#miRNA_selected <- "hsa-miR-7-5p-1-2-3"
#miRNA_selected <- "hsa-miR-103a-3p-1-2"
#miRNA_selected <- "hsa-miR-221-3p"
#miRNA_selected <- "hsa-miR-148a-3p"
#miRNA_selected <- "hsa-miR-25-3p"
#miRNA_selected <- "hsa-let-7a-5p-1-2-3"
#miRNA_selected <- "hsa-let-7f-5p-1-2"
#miRNA_selected <- "hsa-miR-378a-3p"
#miRNA_selected <- "hsa-miR-30d-5p"
#miRNA_selected <- "hsa-miR-10b-5p"
#miRNA_selected <- "hsa-miR-99b-5p"
#miRNA_selected <- "hsa-miR-222-3p"
#miRNA_selected <- "hsa-miR-182-5p"
#miRNA_selected <- "hsa-miR-26a-5p-1-2"
#miRNA_selected <- "hsa-miR-769-5p"
#miRNA_selected <- "hsa-miR-30a-5p"
#miRNA_selected <- "hsa-miR-186-5p"
#miRNA_selected <- "hsa-miR-30e-5p"

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
setwd("~/Desktop/P6 Tiejuan/00_DIS3L2_KO_AG2017/new analysis DIS3L2_KO/individual Tscore/")
write.table(miRNA_selected_T_scores, paste0("T_score_",miRNA_selected,".tsv"), sep="\t", append = FALSE, row.names = F)



miRNA_selected <- "hsa-miR-7-5p-1-2-3"
#miRNA_selected <- "hsa-miR-769-5p"
#miRNA_selected <- "hsa-miR-222-3p"
miRNA_selected_group <- sequences[which(sequences$MIRNA==miRNA_selected & sequences$LEN_TAIL>=2),]
miRNA_selected_T_scores <- aggregate(miRNA_selected_group$CPM, by=list(SAMPLE=miRNA_selected_group$SAMPLE,Tscore=miRNA_selected_group$Gscore), FUN=sum)
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
setwd("~/Desktop/P6 Tiejuan/00_DIS3L2_KO_AG2017/new analysis DIS3L2_KO/individual Tscore/")
write.table(miRNA_selected_T_scores, paste0("G_score_",miRNA_selected,".tsv"), sep="\t", append = FALSE, row.names = F)


