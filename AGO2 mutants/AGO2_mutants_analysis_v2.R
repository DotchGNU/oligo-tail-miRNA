setwd("~/Desktop/P6 Tiejuan/02_F2L3 AGO2IP comparison/AGO2_mutants_TCGA/")
library(ggplot2)

manifest_data <- read.table(file="manifest_all_AGO2_mutants_TCGA2.txt", header = TRUE, fill = TRUE, sep =  "\t")
manifest_data$name2 <- gsub(manifest_data$name,pattern = "_mirna_gdc_realn.bam", replacement = "")
manifest_data <- manifest_data[which(manifest_data$sample_type=="Primary Tumor"),] #only_primary_tumors
variant_des <- unique(manifest_data[,c(26,22,23)])
variant_des <- variant_des[which(variant_des$Variant.Classification=="synonymous_variant" | variant_des$Variant.Classification=="missense_variant"),]
summary(variant_des$Variant.Classification)
summary(variant_des$AGO2_mut)
mut_to_domain <- unique(manifest_data[,c(23,24)])
mut_to_domain <- mut_to_domain[complete.cases(mut_to_domain), ]
mut_to_domain$Domain <- gsub(pattern = "N_terminal", replacement = "N", x = mut_to_domain$Domain)

isomir_sequences <- read.table(file="AGO2_mutants_TCGA_run2_5.isomir.sequence_info.tsv", header = TRUE, sep =  "\t")
most_abundant <- aggregate(isomir_sequences$READS, by=list(SAMPLE=isomir_sequences$SAMPLE, MIRNA=isomir_sequences$MIRNA), FUN=max)
colnames(most_abundant) <- c("SAMPLE","MIRNA","READS")
most_abundant <- most_abundant[which(most_abundant$READS>500),]
most_abundant <- merge(most_abundant, isomir_sequences, by=c("MIRNA","SAMPLE","READS"))
most_abundant$ratio2 <- 100-as.numeric(substr(most_abundant$RATIO, 1, (nchar(as.character(most_abundant$RATIO))-1)))
most_abundant$name2 <- gsub(most_abundant$SAMPLE,pattern = "_mirna_gdc_realn.converted.unpaired.fastq", replacement = "")
rm(isomir_sequences)

isomir_summary <- read.table(file="AGO2_mutants_TCGA_run2_5.isomir.tsv", header = TRUE, sep =  "\t")
MIR143_3P <- isomir_summary[which(isomir_summary$MIRNA=="hsa-miR-143-3p"),]
library(Rtsne) # Load package
set.seed(42) # Sets seed for reproducibility
MIR143_3P$name2 <- gsub(MIR143_3P$SAMPLE,pattern = "_mirna_gdc_realn.converted.unpaired.fastq", replacement = "")
MIR143_3P <- merge(MIR143_3P,manifest_data,by="name2")
MIR143_3P <- MIR143_3P[,c(29,31,8:17)]
MIR143_3P <- unique(MIR143_3P) # Remove duplicates
tsne_out <- Rtsne(as.matrix(MIR143_3P[,3:12])) # Run TSNE
plot(tsne_out$Y,col=MIR143_3P$primary_site,asp=1) # Plot the result

#mutation <- "P295L"
#mutation <- "R315M"
mutation <- "E299K"

avg_expression <- aggregate(most_abundant$CPM, by=list(most_abundant$MIRNA), FUN=mean)
n <- 150
colnames(avg_expression) <- c("miRNA","avg_CPM")
avg_expression <- avg_expression[order(-avg_expression$avg_CPM),]
row.names(avg_expression) <- NULL
avg_expression$synonimous_iso <- NA
avg_expression$missense_iso <- NA
avg_expression$N_iso <- NA
avg_expression$L1_iso <- NA
avg_expression$PAZ_iso <- NA
avg_expression$L2_iso <- NA
avg_expression$MID_iso <- NA
avg_expression$PIWI_iso <- NA
avg_expression$P295L <- NA
avg_expression <- avg_expression[c(1:n),]

i <- 2
for (i in 1:n) {
  miRNA <- as.character(avg_expression$miRNA[i])
  mir21 <- most_abundant[which(most_abundant$MIRNA==miRNA),]
  mir21 <- merge(mir21,variant_des, by="name2")
  isomir_index_synonimous <- mean(mir21[which(mir21$Variant.Classification=="synonymous_variant"),16])
  isomir_index_missense <- mean(mir21[which(mir21$Variant.Classification=="missense_variant"),16])
  avg_expression$synonimous_iso[i] <- isomir_index_synonimous
  avg_expression$missense_iso[i] <- isomir_index_missense
  mut_P295L <- mir21[which(mir21$AGO2_mut==mutation),]
  if(nrow(mut_P295L)>0){
    avg_expression$P295L[i] <- mir21[which(mir21$AGO2_mut==mutation),16]
  }
  mir21 <- merge(mir21,mut_to_domain, by="AGO2_mut")
  isomir_index_N <- mean(mir21[which(mir21$Domain=="N"),17])
  avg_expression$N_iso[i] <- isomir_index_N
  isomir_index_L1 <- mean(mir21[which(mir21$Domain=="L1"),17])
  avg_expression$L1_iso[i] <- isomir_index_L1
  isomir_index_PAZ <- mean(mir21[which(mir21$Domain=="PAZ"),17])
  avg_expression$PAZ_iso[i] <- isomir_index_PAZ
  isomir_index_L2 <- mean(mir21[which(mir21$Domain=="L2"),17])
  avg_expression$L2_iso[i] <- isomir_index_L2
  isomir_index_MID <- mean(mir21[which(mir21$Domain=="MID"),17])
  avg_expression$MID_iso[i] <- isomir_index_MID
  isomir_index_PIWI <- mean(mir21[which(mir21$Domain=="PIWI"),17])
  avg_expression$PIWI_iso[i] <- isomir_index_PIWI
  
  rm(miRNA,mir21,isomir_index_synonimous,isomir_index_missense,isomir_index_N,isomir_index_L1,isomir_index_PAZ,isomir_index_L2,isomir_index_MID,isomir_index_PIWI)
  print(i)
}

avg_expression <- avg_expression[which(avg_expression$synonimous_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$missense_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$N_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$L1_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$PAZ_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$L2_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$MID_iso!="NaN"),]
avg_expression <- avg_expression[which(avg_expression$PIWI_iso!="NaN"),]
avg_expression <- avg_expression[complete.cases(avg_expression), ]
avg_expression <- avg_expression[order(avg_expression$synonimous_iso),]

pval <- t.test(avg_expression$synonimous_iso, avg_expression$P295L, paired=TRUE)[3]
title <- paste0(mutation,"_p-val:",pval)
plot(avg_expression$synonimous_iso, avg_expression$P295L, main=title)
abline(a=1, b=1)
wilcox.test(avg_expression$synonimous_iso, avg_expression$P295L, paired=TRUE, alternative="less")
wilcox.test(avg_expression$missense_iso, avg_expression$P295L, paired=TRUE, alternative="less")


# t.test(avg_expression$synonimous_iso, avg_expression$missense_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$N_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$L1_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$PAZ_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$L2_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$MID_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$PIWI_iso, paired=TRUE, alternative="less")
# t.test(avg_expression$synonimous_iso, avg_expression$P295L, paired=TRUE, alternative="less")
# t.test(avg_expression$missense_iso, avg_expression$PAZ_iso, paired=TRUE, alternative="less")

# t.test(avg_expression$synonimous_iso, avg_expression$missense_iso, paired=TRUE)

# t.test(avg_expression$synonimous_iso, avg_expression$PAZ_iso, paired=TRUE)
# t.test(avg_expression$missense_iso, avg_expression$PAZ_iso, paired=TRUE)

name_file <- paste0("AGO2_missense_synonimous_",mutation,"_miRNA2.txt")
write.table(avg_expression, name_file, sep="\t", append = FALSE, row.names = FALSE)
