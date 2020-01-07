setwd("~/Desktop/P6 Tiejuan/TDMD_target_pred/")
HEK293T <- read.table(file="HEK293T_WT_genes_FPKM.tsv", header = T, fill = TRUE, sep =  "\t")
HEK293T <- unique(HEK293T[,c(2,5)])
HEK293T <- aggregate(HEK293T$FPKM, by=list(SYMBOL=HEK293T$external_gene_name), max)

miR7 <- read.table(file="mir_7_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
#cyrano <- readLines("cyrano_cDNA.txt")
#cyrano <- gsub("T", "U", cyrano)
setwd("~/Desktop/P6 Tiejuan/TDMD_target_pred/energy_pred/")
getDuplexEnergy<-function(ref_seq,target_seq){
  input_data<-c(paste('>ref', sep=''), as.character(ref_seq), paste('>test', sep=''), as.character(target_seq))
  input_file<-'temp-miRNA-seq.fa'
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  energy<-system('RNAduplex < temp-miRNA-seq.fa', intern=T)[3]
  energy <- substr(energy,nchar(energy)-10,nchar(energy))
  energy<-gsub("\\(", "\\*", energy)
  energy<-gsub("\\)", "\\*", energy)
  energy<-strsplit(energy, "*", fixed=T)[[1]][2]
  energy<-as.numeric(energy)
  return(energy)
}
getDuplexEnergy("UGGAAGACUAGUGAUUUUGUUGUU",cyrano)

setwd("~/Desktop/P6 Tiejuan/TDMD_target_pred/")
miR7_FPKM <- merge(miR7,HEK293T, by="SYMBOL")
miR7_FPKM <- miR7_FPKM[which(miR7_FPKM$x>0),]
plot(log10(miR7_FPKM$x), miR7_FPKM$energy*-1)
write.table(x = miR7_FPKM, file = "miR7_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

miR222 <- read.table(file="mir_222_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
miR222_FPKM <- merge(miR222,HEK293T, by="SYMBOL")
miR222_FPKM <- miR222_FPKM[which(miR222_FPKM$x>0),]
plot(log10(miR222_FPKM$x), miR222_FPKM$energy*-1)
write.table(x = miR222_FPKM, file = "miR222_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

miR769 <- read.table(file="mir_769_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
miR769_FPKM <- merge(miR769,HEK293T, by="SYMBOL")
miR769_FPKM <- miR769_FPKM[which(miR769_FPKM$x>0),]
plot(log10(miR769_FPKM$x), miR769_FPKM$energy*-1)
write.table(x = miR769_FPKM, file = "miR769_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

miR21 <- read.table(file="mir_21_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
miR21_FPKM <- merge(miR21,HEK293T, by="SYMBOL")
miR21_FPKM <- miR21_FPKM[which(miR21_FPKM$x>0),]
plot(log10(miR21_FPKM$x), miR21_FPKM$energy*-1)
write.table(x = miR21_FPKM, file = "mi21_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

miR92 <- read.table(file="mir_92a_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
miR92_FPKM <- merge(miR92,HEK293T, by="SYMBOL")
miR92_FPKM <- miR92_FPKM[which(miR92_FPKM$x>0),]
plot(log10(miR92_FPKM$x), miR92_FPKM$energy*-1)
write.table(x = miR92_FPKM, file = "mi92_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

miR10a <- read.table(file="mir_10a_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
miR10a_FPKM <- merge(miR10a,HEK293T, by="SYMBOL")
miR10a_FPKM <- miR10a_FPKM[which(miR10a_FPKM$x>0),]
plot(log10(miR10a_FPKM$x), miR10a_FPKM$energy*-1)
write.table(x = miR10a_FPKM, file = "mi10a_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

miR148a <- read.table(file="mir_148a_7mer.tsv", header = T, fill = TRUE, sep =  "\t")
miR148a_FPKM <- merge(miR148a,HEK293T, by="SYMBOL")
miR148a_FPKM <- miR148a_FPKM[which(miR148a_FPKM$x>0),]
plot(log10(miR148a_FPKM$x), miR148a_FPKM$energy*-1)
write.table(x = miR148a_FPKM, file = "mi148a_FPKM_TDMD_targets.tsv",row.names = F,col.names = T, sep = "\t")

