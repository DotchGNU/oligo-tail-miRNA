setwd("~/Desktop/P6 Tiejuan/TDMD_target_pred/")

human_UTR <- read.table(file="human_3UTR_v2.txt", header = FALSE, fill = TRUE, sep =  "\t")
colnames(human_UTR) <- c("SYMBOL","UTR")
human_UTR$UTR <- toupper(human_UTR$UTR)

#mir7 <- "UGGAAGACUAGUGAUUUUGUUGUU"
#mir21 <- "UAGCUUAUCAGACUGAUGUUGA"
#mir769 <- "UGAGACCUCUGGGUUCUGAGCU"
#mir222 <- "AGCUACAUCUGGCUACUGGGU"
mir92a <- "UAUUGCACUUGUCCCGGCCUGU"
#mir10a <- "UACCCUGUAGAUCCGAAUUUGUG"
#mir148a <- "UCAGUGCACUACAGAACUUUGU"
miRNA <- mir92a
miRNA <- gsub("U", "T", miRNA)

library(seqinr)
seed_mer <- 7
rev_compl <- toupper(c2s(rev(comp(s2c(miRNA)))))
seed <- substr(rev_compl,(nchar(rev_compl)-seed_mer),(nchar(rev_compl)-1))
seed <- gsub("T", "U", seed)

human_UTR$motif <- gregexpr(seed, human_UTR$UTR, perl = TRUE) #analysis only of UTRs with seed
human_UTR$energy <- NA
mer7<- human_UTR[which(human_UTR$motif!="-1"),]
#no_targets<- human_UTR[which(human_UTR$motif=="-1"),] #no_targets
targets <- mer7

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

i <- 1
for (i in 1:nrow(targets)) {
  test_miRNA <- gsub("T", "U", miRNA)
  test_UTR <- targets$UTR[i]
  energy <- getDuplexEnergy(test_miRNA,test_UTR)
  targets$energy[i] <- energy
  print(i/nrow(targets)*100)
  rm(test_miRNA,test_UTR,energy)
}

setwd("~/Desktop/P6 Tiejuan/TDMD_target_pred/")
targets$Nsites <- nchar(as.character(targets$motif))-nchar(gsub(",","",as.character(targets$motif)))+1

targets <- targets[,c(1,4,5)]
targets <- unique(targets)

hist(targets$Nsites)
hist(targets$energy)

write.table(x = targets, file = "mir_148a_7mer.tsv",row.names = F,col.names = T, sep = "\t")
#targets <- read.table("mir_769_no_mer.tsv", header = TRUE)
breaks <-  seq(-40, 10, by=0.5)
energy_dist = c(0,table(cut(targets$energy, breaks, right=FALSE))/nrow(targets)*100)
plot(breaks, energy_dist, type="n", xlab = "Dependency Score", ylab = "Distribution", ylim=c(0,15), xlim=c(-40,10))
lines(breaks, energy_dist)
final_dist <- rbind(breaks,energy_dist)
#write.table(x = final_dist, file = "mir_769_no_mer_dist.tsv",row.names = F,col.names = T, sep = "\t")

target_specific <- mer7[which(mer7$SYMBOL=="METTL8"),2]
getDuplexEnergy(mir92a,target_specific)
