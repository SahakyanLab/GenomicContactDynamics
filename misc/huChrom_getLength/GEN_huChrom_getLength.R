#loaded human chromosomes fasta files to get the length (in bp) and 
#to partition into 40-kb bins 
###################################################################################################

#dependent functions
#readfasta script from alex
source("/Users/liezeltamon/OneDrive\ -\ Nexus365/ProjectBoard/HiC_general/HiC_Human21/TrantoRextr/GEN_readfasta.R")
#readLinesFast script from alex
source("/Users/liezeltamon/OneDrive\ -\ Nexus365/ProjectBoard/HiC_general/HiC_Human21/TrantoRextr/UTIL_readLinesFast.R")

#load GRCh37/hg19 genome using GEN_loadGenome script from Alex
#The directory of human chromosome fasta files 
huChromosome.dir  = "/Users/liezeltamon/OneDrive\ -\ Nexus365/ProjectBoard/HiC_general/human_genome_unmasked_37.73"

chrom.list       <- list()
chromLength.list <- list()
#chromHeader.list <- list()
for(i in c(1:22, "X", "Y", "MT")){
  chrom.list <- readfasta(paste0(huChromosome.dir, "/Homo_sapiens.GRCh37.73.dna.chromosome.", i, ".fa"))
  chromLength.list[[paste0("chr", i)]] <- chrom.list$length
  #chromHeader.list[[paste0("chr", i)]] <- chrom.list$header
}

#check sum of chromosome bp
genome.size <- sum(unlist(chromLength.list))

#instead of files, save as .RData
#makes a file of human chromosome info (length)
info.table <- matrix(unlist(chromLength.list, use.names=TRUE), byrow=TRUE)
colnames(info.table) <- "Length.bp"
write.table(info.table, file="/Users/liezeltamon/OneDrive\ -\ Nexus365/ProjectBoard/HiC_general/human_genome_unmasked_37.73/Homo_sapiens_GRCh37_73_chromosome_info1",
            sep="\t")
