################################################################################
# Check if last bins of chromosome have N's
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
out.dir = paste0(wk.dir, "/out_lastbinN")
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"), sep="")
bin.len = 40000
genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome."
fastafile.ending = ".fa"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 
source(paste0(lib.TrantoRextr, "/GEN_getGenomicSeq.R"))  
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrLen.df <- read.table(file=chrLenfile, as.is=TRUE, header=TRUE)
df <- data.frame(chr=chr.v, lastbin.len=NA, hasN=NA)

for(chr in chr.v){
  chr.len <- chrLen.df[chrLen.df$chromosome==chr, "length.bp"]
  last.bin <- ceiling(chr.len/bin.len)
  last.bin.start <- last.bin*bin.len-40000+1
  seq <- getGenomicSeq(PATH.genome=genome.dir,
                       genome.prefix=genome.prefix,
                       fastafile.ending=fastafile.ending,
                       chr.id=strsplit(chr,"chr")[[1]][2],
                       remove.other.loads=TRUE,
                       silent=FALSE, split=FALSE,
                       borders=c(last.bin.start, chr.len))
  if(length(grep("N", seq))==1){
    df[df$chr==chr,c("lastbin.len", "hasN")] <- c(chr.len-last.bin.start+1, "T")
  } else {
    df[df$chr==chr,c("lastbin.len", "hasN")] <- c(chr.len-last.bin.start+1, "F")
  }
  rm(seq, chr.len, last.bin, last.bin.start); gc()
  print(paste0(chr, " done!"), quote=FALSE)
}
write.table(df, file=paste0(out.dir, "/lastbinN"), quote=FALSE, sep="\t",
            col.names=TRUE, row.names=FALSE)

# rm(list=ls())
