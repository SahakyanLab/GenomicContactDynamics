################################################################################
# Re-calculate k-mer counts per unique contacting bin (BINKMER.MX) using an
# already masked genome (e.g. when treating repeats as mask and using
# repeat-masked genome). Bins with missing sequence (N) populated with NAs unless 
# masking character is N.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/X2_MaskingFeatures")
    data.dir = paste0(home.dir, "/Database")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/12_MaskingFeatures")
    data.dir = paste0(home.dir, "/Database")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
binkmer.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_allbins")
PATH.genome = paste0(data.dir, "/human_genome_repeatmasked_37.73")
out.dir = paste0(wk.dir, "/out_binkmer_repeatAsMask")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"))
kmer.len = 7L
# Unique HiC bins per chr ( getKmerCountsPerInterval() )
nCPU = 2L 
maskingChar = "N"
genome.prefix = "Homo_sapiens.GRCh37.73.dna_rm.chromosome."
fastafile.ending = ".fa"
affix = "_hg19_rm"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(gtools)
library(Biostrings)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 
source(paste0(lib.TrantoRextr, "/GEN_getKmers.R"))                 
source(paste0(lib, "/getMaskedKmerCount/getKmerCountsPerIntervalMASK.R"))  
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){
  
  chr.id <- strsplit(x=chr, split="chr")[[1]][2]
  
  loadGenome(PATH.genome=PATH.genome, genome.prefix=genome.prefix,
             fastafile.ending=fastafile.ending, chr.id=chr.id, 
             silent=FALSE, split=FALSE, remove.other.loads=TRUE)
  
  genome.filename <- paste0(genome.prefix, chr.id, fastafile.ending)
  
  #
  load(paste0(binkmer.dir, "/", chr, "_BinKmer7.RData"))
  ubins <- BINKMER.MX[,"bins"]
  bin.end <- BINKMER.MX[,"endpos"]
  bin.start <- BINKMER.MX[,"startpos"]
  
  rm(BINKMER.MX); gc()
  
  eval(parse(text=paste0(
    "BINKMER.MX <- getKmerCountsPerInterval(chr=chr, sequence=", genome.filename, "$seq, startpos=bin.start, endpos=bin.end, K=kmer.len, PATH.genome=PATH.genome, genome.prefix=genome.prefix, fastafile.ending=fastafile.ending, nCPU=nCPU, maskingChar=maskingChar)"
  )
  ))
  
  BINKMER.MX <- cbind(bins=ubins, startpos=bin.start, endpos=bin.end, BINKMER.MX)
  rm(ubins, bin.end, bin.start); gc()
  
  save(BINKMER.MX, file=paste0(out.dir, "/", chr, "_BinKmer", kmer.len, affix, ".RData"))
  
}

# rm(list=ls()); gc()
