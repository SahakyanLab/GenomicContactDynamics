################################################################################
# Count k-mers for all bins of a chr
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Complementarity"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
out.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer3_divLen_all")
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"), sep="")
bin.len = 40000
kmer.len = 3
nCPU = 2 # ~4G
genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome."
fastafile.ending = ".fa"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 
source(paste0(lib.TrantoRextr, "/GEN_getKmers.R"))  
source(paste0(lib.TrantoRextr, "/GEN_getGenomicSeq.R"))  
source(paste0(wk.dir, "/lib/getKmerCountsPerInterval_new.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=TRUE, header=TRUE)

for(chr in chr.v){
  
  tot.bin <- ceiling(chrLen.df[chrLen.df$chromosome==chr, "length.bp"]/bin.len)
  
  if( tot.bin==chrLen.df[chrLen.df$chromosome==chr, "bins.40kb"] ){
    
    bin.end <- (1:tot.bin)*bin.len
    bin.start <- bin.end-bin.len+1
    
    BINKMER.MX <- getKmerCountsPerInterval(
      chr           = chr,
      startpos      = bin.start, # vector of starting positions
      endpos        = bin.end,   # vector of ending positions
      K             = kmer.len,  # k in the k-mer
      PATH.genome   = genome.dir,
      genome.prefix = genome.prefix,
      fastafile.ending = fastafile.ending,
      silent        = FALSE,
      nCPU          = nCPU
    )
    
    if( tot.bin!=nrow(BINKMER.MX) ){
      stop("Not all bins in BINKMER.MX.")
    }
    
  }
  
  genome.filename <- paste0(genome.prefix, strsplit(x=chr, split="chr")[[1]][2],
                            fastafile.ending)
  # Remove genome
  eval(parse(text=paste0(
    "rm(", genome.filename, "); gc()"
  )))

  # Change last position in bin end to actual end of chromosome
  bin.end.len <- length(bin.end)
  bin.end[bin.end.len] <- chrLen.df[chrLen.df$chromosome==chr, "length.bp"]
  
  BINKMER.MX <- cbind(bins=1:tot.bin, startpos=bin.start, endpos=bin.end, BINKMER.MX)
  
  # Check if the last bin has the correct length
  if(BINKMER.MX[bin.end.len,"numUMChar"]!=(BINKMER.MX[bin.end.len,"endpos"]-BINKMER.MX[bin.end.len,"startpos"]+1)){
    stop("Wrong length of last bin.")
  }
  
  rm(bin.end, bin.start, tot.bin, genome.filename)
  
  save(BINKMER.MX, file=paste0(out.dir, "/", chr, "_BinKmer", kmer.len, ".RData"))
  
  print(paste0(chr, " done!"), quote=FALSE)
  
  rm(BINKMER.MX); gc()
  
}

# rm(list=ls())
