################################################################################
# Measure complementarity in terms of matching short k-mers or 
# long span alignment
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    data.dir =  "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# both
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
out.dir = paste0(wk.dir, "/out_constraints")
persist.dir = paste0(wk.dir, "/out_features")
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
# align
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
# kmer
gfreepar.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binkmer.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_divLen_all")
### OTHER SETTINGS #############################################################
# both
gcb = "min2Mb"
chr.v = "chr18" #paste0("chr", c(22:1, "X"))
bin.len = 4e4L
kmer.len = 7L
type = "kmer" # "kmer" | "align"
# For type=align, nCPU based on number of chunks
# For type=kmer, nCPU based on number of contacts, ~30G for chr1
# chr21 - align - 368511 good contacts - 30G - 2 days
nCPU = 2L # chr1 - 4L, chr21 - 2L
ct = "hg19"
# align
affix.persist = "_ijShuffled"
affix.binkmer = ""
numChunks = 30L # chr1 - 32L, chr21 - 2L
gfreeparfile = paste0(gfreepar.dir, "/Gfree_", kmer.len, "mer.par")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(compiler)
library(Biostrings)
source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 
source(paste0(lib.TrantoRextr, "/GEN_getGenomicSeq.R"))
source(paste0(lib, "/getComplementarity.R"))
source(paste0(lib, "/getComplementarityScore.R"))

if(type=="align"){
  library(Rcpp)
  sourceCpp(file=paste0(lib, "/edlibNW.cpp"))
  source(paste0(lib, "/getComplementarityAlign.R"))
} else if(type=="kmer"){
  source(paste0(lib, "/getComplementarityKmer.R"))
} else {
  stop("Invalid input for type.")
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
getComplementarity(
  # both
  lib.TrantoRextr=lib.TrantoRextr,
  out.dir=out.dir,
  persist.dir=persist.dir,
  # File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
  chrLenfile=chrLenfile,
  # align
  genome.dir=genome.dir,
  # kmer
  gfreepar.dir=gfreepar.dir,
  binkmer.dir=binkmer.dir,
  # both
  gcb=gcb,
  chr.v=chr.v,
  bin.len=bin.len,
  kmer.len=kmer.len,
  type=type, 
  nCPU=nCPU,
  allij=FALSE,
  ct=ct,
  affix.persist=affix.persist,
  affix.binkmer=affix.binkmer, 
  affix.out=affix.persist, # paste0(chr, "_Hyb", kmer.len, "_", gcb, affix.out)
  # align
  numChunks=numChunks, 
  # kmer
  gfreeparfile=gfreeparfile 
)

# rm(list=ls())
