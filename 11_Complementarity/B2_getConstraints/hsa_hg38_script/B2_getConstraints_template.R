################################################################################
# Measure complementarity in terms of matching short k-mers or 
# long span alignment
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/11_Complementarity")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Constraints")
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/11_Complementarity")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

# both
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
out.dir = paste0(wk.dir, "/out_constraints_hg38_GfreeSingleNorm")
persist.dir = NULL #paste0(data.dir, "/HiC_features_GSE87112_RAWpc")  # NULL if not relevant, populate Cp with NAs

# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh38_chr_info.txt")
# align
genome.dir = paste0(data.dir, "/human_genome_unmasked_38.108")
# kmer
gfreepar.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
#binkmer.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/12_MaskingFeatures/out_binkmer_maskingThreshold")
binkmer.dir = paste0(data.dir, "/HiC_features_human_hg38.108/binkmer_allbins")
### OTHER SETTINGS #############################################################
# both
gcb = "min2Mb"
chr.v = "chrarr1.repl"
bin.len = 50000 #2000 #4e4L
kmer.len = 7
type = "arr2.repl" # "kmer" | "align"
# For type=align, nCPU based on number of chunks
# For type=kmer, nCPU based on number of contacts, ~30G for chr1
# chr21 - align - 368511 good contacts - 30G - 2 days
nCPU = 2 # chr1 - 4L (~40G), chr22 - 2L (~4G)
allij = TRUE
# align
numChunks = 2 # human chr1 - 32L, chr21 - 2L
gfreeparfile = paste0(gfreepar.dir, "/Gfree_", kmer.len, "mer.par")
genome.prefix = "Homo_sapiens.GRCh38.dna.chromosome." #"Homo_sapiens.GRCh37.73.dna.chromosome." 
fastafile.ending = ".fa"
affix.binkmer = ""
affix.persist = ""
affix.out = ""
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(compiler)
library(gtools)
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
  allij=allij,
  ct=NULL,
  affix.persist=affix.persist,
  affix.binkmer=affix.binkmer, # paste0(chr, "_Hyb", kmer.len, "_", gcb, affix)
  affix.out=affix.out, 
  genome.prefix=genome.prefix,
  fastafile.ending=fastafile.ending,
  # align
  numChunks=numChunks, 
  # kmer
  gfreeparfile=gfreeparfile 
)

# rm(list=ls()); gc()
