################################################################################
# R API to edlib, Lightweight, super fast C/C++ (& Python) library for sequence 
# alignment using edit (Levenshtein) distance
# https://github.com/Martinsos/edlib

# Dependencies
#- edlib.h - header file
#- edlib.cpp - edlib function source 
#- edlib.NW - Rcpp code calling edlib function
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
### OTHER SETTINGS #############################################################
HiC.res = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(Rcpp)
sourceCpp(paste0(wk.dir, "/lib/edlibNW.cpp"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
#load(file=paste0(wk.dir, "/sample_ij_seq.RData"))

q.v <- c("AACTGGC", "hello", "bag", "AACT", "telephone")
ql.v <- unname(sapply(X=q.v, FUN=nchar))
t.v <- c("AACT", "world!", "bag", "AACTGGC", "elephant")
tl.v <- unname(sapply(X=t.v, FUN=nchar))

edlibNW(query=q.v, target=t.v, 
        queryLength=ql.v, targetLength=tl.v)

#-----
edlibNW(query="AACTGGC", target="AACT", 
        queryLength=7, targetLength=4)

edlibNW(query="AACT", target="AACTGGC", 
        queryLength=4, targetLength=7)

edlibNW(query="telephone", target="elephant", 
        queryLength=9, targetLength=8)

contact.mx=contact.mx[,1:2];
binlength=bin.len;
type=type;
# For type=align, nCPU based on number of chunks
# For type=kmer, nCPU based on number of contacts 
nCPU=nCPU;

# Alignment-specific arguments
PATH.genome=genome.dir;
genome.prefix="Homo_sapiens.GRCh37.73.dna.chromosome.";
fastafile.ending=".fa";
silent=TRUE;
# Number of chunks to split the contacts for alignment
numChunks=numChunks;

# K-mer-method-specific arguments
out.dir=out.dir;
out.id=paste0(chr, "_Hyb", kmer.len, "_", gcb, affix);
binkmer.mx=BINKMER.MX;
gfreeparfile=gfreeparfile;
# Kmer length
kmerlength=kmer.len;
saveOut=TRUE

rm(list=ls())
