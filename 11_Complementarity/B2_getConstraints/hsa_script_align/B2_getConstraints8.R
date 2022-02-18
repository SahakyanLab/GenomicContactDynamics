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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# both
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
out.dir = paste0(wk.dir, "/out_constraints")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
# align
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
# kmer
gfreepar.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binkmer.dir = paste0(persist.dir, "/binkmer_divLen_all")
### OTHER SETTINGS #############################################################
# both
gcb = "min2Mb"
chr.v = "chr8"
bin.len = 4e4L
kmer.len = 7L
type = "align" # "kmer" | "align"
# For type=align, nCPU based on number of chunks
# For type=kmer, nCPU based on number of contacts, ~30G for chr1
# chr21 - align - 368511 good contacts - 30G - 2 days
nCPU = 2L # chr1 - 4L, chr21 - 2L
allij = TRUE
# align
numChunks = 30L # chr1 - 32L, chr21 - 2L
# kmer
affix = "" # paste0(chr, "_Hyb", kmer.len, "_", gcb, affix) 
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
source(paste0(wk.dir, "/lib/getComplementarityScore.R"))

if(type=="align"){
  
  library(Rcpp)
  sourceCpp(file=paste0(wk.dir, "/lib/edlibNW.cpp"))
  source(paste0(wk.dir, "/lib/getComplementarityAlign.R"))
  
} else if(type=="kmer"){
  
  source(paste0(wk.dir, "/lib/getComplementarityKmer.R"))
  
} else {
  
  stop("Invalid input for type.")
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=TRUE, header=TRUE)

for(chr in chr.v){
  
  print(paste0(chr, "..."), quote=FALSE)
  
  start.time <- Sys.time()
  
  # Load from PERSIST.MX to get Cp data; rowname of contact in PERSIST.MX$hits 
  # corresponds to the index of that contact in a matrix of all ij contacts for a
  # given chr obtained using expand.grid() with self and duplicated pairs not yet 
  # removed. The rowname therefore can be used to assign the Cp value of some 
  # contacts in the all ij contact matrix.
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))

  # Make mx of contacts
  if(allij){
    
    print("All ij contacts, where i<j...", quote=FALSE)
    
    tot.bin <- ceiling(chrLen.df[chrLen.df$chromosome==chr, "length.bp"]/bin.len)
    ubins <- 1:tot.bin
    
    # Create matrix of all pairs 
    contact.mx <- data.matrix(expand.grid(ubins, ubins))
    rm(ubins, tot.bin)
    
    # Add Cp using the rownames of PERSIST.MX$hits
    contact.mx <- cbind(contact.mx, Cp=NA)
    contact.mx[as.numeric(rownames(PERSIST.MX$hits)), "Cp"] <- PERSIST.MX$ntis
    
    # To remove self pairs and duplicates like {1,2}(keep) and {2,1}
    contact.mx <- contact.mx[contact.mx[,2]-contact.mx[,1]>0,]
    
  } else {
    
    print("ij contacts from PERSIST.MX...", quote=FALSE)
    contact.mx <- cbind( data.matrix( PERSIST.MX$hits[, c("i", "j")] ),
                         Cp=PERSIST.MX$ntis )
    rownames(contact.mx) <- NULL
    
  }
  
  rm(PERSIST.MX); gc()
  
  if(type=="kmer"){
    # Load BINKMER.MX
    load(file=paste0(binkmer.dir, "/", chr, "_BinKmer", kmer.len, ".RData"))
  }
  
  # Calculate complementarity score
  CII.MX <- cbind(contact.mx,
                  getComplementarityScore(
                    
                    contact.mx=contact.mx[,1:2],
                    binlength=bin.len,
                    type=type,
                    # For type=align, nCPU based on number of chunks
                    # For type=kmer, nCPU based on number of contacts 
                    nCPU=nCPU,
                    
                    # Alignment-specific arguments
                    PATH.genome=genome.dir,
                    genome.prefix="Homo_sapiens.GRCh37.73.dna.chromosome.",
                    fastafile.ending=".fa",
                    silent=TRUE,
                    # Number of chunks to split the contacts for alignment
                    numChunks=numChunks,
                    
                    # K-mer-method-specific arguments
                    out.dir=out.dir,
                    out.id=paste0(chr, "_Hyb", kmer.len, "_", gcb, affix),
                    binkmer.mx=BINKMER.MX,
                    gfreeparfile=gfreeparfile ,
                    # Kmer length
                    kmerlength=kmer.len,
                    saveOut=TRUE
                    
                  )
  )
  rm(contact.mx); gc()
  
  dimnames(CII.MX)[[2]] <- c("i", "j", "Cp", "C||")
  save(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb, ".RData"))
  
  rm(CII.MX); gc()
  
  end.time <- Sys.time()
  print(end.time-start.time, quote=FALSE)
  
  print(paste0(chr, " done!"), quote=FALSE)
  
}

# rm(list=ls())
