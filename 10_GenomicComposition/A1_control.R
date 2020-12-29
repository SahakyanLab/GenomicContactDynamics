################################################################################
# Make control matrix containing k-mer fraction of all unique HiC bins (strand 
# invariant) using BINKMER.MX.
# Output: KMERfrBIN.HiCAll (k-mer X ubins)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/9_GenomicComposition"
    data.dir = "/Users/ltamon/Database/"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/9_GenomicComposition"
    data.dir = "/t1-data/user/ltamon/Database/"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
binkmer.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_orig_persist")
out.dir = paste0(wk.dir, "/out_control")
### OTHER SETTINGS #############################################################
chr.v = paste("chr", c(17, 22:18, 16:1, "X"), sep="")
gcb = "min05Mb"
kmer.len = 7
HiC.res = 4e4L
# Unique bins per chromosome
nCPU = 3L 
# Identifier of BINKMER.MX (could be the shuffled set so I added this)
affix = ""
saveStrandInvarMx = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(Biostrings)
library(compiler)
source(paste0(wk.dir, "/lib/makeKmerStrandInvar.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)

# Multiply by 2 because your counting kmers of both strands
totcountPBin <- ( HiC.res-(kmer.len-1L) )*2L

toExport <- c("binkmer.dir", "kmer.len", "gcb", "totcountPBin")

#### PARALLEL EXECUTION #########

KMERfrBIN.HiCAll <- foreach(chr=chr.v, .combine="cbind"
                            
) %do% {
  
  # Load BINKMER.MX 
  data.nme <- load(paste0(binkmer.dir, "/", chr, "_BinKmer", kmer.len,"_", gcb, 
                          affix, ".RData"))
  # Because old naming specifies kmer length i.e. BINKMER7.MX
  eval(parse( text=paste0("BINKMER.MX <-", data.nme) ) )
  if(data.nme!="BINKMER.MX"){
    eval(parse( text=paste0("rm(", data.nme, ")") ) )
  }
  
  if(chr=="chr17"){
    ind <- which( !rowSums(BINKMER.MX[,-(1:3)], na.rm=TRUE)%in%c(totcountPBin/2, 0) )
    print( paste0(chr, ": Removed ", ind, "th bin"), quote=FALSE )
    BINKMER.MX <- BINKMER.MX[-(ind),]
  }
  
  # Convert kmer counts to account for both strand
  BINKSTRINV.MX <- makeKmerStrandInvar(mx=BINKMER.MX, nCPU=nCPU)
  
  # Added values after 0 are for chr17
  if( any( !rowSums(BINKSTRINV.MX, na.rm=TRUE)%in%c(totcountPBin, 0) ) | nrow(BINKSTRINV.MX)!=nrow(BINKMER.MX) ){
    stop("Checkpoint 1.")
  } else {
    
    if(saveStrandInvarMx==TRUE){
      BINKSTRINV.MX <- cbind(BINKMER.MX[,1:3], BINKSTRINV.MX)
      save(BINKSTRINV.MX, file=paste0(binkmer.dir, "/", chr, "_BinKStrInv", 
                                      kmer.len, "_", gcb, affix, ".RData"))
      BINKSTRINV.MX <- BINKSTRINV.MX[,-(1:3)]
    }
    
  }
  
  ubins <- BINKMER.MX[,"bins"]
  rm(BINKMER.MX); gc()
  
  # Convert kmer counts to fraction of kmer count per bin
  BINKSTRINV.MX <- BINKSTRINV.MX/rowSums(BINKSTRINV.MX, na.rm=FALSE)
  chr.num <- gsub(x=chr, pattern="chr", replacement="", fixed=TRUE)
  rownames(BINKSTRINV.MX) <- paste(chr.num, ubins, sep="_")
  rm(ubins); gc()
  
  # Transpose to kmer X bin
  BINKSTRINV.MX <- t(BINKSTRINV.MX)
  
  print(paste0(chr, " done!"))
  
  return(BINKSTRINV.MX)
  
} 

### END OF PARALLEL EXECUTION ###

save(KMERfrBIN.HiCAll, file=paste0(out.dir, "/HiCAll_KMERfrBINPCP", kmer.len,
                                   "_", gcb, ".RData"))

# rm(list=ls()) 