################################################################################
# Function to make kmer counts per bin (in BINKMER.MX) account for both
# strands
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(foreach)
# library(doParallel)
# library(itertools)
# library(Biostrings)
# library(compiler)
### FUNCTION ###################################################################
makeKmerStrandInvar <- function(
  mx = BINKMER.MX,
  nCPU = 2L
){
  
  # Initialise the nCPU-driven (number of CPUs) setting for do/dopar
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  # Unique bins  
  ubin.v <- mx[,"bins"]
  ubin.v.len <- length(ubin.v)
  print(ubin.v.len)
  
  kmer.vec <- dimnames(mx)[[2]][-c(1:3)]
  kmer.vec.revcomp <- as.vector( reverseComplement(DNAStringSet(kmer.vec)) )
  kmer.revcomp.ind <- match(kmer.vec, kmer.vec.revcomp)
  rm(kmer.vec, kmer.vec.revcomp); gc()
  
  toExport <- c("mx", "kmer.revcomp.ind", "ubin.v")
  
  #### PARALLEL EXECUTION #########
  si.mx <- foreach(itr=isplitVector(1:ubin.v.len, chunks=nCPU),
                           .combine="rbind", .inorder=TRUE,
                           .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      bin <- ubin.v[i]
      bin.row <- mx[mx[,"bins"]==bin,-c(1:3)]
      # kmerBin + kmerBinrc to account for both strands
      kmerBinAll <- bin.row + bin.row[kmer.revcomp.ind]
    })
    return(do.call("rbind", chunk))
  }
  ### END OF PARALLEL EXECUTION ###
 
  return(si.mx) 
}
################################################################################
makeKmerStrandInvar <- cmpfun(makeKmerStrandInvar, options=list(suppressUndefined=TRUE))
################################################################################
