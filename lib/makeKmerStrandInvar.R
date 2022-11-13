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
  nCPU = 1L,
  removeRevcomp = F
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

  kmer.vec <- setdiff( dimnames(mx)[[2]], c("bins", "startpos", "endpos", "numUMChar") ) #kmer.vec <- dimnames(mx)[[2]][-c(1:3)]
  kmer.vec.revcomp <- as.vector( reverseComplement(DNAStringSet(kmer.vec)) )
  kmer.revcomp.ind <- match(kmer.vec, table=kmer.vec.revcomp)
  rm(kmer.vec.revcomp)
  
  toExport <- c("ubin.v", "mx", "kmer.vec", "kmer.revcomp.ind")
  
  #### PARALLEL EXECUTION #########
  si.mx <- foreach(itr=isplitVector(1:ubin.v.len, chunks=nCPU), .combine="rbind", .inorder=TRUE,
                   .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      bin <- ubin.v[i]
      bin.row <- mx[ mx[,"bins"] == bin, kmer.vec ] #bin.row <- mx[ mx[,"bins"]==bin, -c(1:3) ]
      kmerBinAll <- bin.row + bin.row[kmer.revcomp.ind] # kmerBin + kmerBinrc to account for both strands
    })
    return(do.call("rbind", chunk))
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  mx[,kmer.vec] <- si.mx[,kmer.vec]
  rm(si.mx)
  
  if(removeRevcomp){
    
    ind.mx <- cbind(1:length(kmer.revcomp.ind), kmer.revcomp.ind)
    is_switch <- (ind.mx[,1] - ind.mx[,2]) > 0
    ind.mx[is_switch,] <- ind.mx[is_switch,2:1]
    kmer.vec.noRC <- kmer.vec[ !duplicated(ind.mx, fromLast=T) ]
    mx <- mx[, !dimnames(mx)[[2]] %in% kmer.vec.noRC ]
    
    # Checks 
    kmer.num.noRC <- (4^unique(nchar(kmer.vec))) / 2
    kmer.vec.noRC.rc <- as.vector( reverseComplement(DNAStringSet(kmer.vec.noRC)) )
    is_wrong <- c(any(duplicated(dimnames(mx)[[2]])), 
                  length(kmer.vec.noRC) != kmer.num.noRC,
                  length(intersect(kmer.vec.noRC, kmer.vec.noRC.rc)) != 0
                  )
    if( any(is_wrong) ){
      rm(mx)
      stop("makeKmerStrandInvar(): Unexpected behaviours.")
    }
  
  }
  
  return(mx) 
  
}
################################################################################
makeKmerStrandInvar <- cmpfun(makeKmerStrandInvar, options=list(suppressUndefined=TRUE))
################################################################################
