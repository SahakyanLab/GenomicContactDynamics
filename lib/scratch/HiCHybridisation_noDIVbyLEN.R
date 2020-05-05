################################################################################
# Function to calculate kmer-based discordance score of HiC contacts for each
# chromosome
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(Biostrings)
# library(foreach)
# library(doParallel)
# library(itertools)
# library(compiler)
### FUNCTION ###################################################################
HiCHybridisation <- function(
  out.dir = "/dir",
  out.id = "_ij_shuffled",
  PERSIST.MX = "obj",
  BINKMER.MX = "obj",
  # gfree.par file
  gfreeparfile = paste0(gfreepar.dir, "/Gfree_", k, "mer.par"),
  gcb = "min2Mb",
  chr = "chr21",
  # Contacts per chromosome
  nCPU = 2L, 
  # Kmer length
  k = 7L,
  # Length of DNA/Hi-C resolution for Gfree scaling to per-nt
  scale = 40000L,
  saveOut = TRUE
){
  
  id <- paste0(chr, "_" , gcb, "_kmer", k)
  
  gfree <- read.table(file=gfreeparfile, header=TRUE, stringsAsFactors=FALSE)
  
  ############################
  kmer.vec <- dimnames(BINKMER.MX)[[2]][-c(1:3)]
  kmer.vec.revcomp <- as.vector(reverseComplement(DNAStringSet(kmer.vec)))
  kmer.revcomp.ind <- match(kmer.vec, kmer.vec.revcomp)
  rm(kmer.vec, kmer.vec.revcomp)
  ############################
  hits.len <- length(PERSIST.MX$hits[,1])
  ij.mx <- cbind(i=PERSIST.MX$hits[,"i"], j=PERSIST.MX$hits[,"j"])
  rm(PERSIST.MX); gc()
  ############################

  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  toExport <- c("BINKMER.MX", "ij.mx", "kmer.revcomp.ind", "gfree", "scale")
  
  #### PARALLEL EXECUTION #########
  
  HYB.MX <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                    .combine="cbind", .inorder=TRUE,
                    .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    hyb.data <- sapply(X=itr,
                       FUN=function(i){
                         
                         kmer.mx <- BINKMER.MX[ (BINKMER.MX[,"bins"] %in% ij.mx[i,]) , -c(1:3) ]
                         #kmerA <- kmer.mx[1,]
                         #kmerB <- kmer.mx[2,]
                         #kmerArc <- kmer.mx[1,kmer.revcomp.ind]
                         #kmerBrc <- kmer.mx[2,kmer.revcomp.ind]
                         kmerAall <- (kmer.mx[1,] + kmer.mx[1,kmer.revcomp.ind]) #(kmerA + kmerArc)
                         kmerBall <- (kmer.mx[2,] + kmer.mx[2,kmer.revcomp.ind]) #(kmerB + kmerBrc)
                         difference <- kmerAall - kmerBall
                         return( c(
                           # everything is counted twice, hence 2 in the division
                           sum(pmin(kmerAall, kmerBall)*gfree[,2])/(2*scale),
                           # everything is twice in frequency
                           #hist(difference, breaks=100, col="navy")
                           sd( difference ),
                           #mean( difference ) # mean is alws 0; double checked numerically
                           sum(abs(difference))
                         ) )
                       }, simplify=FALSE) # USE.NAMES=FALSE
    
    #return(hyb.data)
    return( do.call("cbind", hyb.data) )
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  print("HYB.MX generated!")
  
  ############################
  
  ## Add Cp 
  dimnames(HYB.MX)[[1]] <- c("Gfree", "Discordance", "Sumabsdif")
  #dimnames(HYB.MX)[[1]] <- c("Gfree", "Discordance", "Sumabsdif", "Cp")
  
  if(saveOut){
    save(HYB.MX, file=paste0(out.dir, "/", chr, "_Hyb", k, "_", gcb, out.id, ".RData"))
  }
  
  return(HYB.MX)

}
################################################################################
HiCHybridisation <- cmpfun(HiCHybridisation, options=list(suppressUndefined=TRUE))
################################################################################

# rm(list=ls())
 