################################################################################
# Function to calculate kmer-based complementarity score of Hi-C contacting bins
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(Biostrings)
# library(foreach)
# library(doParallel)
# library(itertools)
# library(compiler)
### FUNCTION ###################################################################
getComplementarityKmer <- function(
  out.dir = "/dir",
  out.id = paste0(chr, "_Hyb", k, "_", gcb),
  # Should be 2 columns, i and j
  contact.mx = PERSIST.MX$hits[,c("i", "j")],
  # Contacts per chromosome
  nCPU = 2L, 
  binkmer.mx = BINKMER.MX,
  # gfree.par file
  gfreeparfile = paste0(gfreepar.dir, "/Gfree_", k, "mer.par"),
  # Kmer length
  kmerlength = 7L,
  # Length of DNA/Hi-C resolution for Gfree scaling to per-nt
  binlength = 40000L,
  saveOut = TRUE
){
 
  gfree <- read.table(file=gfreeparfile, header=TRUE, stringsAsFactors=FALSE)
  
  ############################
  kmer.vec <- dimnames(binkmer.mx)[[2]][-c(1:4)]
  kmer.vec.revcomp <- as.vector(reverseComplement(DNAStringSet(kmer.vec)))
  kmer.revcomp.ind <- match(kmer.vec, kmer.vec.revcomp)
  rm(kmer.vec, kmer.vec.revcomp)
  ############################
  contact.mx <- data.matrix(contact.mx)
  hits.len <- length(contact.mx[,1])
  ############################

  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  toExport <- c("binkmer.mx", "contact.mx", "kmer.revcomp.ind", "gfree", "scale")
  
  #### PARALLEL EXECUTION #########
  
  HYB.MX <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                    .combine="cbind", .inorder=TRUE,
                    .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    hyb.data <- sapply(X=itr,
                       FUN=function(i){
                         
                         kmer.mx <- binkmer.mx[ (binkmer.mx[,"bins"] %in% contact.mx[i,]) , -c(1:4) ]
                         
                         if( length(kmer.mx[,1])!=2 ){
                           return( c(NA, NA, NA) )
                         } else {
                           #kmerA <- kmer.mx[1,]
                           #kmerB <- kmer.mx[2,]
                           #kmerArc <- kmer.mx[1,kmer.revcomp.ind]
                           #kmerBrc <- kmer.mx[2,kmer.revcomp.ind]
                           kmerAall <- (kmer.mx[1,] + kmer.mx[1,kmer.revcomp.ind]) #(kmerA + kmerArc)
                           kmerBall <- (kmer.mx[2,] + kmer.mx[2,kmer.revcomp.ind]) #(kmerB + kmerBrc)
                           difference <- kmerAall - kmerBall
                           return( c(
                             # everything is counted twice, hence 2 in the division
                             #sum(pmin(kmerAall, kmerBall)*gfree[,2])/(2*binlength),
                             sum(pmin(kmerAall, kmerBall)*gfree[,2])/2, # kmer counts already normalised by binlength
                             # everything is twice in frequency
                             #hist(difference, breaks=100, col="navy")
                             sd( difference ),
                             #mean( difference ) # mean is alws 0; double checked numerically
                             # negate to get value directly proportional to complementarity
                             -( sum(abs(difference)) )
                           ) )
                         }
                         
                       }, simplify=FALSE) # USE.NAMES=FALSE
    
    #return(hyb.data)
    return( do.call("cbind", hyb.data) )
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  print("HYB.MX generated!")
  
  ############################
 
  dimnames(HYB.MX)[[1]] <- c("Gfree", "sdDifference", "NegSumAbsDiff")
 
  if(saveOut){
    save(HYB.MX, file=paste0(out.dir, "/", out.id, ".RData"))
  }
  
  return(HYB.MX["NegSumAbsDiff",])

}
################################################################################
getComplementarityKmer  <- cmpfun(getComplementarityKmer, options=list(suppressUndefined=TRUE))
################################################################################

# rm(list=ls())
 