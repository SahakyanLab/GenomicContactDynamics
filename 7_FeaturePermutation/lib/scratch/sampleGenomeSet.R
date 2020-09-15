################################################################################
# Sample genome set provided as set of ranges
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(regioneR)
### FUNCTION ###################################################################
sampleGenomeSet <- function(A, genome=genome, per.chromosome=TRUE){
  
  A <- regioneR::toGRanges(A)
  A.len <- length(A)
  genome <- regioneR::toGRanges(genome)
  g.len <- length(length(genome))
  
  if(per.chromosome){
    
    if( any(!seqnames(A)%in%seqnames(genome)) ){
      stop("Some chr in A not in genome.")
    } else {
      chr.v <- as.character( unique(seqnames(A)) )
      smpld <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
        chr.genome <- genome[seqnames(genome)==chr]
        chr.A.len <- sum(seqnames(A)==chr)
        chr.g.len <- length(chr.genome)
        ind.smpld <- sample(x=1:chr.g.len, size=chr.A.len, replace=FALSE)
        return(chr.genome[ind.smpld])
      })
      return( do.call("append", smpld) )
    }
    
  } else if(per.chromosome==FALSE & A.len<=g.len){
    ind.smpld <- sample(x=1:g.len, size=A.len, replace=FALSE)
    return(genome[ind.smpld])
  } else {
    stop("Length of A > length of genome.")
  }
  
}
################################################################################
