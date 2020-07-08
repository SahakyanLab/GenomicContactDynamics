################################################################################
# Check overlap of features (provided as bed file) at certain genomic positions
# relative to supplied midpoints
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
#library(compiler)
#library(foreach)
#library(itertools)
#library(IRanges)
#source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
################################################################################
olapPerGenomicPos <- function(
  foi.id = foi,
  foi.df = "df bed file format",
  start = 1,
  end = 40000,
  res = 40000,
  pos.v = -2:2,
  nCPU = 3,
  minPos = 1,
  maxPos = "end of chromosome",
  # Overlap parameters
  min.olap = 1L,
  max.gap = -1L, # for Mac, 0L for Linux
  type.olap = "any"
){
  
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    #print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  len <- length(start)
  toExport <- c("foi.df", "start", "end", "pos.v", "minPos", "maxPos",
                "min.olap", "max.gap", "type.olap")
  #### PARALLEL EXECUTION #########
  MX <- foreach(itr=isplitVector(1:len, chunks=nCPU), 
                .inorder=FALSE, .combine="rbind",
                .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      start.actual <- start[i]+(pos.v*res)
      end.actual <- end[i]+(pos.v*res)
      incl.TF <- !(start.actual<minPos | start.actual>maxPos | end.actual<minPos | end.actual>maxPos)
      toOlap <- sum(incl.TF)
      
      if(toOlap!=0L){
        # Query   <- positions
        # Subject <- features 
        # Note that WhichOverlap() only gives the indices of positions with overlap
        olap <- WhichOverlap(start.query=start.actual[incl.TF], 
                             end.query=end.actual[incl.TF], 
                             space.query=rep("a", toOlap),
                             start.subject=foi.df[,2], 
                             end.subject=foi.df[,3], 
                             space.subject=rep("a",length(foi.df[,2])),
                             maxgap=max.gap, minoverlap=min.olap,
                             type=type.olap)
      }
      foicount <- rep(NA, times=length(pos.v))
      foicount[incl.TF] <- 0L
      
      if( nrow(olap)!=0L ){
        olap <- table(olap[,"query"])
        # Positions with overlap
        foicount[incl.TF][ as.numeric(names(olap)) ] <- olap
      } 
      return(foicount)
    })
    return(do.call("rbind", chunk))
  }
  ### END OF PARALLEL EXECUTION ###
  dimnames(MX)[[2]] <- as.character(format(x=pos.v, scipen=FALSE))
  return(MX)

}
################################################################################
olapPerGenomicPos <- cmpfun(olapPerGenomicPos, options=list(suppressUndefined=TRUE))
################################################################################
