################################################################################
# Draw grid estimates based on kernel estimate error
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(foreach)
#library(doParallel)
#library(itertools)
#source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
drawGridEst <- function(est.v=c(), err.v=c(), nboot=1000, SEED=4521, nCPU=1){
  
  set.seed(SEED)
  seed.v <- ceiling(runif(n=nboot, min=1, max=1000))
  
  toExport <- c("seed.v", "est.v", "err.v")
  #### PARALLEL EXECUTION #########
  BOOT <- foreach(itr=isplitVector(1:nboot, chunks=nCPU), 
                  .inorder=TRUE, .combine="cbind",
                  .export=toExport, .noexport=ls()[!ls()%in%toExport]
               
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      boot <- mapply(MEAN=est.v, SD=err.v, FUN=function(MEAN, SD){
        set.seed(seed.v[i])
        rnorm(n=1, mean=MEAN, sd=SD)
      })
      return(boot)
    })
    return( do.call("cbind", chunk) )
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  return(BOOT)
  
}
################################################################################
