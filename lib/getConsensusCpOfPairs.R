################################################################################
# Function to get consensus Cp (e.g., max Cp, mean Cp) of a pair of features 
# based on contacts linking them.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(foreach)
#library(doParallel)
#library(itertools)
#source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
getConsensusCpOfPairs <- function(bin.len, nCPU,
                                  ij.mx, # [i, j, Cp] 
                                  p.coord.mx, # 1-based CS, [start.a, end.a, start.b, end.b]
                                  consensus.FUN){
  
  
  min.gap <- min(ij.mx[,2] - ij.mx[,1])
  
  dimnames(p.coord.mx)[[2]] <- c("start.a", "end.a", "start.b", "end.b")
  
  # Check for negative lengths
  if( any( (p.coord.mx[,2] < p.coord.mx[,1]) |
           (p.coord.mx[,4] < p.coord.mx[,3])
         ) ){
    
    rm(p.coord.mx)
    stop("Negative feature lengths.")
  
  }
  
  # Convert feature coordinates to bins
  p.coord.mx <- ceiling(p.coord.mx/bin.len)
  
  toExport <- c("p.coord.mx", "ij.mx")
  pair.len <- length(p.coord.mx[,1])
  
  #### PARALLEL EXECUTION #########
  
  PAIRCP.MX <- foreach(pair.ind.v=isplitVector(1:pair.len, chunks=nCPU), 
                       .inorder=F, .combine="rbind",
                       .export=toExport, .noexport=ls()[!ls()%in%toExport]
                       
  ) %op% {
    
    chunk <- sapply(X=pair.ind.v, simplify=F, FUN=function(pair.ind){
      
      print(pair.ind, quote=F)
      
      a.bins <- p.coord.mx[pair.ind,1]:p.coord.mx[pair.ind,2]
      b.bins <- p.coord.mx[pair.ind,3]:p.coord.mx[pair.ind,4]
      
      max.gap.p <- diff( range(c(a.bins, b.bins)) )
      
      if(max.gap.p < min.gap){
        
        consCp <- NA
        
      } else {
        
        test.ab <- ij.mx[,"i"]%in%a.bins & ij.mx[,"j"]%in%b.bins
        test.ba <- ij.mx[,"i"]%in%b.bins & ij.mx[,"j"]%in%a.bins
        testij <- test.ab | test.ba
        rm(test.ab, test.ba)
        
        consCp <- ifelse( sum(testij)==0,  NA,
                          sapply(X=ij.mx[testij,"Cp"], simplify=T, FUN=consensus.FUN))
      }
      
      return( c(pair.ind, consCp) )
      
    })
    
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
}

################################################################################

# rm(list=ls()); gc()