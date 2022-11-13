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
                                  coord.system = 1,
                                  p.coord.mx, # 1-based CS, [start.a, end.a, start.b, end.b]
                                  consensus.FUN){
  
  dimnames(p.coord.mx)[[2]] <- c("start.a", "end.a", "start.b", "end.b")
  
  #
  if(coord.system == 1){
    print("getConsensusCpOfPairs(): p.coord.mx is 1-based.")
  }  else if(coord.system == 0){
    p.coord.mx[,"start.a"] <- p.coord.mx[,"start.a"] +  1
    p.coord.mx[,"start.b"] <- p.coord.mx[,"start.b"] +  1
    print("getConsensusCpOfPairs(): p.coord.mx converted to 1-based.")
  } else {
    rm(p.coord.mx)
    stop("getConsensusCpOfPairs(): coord.system can only be either 0 or 1 (default).")
  }
    
  # Check for missing coordinates
  is.noNA <- !is.na(rowSums(p.coord.mx, na.rm=F))
  if( any(!is.noNA) ){
    warning("getConsensusCpOfPairs(): Missing coordinates in p.coord.mx.
             NA will be returned for these rows/cases.")
  }
  
  if( all(!is.noNA) ){
    rm(p.coord.mx)
    stop("getConsensusCpOfPairs(): All pairs have missing coordinates.")
  }
  
  # Check for negative lengths
  if( any( (p.coord.mx[is.noNA,2] < p.coord.mx[is.noNA,1]) |
           (p.coord.mx[is.noNA,4] < p.coord.mx[is.noNA,3])
         ) ){
    
    rm(p.coord.mx)
    stop("getConsensusCpOfPairs(): Negative feature lengths.")
  
  }
  
  # Convert feature coordinates to bins
  p.coord.mx <- ceiling(p.coord.mx / bin.len)
  
  min.gap <- min(ij.mx[,2] - ij.mx[,1])
  pair.len <- length(p.coord.mx[,1])
  noNApairs.ind <- (1:pair.len)[is.noNA]
  funx.len <- length(consensus.FUN)
  
  toExport <- c("p.coord.mx", "funx.len", "ij.mx", "min.gap", "consensus.FUN")
  
  #### PARALLEL EXECUTION #########
  
  noNA.PAIRCP.MX <- foreach(pair.ind.v=isplitVector(noNApairs.ind, chunks=nCPU), 
                            .inorder=F, .combine="rbind",
                            .export=toExport, .noexport=ls()[!ls()%in%toExport]
                       
  ) %op% {
    
    chunk <- sapply(X=pair.ind.v, simplify=F, FUN=function(pair.ind){
      
      print(pair.ind, quote=F)
      
      a.bins <- p.coord.mx[pair.ind,1]:p.coord.mx[pair.ind,2]
      b.bins <- p.coord.mx[pair.ind,3]:p.coord.mx[pair.ind,4]
      
      max.gap.p <- diff( range(c(a.bins, b.bins)) )
      
      consCp <- rep(NA, times=funx.len)
      if(max.gap.p >= min.gap){
        
        test.ab <- ij.mx[,"i"]%in%a.bins & ij.mx[,"j"]%in%b.bins
        test.ba <- ij.mx[,"i"]%in%b.bins & ij.mx[,"j"]%in%a.bins
        test.ij <- test.ab | test.ba
        vals <- ij.mx[test.ij,"Cp"]
        
        if( sum(test.ij) > 0 ){
          consCp <- sapply(X=consensus.FUN, simplify=T, USE.NAMES=F, FUN=function(fnx){
            do.call(what=fnx, list(vals))
          })
        }
        
        rm(test.ab, test.ba, vals)
        
      }
      
      return(c(pair.ind, consCp))
      
    })
    
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  value.cols <- paste0(consensus.FUN, ".consCp")
  dimnames(noNA.PAIRCP.MX)[[2]] <- c("pair.ind", value.cols)
  
  PAIRCP.MX <- matrix( data=NA, nrow=pair.len, ncol=length(noNA.PAIRCP.MX[1,]),
                       dimnames=dimnames(noNA.PAIRCP.MX) ) 
  PAIRCP.MX[,"pair.ind"] <- as.integer(1:pair.len)

  PAIRCP.MX[noNA.PAIRCP.MX[,"pair.ind"],] <- noNA.PAIRCP.MX
  
  return(PAIRCP.MX)
  
}

################################################################################

# rm(list=ls()); gc()