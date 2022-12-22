################################################################################
# Function to get consensus Cp (e.g., max Cp, mean Cp) of a pair of features 
# based on contacts linking them. 
# The new version offers two ways to do calculation 1) slow per-contact method 
# and 2) faster data.table method but currently assumes that p.coord.mx 
# follows ij.mx contact where i < j i.e. it consists of non-overlapping region 
# pairs and that the region defined by the first two columns should be upstream 
# of the other region. They also differ on the format of consensus.FUN argument
# they can accept. 1) can accept a string vector of function names while for 2)
# consensus.FUN should be a function object defined by function(). Method one is
# stable but shaky and it should not be used inside a foreach parallel 
# implementation as it needs foreach parallel computing two. 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
## Per-contact method
#library(foreach)
#library(doParallel)
#library(itertools)
#source(paste0(lib, "/UTL_doPar.R"))

## data.table method
#library(data.table)
#library(GenomicRanges)
#source(paste0(lib, "/TranToRextr/GEN_WhichOverlap.R"))
### FUNCTION ###################################################################
getConsensusCpOfPairs <- function(bin.len,
                                  ij.mx, # [i, j, Cp] 
                                  coord.system = 1,
                                  p.coord.mx, # 1-based CS, [start.a, end.a, start.b, end.b]
                                  consensus.FUN, # Depends on method, see description above
                                  methd, # INT, 1 or 2
                                  nCPU=1
                                  ){
  
  dimnames(p.coord.mx)[[2]] <- c("start.a", "end.a", "start.b", "end.b")
  
  # Initialise output
  
  pair.len <- length(p.coord.mx[,1])
  PAIRCP.MX <- data.frame(pair.ind=1:pair.len, value=rep(NA))
  
  # Check for missing coordinates
  is.noNA <- !is.na(rowSums(p.coord.mx, na.rm=F))
  noNApairs.ind <- (1:pair.len)[is.noNA]
  
  if( all(!is.noNA) ){
    
    warning("getConsensusCpOfPairs(): All pairs have missing coordinates. 
            Returning output with all NAs.")
    
  } else {
    
    if( any(!is.noNA) ){
      warning("getConsensusCpOfPairs(): Missing coordinates in p.coord.mx.
              NA will be returned for these rows/cases.")
    }
    
    # Check coordinate system
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
    
    # Check for negative lengths
    if( any( (p.coord.mx[is.noNA,2] < p.coord.mx[is.noNA,1]) |
             (p.coord.mx[is.noNA,4] < p.coord.mx[is.noNA,3])
    ) ){
      rm(p.coord.mx)
      stop("getConsensusCpOfPairs(): Negative feature lengths.")
    }
    
    #
    
    if(methd == 1){
      
      # Convert feature coordinates to bins
      p.coord.mx <- ceiling(p.coord.mx / bin.len)
      min.gap <- min(ij.mx[,2] - ij.mx[,1])
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
            vals <- unname( ij.mx[test.ij,"Cp"] )
            
            if( sum(test.ij) > 0 ){
              consCp <- sapply(X=consensus.FUN, simplify=T, USE.NAMES=F, FUN=function(fnx){
                
                out <- do.call(what=fnx, list(vals))
                if( is.atomic(out) & length(out) == 1 ){
                  return(out)
                } else if( is.atomic(out) & length(out) > 1 ){
                  out <- paste(out, collapse=";")
                  return(out)
                } else {
                  rm(out)
                  stop("getConsensusCpOfPairs(): Checkpoint in consensus Cp return.")
                }
                
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
      
      # Initialise output object
      PAIRCP.MX <- matrix( data=NA, nrow=pair.len, ncol=length(noNA.PAIRCP.MX[1,]),
                           dimnames=dimnames(noNA.PAIRCP.MX) ) 
      PAIRCP.MX <- as.data.frame(PAIRCP.MX)
      PAIRCP.MX[,"pair.ind"] <- as.integer(1:pair.len)
      
      # 
      noNA.PAIRCP.MX <- as.data.frame(noNA.PAIRCP.MX)
      noNA.PAIRCP.MX[,"pair.ind"] <- as.integer(noNA.PAIRCP.MX[,"pair.ind"])
      
      dim.obj <- list(dim(PAIRCP.MX[noNA.PAIRCP.MX[,"pair.ind"],]),
                      dim(noNA.PAIRCP.MX))
      if( !identical(dim.obj[[1]], dim.obj[[2]]) ){
        rm(noNA.PAIRCP.MX)
        stop("getConsensusCpOfPairs(): Checkpoint dimension.")
      } else {
        PAIRCP.MX[noNA.PAIRCP.MX[,"pair.ind"],] <- noNA.PAIRCP.MX
      }
      
    } else if(methd == 2){
      
      if( !is.function(consensus.FUN) ){
        rm(ij.mx)
        stop("getConsensusCpOfPairs(): consensus.FUN should be a function for method 2.")
      }
      
      #
      
      ## GRanges way
      
      # end.i <- ij.mx[,"i"] * bin.len
      # olap.i <- WhichOverlap(start.query=p.coord.mx[is.noNA,1],
      #                        end.query=p.coord.mx[is.noNA,2],
      #                        space.query="a",
      #                        start.subject=end.i - bin.len + 1,
      #                        end.subject=end.i, space.subject="a",
      #                        maxgap=-1L, minoverlap=1L, type="any")
      # 
      # end.j <- ij.mx[,"j"] * bin.len
      # olap.j <- WhichOverlap(start.query=p.coord.mx[is.noNA,3],
      #                        end.query=p.coord.mx[is.noNA,4],
      #                        space.query="b",
      #                        start.subject=end.j - bin.len + 1,
      #                        end.subject=end.j, space.subject="b",
      #                        maxgap=-1L, minoverlap=1L, type="any")
      # rm(end.i, end.j)

      ## data.table way
      
      # query / x
      i.pc.DT <- data.table(start=p.coord.mx[is.noNA,"start.a"], end=p.coord.mx[is.noNA,"end.a"])
      j.pc.DT <- data.table(start=p.coord.mx[is.noNA,"start.b"], end=p.coord.mx[is.noNA,"end.b"])
      
      # subject / y
      end.i <- ij.mx[,"i"] * bin.len
      end.j <- ij.mx[,"j"] * bin.len
      ij.len <- length(ij.mx[,1])
      i.ij.DT <- data.table(start=end.i - bin.len + 1, end=end.i, orig.ind=1:ij.len)
      j.ij.DT <- data.table(start=end.j - bin.len + 1, end=end.j, orig.ind=1:ij.len)
      
      setkey(i.ij.DT, start, end)
      olap.i <- foverlaps(x=i.pc.DT, y=i.ij.DT, type="any", which=TRUE, nomatch=NULL)
      
      setkey(j.ij.DT, start, end)
      olap.j <- foverlaps(x=j.pc.DT, y=j.ij.DT, type="any", which=TRUE, nomatch=NULL)
      
      # Convert subject / y keyed index (sorted) to original index (orig.ind)
      olap.i[["yid"]] <- i.ij.DT[["orig.ind"]][ olap.i[["yid"]] ]
      olap.j[["yid"]] <- j.ij.DT[["orig.ind"]][ olap.j[["yid"]] ]
      
      rm(i.pc.DT, j.pc.DT, i.ij.DT, j.ij.DT, end.i, end.j)
      
      # Get contacts with both i and j regions overlapping between two region sets
      
      olap.ij <- setDT(rbind.data.frame(olap.i, olap.j))
      with.olap.ij <- olap.ij[duplicated(olap.ij),]
      rm(olap.ij, olap.i, olap.j)
      
      #
      
      colnames(with.olap.ij) <- c("query", "subject")
      with.olap.ij[["Cp"]] <- ij.mx[with.olap.ij[["subject"]], "Cp"]
      with.olap.ij[["subject"]] <- NULL
      with.olap.ij[["query"]] <- noNApairs.ind[ with.olap.ij[["query"]] ]
      setnames(x=with.olap.ij, old="query", new="pair.ind")
      
      # Apply consensus.FUN on overlapping Cps
      
      noNA.PAIRCP.MX <- with.olap.ij[, consensus.FUN(Cp), by=pair.ind]
      colnames(noNA.PAIRCP.MX) <- colnames(PAIRCP.MX)
      PAIRCP.MX[ noNA.PAIRCP.MX[["pair.ind"]], ] <- noNA.PAIRCP.MX
      
    } else {
      stop("getConsensusCpOfPairs(): Invalid methd argument. Should either be 
            integer 1 or 2 currently.")
    }
    
  }
  
  return(PAIRCP.MX)
  
}

################################################################################

# rm(list=ls()); gc()
