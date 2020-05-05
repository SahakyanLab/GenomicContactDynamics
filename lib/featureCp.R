################################################################################
# Determine cp of feature
# Uses FEATURE.BIN.MX, the output of mapToHiCcontactBins() and PERSIST.MX
# Dependencies:
## library(foreach)
## library(itertools)
## source(paste0(lib, "/UTL_doPar.R"))
################################################################################
featureCp <- function(
  # FEATURE.BIN.MX (should only contain data from one chromosome)
  FEATUREBINMXobj = feat.bin.mx,
  feature = "name", 
  # PERSIST.MX
  PERSISTMXdir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc",
  OUTPUTdir = paste0(objective.dir, "/out_topoCp"),
  chr = "chr21",
  gcb = "min2Mb",
  out.name = "MCF7_TOP2B_within",
  # Max is number of contacts per chromosome
  # No parallel execution if nCPU=1L
  nCPU=nCPU,
  # Return COUNT.MX (final output)?
  returnObj=TRUE
){
  ################################################################################
  # LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
  ################################################################################
  library(foreach)
  library(itertools)
  
  if(nCPU > 1){
    suppressWarnings(library("doParallel"))
    registerDoParallel(cores = nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  ################################################################################
  # MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
  ################################################################################
  affix <- paste0(chr, "_", gcb)
  
  # Remove strings in feature column (containing the uniqueID of each feature)
  # of FEATURE.BIN.MX
  ind <- grep(x=FEATUREBINMXobj$name, pattern=";", fixed=TRUE)
  if(length(ind)!=0){
    lst <- strsplit(x=FEATUREBINMXobj[ind, feature], split=";")
    names(lst) <- ind
    str.sep <- foreach(i=ind, .combine="rbind.data.frame")%do%{
      
      data.frame( chr=rep(FEATUREBINMXobj[i,"chr"]) , 
                  bin=as.numeric( rep(FEATUREBINMXobj[i,"bin"]) ),
                  name=as.numeric( lst[[as.character(i)]] ), stringsAsFactors=FALSE )
    }
    FEATUREBINMXobj <- rbind(FEATUREBINMXobj[-(ind), c("chr", "bin", feature)], str.sep)
    rm(lst, ind, str.sep); gc()
  }
 
  # Load PERSIST.MX
  load(paste0(PERSISTMXdir, "/", chr, "_Persist_", gcb, ".RData"))
  ij.mx <- cbind(PERSIST.MX$hits[,c("i", "j")], cp=PERSIST.MX$ntis)
  rm(PERSIST.MX); gc()
  
  # Feature-centric
  feat.v <- unique(FEATUREBINMXobj[[feature]])
  feat.v.len <- length(feat.v)
  
  #### PARALLEL EXECUTION #########

  FEATCP.DF <- foreach(itr=isplitVector(1:feat.v.len, chunks=nCPU),
                       .combine="rbind", .inorder=TRUE,
                       .export=c("ij.mx", "FEATUREBINMXobj", 
                                 "feat.v", "chr"),
                       .noexport=ls()[!ls()%in%c("ij.mx", "FEATUREBINMXobj", 
                                                 "feat.v", "chr")]
  )%op%{
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      feat <- feat.v[i]
      feat.bins <- FEATUREBINMXobj[ FEATUREBINMXobj[[feature]]==feat, "bin" ]
      cp.v <- ij.mx[ ij.mx[,"i"]%in%feat.bins | ij.mx[,"j"]%in%feat.bins, "cp"]
      
      countPerCp <- rep(0L, times=21)
      
      for(cp in 1:21){
        count <- sum(cp.v==cp)
        countPerCp[cp] <- count
      }
         
      rm(feat.bins, cp.v); gc()
      return( 
        data.frame(matrix( countPerCp, ncol=21,
                           dimnames=list(NULL, as.character(1:21)) ) ) 
              )
      
    }) 
    
    return(do.call("rbind", chunk)) 
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  
  FEATCP.DF <- cbind(chr=rep(chr), feature=feat.v, FEATCP.DF)
  rownames(FEATCP.DF) <- NULL
  save(FEATCP.DF, file=paste0( OUTPUTdir, "/", affix, "_", 
                               out.name, "_featureCp.RData") )
  
  if(returnObj==TRUE){ return(FEATCP.DF) }
  
}
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
featureCp <- cmpfun(featureCp, options=list(suppressUndefined=TRUE))
################################################################################