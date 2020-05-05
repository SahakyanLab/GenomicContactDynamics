################################################################################
# Count number of features per contact
# Uses FEATURE.BIN.MX, the output of mapToHiCcontactBins() and PERSIST.MX
# Dependencies:
## library(foreach)
## library(itertools)
## source(paste0(lib, "/UTL_doPar.R"))
################################################################################
CountPerHiCcontact <- function(
  # FEATURE.BIN.MX (should only contain data from one chromosome)
  FEATUREBINMXobj = feat.bin.mx,
  # PERSIST.MX
  PERSISTMXdir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc",
  OUTPUTdir = output.dir,
  chr = "chr21",
  gcb = "min2Mb",
  out.name = "MCF7_TOP2B_within",
  # Max is number of contacts per chromosome
  nCPU=4L,
  # Return COUNT.MX (final output)?
  returnObj=TRUE
){
  
  ################################################################################
  # LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
  ################################################################################
  library(foreach)
  library(itertools)
  source(paste0(lib, "/UTL_doPar.R"))
  ################################################################################
  # MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
  ################################################################################
  affix <- paste0(chr, "_", gcb)
  
  # Load PERSIST.MX
  load(paste0(PERSISTMXdir, "/", chr, "_Persist_", gcb, ".RData"))
  ij.mx <- cbind(PERSIST.MX$hits[,c("i", "j")], ntis=PERSIST.MX$ntis)
  hits.len <- nrow(ij.mx)
  rm(PERSIST.MX)
  
  #### PARALLEL EXECUTION #########
  
  COUNT.MX <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                      .combine="c", .inorder=FALSE,
                      .export=c("ij.mx", "FEATUREBINMXobj"),
                      .noexport=ls()[!ls()%in%c("ij.mx", "FEATUREBINMXobj")]
                      
  )%op%{
    
    chunk <- sapply(X=itr, FUN=function(i){
      
      test <- FEATUREBINMXobj$bin%in%unlist(ij.mx[i,c("i","j")])
      sum( FEATUREBINMXobj[test, "countPerBin"] )
      
    })
    return(chunk)
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  COUNT.MX <- cbind(cp=ij.mx[,"ntis"], count=COUNT.MX)
  save(COUNT.MX, file=paste0(OUTPUTdir, "/", affix, "_", 
                             out.name, "_countPerContact.RData"))
  
  if(returnObj==TRUE){ return(COUNT.MX) }

}
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
CountPerHiCcontact <- cmpfun(CountPerHiCcontact, options=list(suppressUndefined=TRUE))
################################################################################


