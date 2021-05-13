################################################################################
# Add Cp data to MUTBIN.DF using BIN.MX to get MUTCP.DF, which is used for
# downstream plotting.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
# library(foreach)
# library(doParallel)
# library(itertools)
# source(paste0(lib, "/UTL_doPar.R"))
# source(paste0(wk.dir, "/lib/getRelevantBins.R"))
# source(paste0(wk.dir, "/lib/funxv.R"))
### FUNCTION ###################################################################
makeMUTBINDFtoContactMUTCPDFperMUT <- function(persist.dir = "directory of PERSIST.MX",
                                               basecont.dir = 'directory of BINKMER.MX', 
                                               basecont.affix = 'affix of BINKMER.MX indicating location',
                                               MUTBIN.DF = 'all chr', Cp.v,
                                               tot.bin.v = 'named vector of total bins per chr',
                                               gcb = 'gcb of PERSIST.MX and BINKMER.MX',
                                               nCPU.chr = 'parallel - chromosomes',
                                               nCPU.ij = 'parallel - contacts',
                                               funx.v = 'function for combining contact values',
                                               funx.id.v = 'same order as fun.v'){
  
  calc.id <- colnames(MUTBIN.DF)[ !colnames(MUTBIN.DF)%in%c("chr", "bin") ]
                                        
  chr.v <- unique(MUTBIN.DF$chr)
  chr.v.len <- length(chr.v)

  toExport <- c("chr.v", "MUTBIN.DF", "calc.id", "persist.dir", "gcb", "tot.bin.v",
                "basecont.dir", "basecont.affix")
  
  #### PARALLEL EXECUTION #########
  MUTCP <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                   .inorder=F, .combine="rbind",
                   .export=toExport, .noexport=ls()[!ls()%in%toExport]
                   
                   ) %op% {

    chunk <- sapply(X=itr, simplify=F, FUN=function(i){
      
      chr <- chr.v[i]
      chr.TF <- MUTBIN.DF$chr==chr
      val <- MUTBIN.DF[chr.TF, calc.id]
      names(val) <- MUTBIN.DF[chr.TF, "bin"]
      
      MUTCP.DF <- binToContact(v=val, binval.str.TF=F, missing.val=0,
                               funx.v=funx.v, funx.id.v=funx.id.v,
                               mxcalc.TF=T, persist.dir=persist.dir, gcb=gcb, 
                               chr=chr, ct="All", tot.bin=unname(tot.bin.v[chr])
                               )
      
      mutbin.v <- as.numeric(names(val))
      rm(val, chr.TF)
      
      # Take only contacts formed both by relevant bin
      relbin.v <- getRelevantbins(basecont.dir=basecont.dir, 
                                  basecont.affix=basecont.affix,
                                  gcb=gcb, chr=chr)
      relbin.v <- relbin.v[relbin.v[,"relevant"]==1,"bins"]
      MUTCP.DF <- MUTCP.DF[MUTCP.DF$i%in%relbin.v & MUTCP.DF$j%in%relbin.v, ]
      rm(relbin.v)
      
      MUTCP.DF$mutbin <- 0
  
      # 1 - Contacts formed by at least one mutated bin
      mut.TF <- MUTCP.DF$i%in%mutbin.v | MUTCP.DF$j%in%mutbin.v
      MUTCP.DF$mutbin[mut.TF] <- 1 
      rm(mut.TF)
      
      # Remove i and j columns
      #MUTCP.DF$i <- MUTCP.DF$j <- NULL
      #MUTCP.DF$chr <- chr # REMOVE
      rownames(MUTCP.DF) <- NULL
      
      return(MUTCP.DF)
      
    })
    
    return(do.call("rbind", chunk))
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  # Format for consistency with Bin-wise version
  colnames(MUTCP)[colnames(MUTCP)=="Cp"] <- "ind"
  MUTCP$ind <- factor(x=as.character(MUTCP$ind), levels=Cp.v)
  return(MUTCP)
  
}
################################################################################
makeMUTBINDFtoContactMUTCPDFperMUT <- cmpfun(makeMUTBINDFtoContactMUTCPDFperMUT,
                                             options=list(suppressUndefined=T))
################################################################################

# rm(list=ls()); gc()