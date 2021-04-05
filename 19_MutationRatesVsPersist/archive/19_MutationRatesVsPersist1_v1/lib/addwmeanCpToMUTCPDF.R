############################################################################### 
# Specific function add weighted mean Cp values to MUTCP.DF.
###############################################################################
addwmeanCpToMUTCPDF <- function(MUTCP.DF, meanCp.dir, gcb){

  #load(file=paste0(mutCp.dir, "/", data.id, mut.id, "_", src.id, "_mutCalcPerCp.RData"))
  rownames(MUTCP.DF) <- paste0(MUTCP.DF$chr, ".", MUTCP.DF$bins)
  MUTCP.DF$wmeanCp0 <- NA
  MUTCP.DF$wmeanCp <- NA
  chr.v <- unique(MUTCP.DF$chr)
  
  # Add mean Cp of bins in MUTCP.DF
  for(chr in chr.v){
    
    chr.TF <- MUTCP.DF$chr==chr
    
    # BINWMEANCP.DF
    load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
    rownames(BINWMEANCP.DF) <- paste0(chr, ".", BINWMEANCP.DF$bin)
    
    if( !identical(as.numeric(MUTCP.DF[chr.TF,"bins"]), as.numeric(BINWMEANCP.DF$bin)) ){
      stop(paste0(chr, ": Checkpoint 1."))
    }
    
    if(!all(rownames(BINWMEANCP.DF)%in%rownames(MUTCP.DF))){
      stop(paste0(chr, ": Checkpoint 2."))
    }
    MUTCP.DF[rownames(BINWMEANCP.DF), "wmeanCp0"] <- BINWMEANCP.DF[rownames(BINWMEANCP.DF), "wmeanCp0"]
    MUTCP.DF[rownames(BINWMEANCP.DF), "wmeanCp"] <- BINWMEANCP.DF[rownames(BINWMEANCP.DF), "wmeanCp"]
    rm(BINWMEANCP.DF); gc()
    
  }
  
  return(MUTCP.DF)
  
}
