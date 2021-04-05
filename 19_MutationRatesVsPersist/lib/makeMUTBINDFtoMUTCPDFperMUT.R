################################################################################
# Add Cp data to MUTBIN.DF using BIN.MX to get MUTCP.DF, which is used for
# downstream plotting.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
### FUNCTION ###################################################################
makeMUTBINDFtoMUTCPDFperMUT <- function(binmx.dir="directory of bin.mx",
                                        gcb="gcb of BIN.MX",
                                        Cp.v='1:21, Cp values in BIN.MX',
                                        bin.len='chr bin length',
                                        MUTBIN.DF="all chr"){
  
  Cp.v <- sort(Cp.v, decreasing=F)
  
  chr.v <- unique(MUTBIN.DF$chr)
  MUTCP.DF <- list()
  for(chr in chr.v){
    
    load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
    BIN.MX <- BIN.MX[,-(1:2)]
    
    # Initiate MUTCP.DF[[chr]]
    col.v <- colnames(MUTBIN.DF)[ !colnames(MUTBIN.DF)%in%c("chr", "bin") ]
    MUTCP.DF[[chr]] <- matrix(data=0, nrow=nrow(BIN.MX), ncol=1+length(col.v)+length(Cp.v),
                              dimnames=list(rownames(BIN.MX), c("mutbin", col.v, Cp.v))
    )
    
    if( is.unsorted(as.numeric(rownames(BIN.MX))) | any(duplicated(rownames(BIN.MX))) ){
      stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": Checkpoint 1."))
    }
    MUTCP.DF[[chr]] <- cbind.data.frame(chr=chr, bins=as.numeric(rownames(BIN.MX)), 
                                        MUTCP.DF[[chr]], stringsAsFactors=F)
    
    # Assign bin to Cp (involved in at least 1 contact for that Cp in any tissue)
    for(Cp in Cp.v){
      
      Cp.ind <- grep(x=colnames(BIN.MX), pattern=paste0("s_Cp_", Cp, "_ct"))
      if( length(Cp.ind)!=length(Cp.v) ){
        stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": Checkpoint 2."))
      }
      # Rowsums only for at least 2 rows
      MUTCP.DF[[chr]][,as.character(Cp)] <- as.numeric(rowSums(x=BIN.MX[,Cp.ind])>0)
      rm(Cp.ind)
      
    }
    
    chr.TF <- MUTBIN.DF$chr==chr
    bin.mut <- as.character(MUTBIN.DF$bin[chr.TF])
    MUTCP.DF[[chr]][bin.mut,col.v] <- MUTBIN.DF[chr.TF,col.v]
    # Mark bins with >=1 mutation
    MUTCP.DF[[chr]][bin.mut,"mutbin"] <- 1
    
    rm(col.v, bin.mut, chr.TF); gc()
    
  } # chr.v for loop end
  
  MUTCP.DF <- do.call("rbind.data.frame", c(MUTCP.DF, stringsAsFactors=F))
  rownames(MUTCP.DF) <- NULL
  
  return(MUTCP.DF)
  
}
################################################################################
makeMUTBINDFtoMUTCPDFperMUT <- cmpfun(makeMUTBINDFtoMUTCPDFperMUT,
                                      options=list(suppressUndefined=T))
################################################################################

# rm(list=ls()); gc()