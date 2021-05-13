################################################################################
# Add Cp data to MUTBIN.DF using BIN.MX to get MUTCP.DF, which is used for
# downstream plotting.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
# library(reshape2)
### FUNCTION ###################################################################
makeMUTBINDFtoBinMUTCPDFperMUT <- function(binmx.dir="directory of bin.mx",
                                        gcb="gcb of BIN.MX",
                                        Cp.v='1:21, Cp values in BIN.MX',
                                        bin.len='chr bin length',
                                        MUTBIN.DF="all chr",
                                        basecont.dir, basecont.affix){
  
  calc.id <- colnames(MUTBIN.DF)[ !colnames(MUTBIN.DF)%in%c("chr", "bin") ]
  
  # Check MUTBIN.DF; should only contain bins with at least 1 mutation.
  TF.v <- c( calc.id%in%c("Tmut", "Nmsite", "TmutDIVNmsite") & 
             any( MUTBIN.DF[[calc.id]]<=0 | !is.finite(MUTBIN.DF[[calc.id]]) ),
             calc.id=="numWTSEQ" & 
             any( MUTBIN.DF[[calc.id]]<1 & !is.na(MUTBIN.DF[[calc.id]]) ),
             calc.id%in%c("Nmsitenorm", "Tmutnorm") & 
             !identical(is.na(MUTBIN.DF[[calc.id]]), !is.finite(MUTBIN.DF[[calc.id]]))
            )
  if( any(TF.v) ){
    stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": Invalid MUTBIN.DF"))
  }
  rm(TF.v)
  
  #-------------------Get MUTCP.DF
  Cp.v <- sort(Cp.v, decreasing=F)
  
  chr.v <- unique(MUTBIN.DF$chr)
  MUTCP.DF <- list()
  for(chr in chr.v){
    
    load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
    BIN.MX <- BIN.MX[,-(1:2)]
    
    # Initiate MUTCP.DF[[chr]]
    MUTCP.DF[[chr]] <- matrix(data=0, nrow=nrow(BIN.MX), ncol=1+length(calc.id)+length(Cp.v),
                              dimnames=list(rownames(BIN.MX), c("mutbin", calc.id, Cp.v))
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
    MUTCP.DF[[chr]][bin.mut,calc.id] <- MUTBIN.DF[chr.TF,calc.id]
    
    # Mark bins with >=1 mutation
    MUTCP.DF[[chr]][bin.mut,"mutbin"] <- 1
    
    #-------------------Mark NA those bins not overlapping with given location
    load(file=paste0(basecont.dir, "/", chr, "_BinKmer1", basecont.affix, ".RData"))
    if( !identical(rownames(MUTCP.DF[[chr]]), as.character(BINKMER.MX[,"bins"])) ){
      stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": MUTCP.DF and BINKMER.MX don't match."))
    }
    if( length(BIN.MX[,1])!=length(BINKMER.MX[,1]) ){
      stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": BIN.MX and BINKMER.MX different row length."))
    }
    rm(BIN.MX)
    
    # Update numWTSEQ for all bins of chr, MUTBIN.DF only contains info for mutated bins
    a <- rowSums(x=BINKMER.MX[,c("A","C","G","T")], na.rm=F)
    if("numWTSEQ"%in%calc.id){
      MUTCP.DF[[chr]][,"numWTSEQ"] <- a
    }
    
    # Mark NA for mutbin, those bins not overlapping with given location and ...
    drop.TF <- (a==0 & !is.na(a))
    if( any(drop.TF & MUTCP.DF[[chr]][,"mutbin"]==1) ){
      stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": Checkpoint 3."))
    }
    # ...those with missing sequence
    drop.TF <- drop.TF | is.na(a)
    MUTCP.DF[[chr]][drop.TF,"mutbin"] <- NA
    
    if( !identical( is.na(MUTCP.DF[[chr]]$mutbin), (is.na(a) | a==0) ) ){
      stop(paste0("MUTCP.DF generation ", gcb, " ", chr, ": Error in marking invalid bins."))
    }
    
    # By default, exclude last bin in analyses because of <Hi-C bin length
    MUTCP.DF[[chr]][length(a),"mutbin"] <- NA
    
    rm(a, drop.TF, bin.mut, chr.TF, BINKMER.MX, chr); gc()
    
  } # chr.v for loop end
  
  MUTCP.DF <- do.call("rbind.data.frame", c(MUTCP.DF, stringsAsFactors=F))
  rownames(MUTCP.DF) <- NULL
  
  MUTCP.DF <- reshape2::melt(data=MUTCP.DF, id.vars=c("chr", "bins", "mutbin", calc.id),
                             stringsAsFactors=F)
  colnames(MUTCP.DF)[colnames(MUTCP.DF)=="variable"] <- "ind"
  colnames(MUTCP.DF)[colnames(MUTCP.DF)=="value"] <- "hasInd"
  
  return(MUTCP.DF)
  
}
################################################################################
makeMUTBINDFtoBinMUTCPDFperMUT <- cmpfun(makeMUTBINDFtoBinMUTCPDFperMUT,
                                         options=list(suppressUndefined=T))
################################################################################

# rm(list=ls()); gc()