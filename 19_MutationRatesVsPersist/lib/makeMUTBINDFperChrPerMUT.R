################################################################################
# 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
### FUNCTION ###################################################################
makeMUTBINDFperChrPerMUT <- function(ncv.df='mutation data per chr per mutation type', 
                                     BINKMER.MX='base content of bins per chr',
                                     chr='chromosome', mut='mutation type',
                                     chr.len='chr length',
                                     bin.len='chr bin length'){
                                    
  
  if( nrow(ncv.df)==0 ){
    
    print(paste0("MUTBIN.DF not generated for ", chr, " ", mut, "..."), quote=F)
    return(NULL)
    
  }
  
  if( !identical(chr,unique(ncv.df$chr)) ){
    stop("MUTBIN.DF generation: ncv.df data not from input, ", chr, ".")
  }
    
  if(mut=="All"){
    
    #WT_SEQ <- "numUMChar"
    # For Nmsitenorm, consider only bins without missing sequence
    WT_SEQ <- "rSum"
    
  } else if( mut%in%c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") ){
    
    if( !identical(mut,unique(ncv.df$MUT)) ){
      stop("MUTBIN.DF generation: ncv.df data not from input, ", mut, ".")
    }
    WT_SEQ <- strsplit(x=mut, split=">", fixed=T)[[1]][1]
    
  } else {
    stop("MUTBIN.DF generation: ", mut, " invalid mutation type notation...")
  }
  
  #-------------------Mutation calculations per bin --> MUTBIN.DF
  
  # Assign mutations to bin
  ncv.df$bin <- ceiling(ncv.df$start/bin.len) 
  
  chr.bin <- ceiling(chr.len/bin.len)
  if( any(ncv.df[,"bin"]>chr.bin) ){
    stop(paste0("MUTBIN.DF generation ", chr, " ", mut, ": Bin exceeding chr."))
  }
  rm(chr.bin, chr.len)
  
  # Total mutation per bin across samples
  Tmut <- table(ncv.df[,"bin"])
  
  # Total mutation per bin divided by number of mutated sites
  x <- unique(paste0(ncv.df$start, "..", ncv.df$bin))
  x <- unlist(lapply(X=strsplit(x=x, split="..", fixed=T), FUN=function(x)x[2]))
  Nmsite <- table(x) # Number of mutated sites
  bin.srt <- as.character(sort(as.numeric(names(Tmut)), decreasing=F))
  if( !all(bin.srt%in%names(Nmsite)) ){
    stop(paste0("MUTBIN.DF generation ", chr, " ", mut, ": Checkpoint 1."))
  }
  TmutDIVNmsite <- Tmut[bin.srt]/Nmsite[bin.srt]
  
  # Number of mutated sites normalised to number of original/WT base per bin
  dimnames(BINKMER.MX)[[1]] <- BINKMER.MX[,"bins"]
  
  # Check if rows in BINKMER.MX don't have rows with both missing and defined values.
  a <- rowSums(x=BINKMER.MX[,c("A","C","G","T")], na.rm=F)
  b <- rowSums(x=BINKMER.MX[,c("A","C","G","T")], na.rm=T)
  if( sum(is.na(a) & b!=0)!=0 ){
    stop(paste0("MUTBIN.DF generation ", chr, " ", mut, ": Invalid BINKMER.MX."))
  }
  
  # Add number of complementary bases because I've already converted equivalent mutation types
  # (e.g. T>G == A>C so Nmsite should be normalised to A+T to get Nmsitenorm)
  AT.sum <- as.vector(apply(X=BINKMER.MX[,c("A", "T")], MARGIN=1, FUN=sum, na.rm=F))
  CG.sum <- as.vector(apply(X=BINKMER.MX[,c("C", "G")], MARGIN=1, FUN=sum, na.rm=F))
  BINKMER.MX[,"A"] <- AT.sum; BINKMER.MX[,"T"] <- AT.sum
  BINKMER.MX[,"C"] <- CG.sum; BINKMER.MX[,"G"] <- CG.sum
  
  if( !identical(as.numeric(a), as.numeric(BINKMER.MX[,"A"] + BINKMER.MX[,"C"])) ){
    stop(paste0("MUTBIN.DF generation ", chr, " ", mut, ": Problem with summing complementary bases."))
  }
  
  BINKMER.MX <- cbind(BINKMER.MX, rSum=a)
  rm(a, b, AT.sum, CG.sum)
  numWTSEQ <- BINKMER.MX[bin.srt,WT_SEQ]
  Nmsitenorm <- Nmsite[bin.srt]/numWTSEQ
  rm(BINKMER.MX); gc()
  
  MUTBIN.DF <- data.frame(chr=chr, bin=as.numeric(bin.srt), 
                          Tmut=as.numeric(Tmut[bin.srt]), 
                          Nmsite=as.numeric(Nmsite[bin.srt]),
                          TmutDIVNmsite=as.numeric(TmutDIVNmsite[bin.srt]),
                          Nmsitenorm=as.numeric(Nmsitenorm),
                          numWTSEQ=as.numeric(numWTSEQ),
                          stringsAsFactors=F)

  return(MUTBIN.DF)
  
}
################################################################################
makeMUTBINDFperChrPerMUT <- cmpfun(makeMUTBINDFperChrPerMUT, options=list(suppressUndefined=T))
################################################################################

# rm(list=ls()); gc()