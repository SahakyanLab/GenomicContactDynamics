################################################################################
# Function to combine masks and reduce intervals
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(data.table)
# library(IRanges)
################################################################################
getMergedReducedMask <- function(maskbed.dir="/dir",
                                 mask="ct_ESC_foi_ATF2_desc_TF;ct_ESC_foi_ASH2L_desc_TF",
                                 splitChar=";"){
  
  mask <- strsplit(x=mask, split=";")[[1]]
  
  # Merge combination of masks
  mergedmask <- sapply(X=mask, simplify=FALSE, FUN=function(x){
    x <- list.files(path=maskbed.dir, pattern=x, full.names=TRUE, recursive=FALSE,
                    include.dirs=FALSE)
    x <- fread(file=x,colClasses=list(character=1, integer=2:3),
               select=1:3, col.names=c("chr", "start", "end"),
               stringsAsFactors=FALSE, data.table=FALSE, header=FALSE)
  })
  mergedmask <- do.call("rbind.data.frame", mergedmask)
  rownames(mergedmask) <- NULL
  
  chr.v <- unique(mergedmask$chr)
  
  # Reduce overlapping ranges per chr
  mergedmask <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
    
    x <- mergedmask[mergedmask$chr==chr,]
    x <- data.frame(reduce(
      IRanges( start=as.numeric(x[,"start"]),
               end=as.numeric(x[,"end"]) )
    ), stringsAsFactors=FALSE)
    # Remove width column
    return( cbind.data.frame(chr=chr, x[,-3], stringsAsFactors=FALSE) )
    
  })
  
  return(mergedmask)
  
}
