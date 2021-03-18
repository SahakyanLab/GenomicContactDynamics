################################################################################
# Function to combine masks and reduce intervals. By default strand is "*"; 
# modify if needed. 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(data.table)
# library(GenomicRanges)
################################################################################
getMergedReducedMask <- function(maskbed.dir, #"/dir"
                                 mask, #"ct_ESC_foi_ATF2_desc_TF;ct_ESC_foi_ASH2L_desc_TF",
                                 splitChar #";"
                                 ){ 
                                   
  mask <- strsplit(x=mask, split=";")[[1]]
  
  # Merge combination of masks
  mergedmask <- sapply(X=mask, simplify=FALSE, FUN=function(x){
    x <- list.files(path=maskbed.dir, pattern=x, full.names=TRUE, recursive=FALSE,
                    include.dirs=FALSE)
    x <- fread(file=x, colClasses=list(character=1, integer=2:3),
               select=1:3, col.names=c("chr", "start", "end"),
               stringsAsFactors=FALSE, data.table=FALSE, header=FALSE)
  })
  mergedmask <- do.call("rbind.data.frame", c(mergedmask, stringsAsFactors=FALSE))
  rownames(mergedmask) <- NULL
  
  # Reduce overlapping ranges per chr
  mergedmask <- GenomicRanges::GRanges(seqnames=mergedmask$chr,
                                       IRanges(start=mergedmask$start,
                                               end=mergedmask$end),
                                       strand=rep(x="*", times=length(mergedmask[,1]))
                                       )
  mergedmask <- GenomicRanges::reduce(mergedmask)
  mergedmask <- as.data.frame(mergedmask)
  mergedmask <- mergedmask[,c("seqnames", "start", "end", "strand")]
  data.table::setnames(x=mergedmask, old="seqnames", new="chr")
  
  mergedmask$chr <- as.character(mergedmask$chr)
  mergedmask$strand <- as.character(mergedmask$strand)
  
  return(mergedmask)
  
}
