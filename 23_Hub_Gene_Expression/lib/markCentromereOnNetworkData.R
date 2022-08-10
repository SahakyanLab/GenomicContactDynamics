################################################################################
# Modify network data object to mark centromeres.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# centrobed.file = paste0(wk.dir, "/txTable/ct_hg19_foi_centromoreonly_desc_DNA")
# library(GenomicRanges)
# source(paste0(wk.dir, "/lib/GEN_WhichOverlap.R"))
### FUNCTION ###################################################################
markCentromereOnNetworkData <- function(NETWRK = 'list(edges, nodes)',
                                        chr = 'chromosome',
                                        bin.len = 'resolution of regions represented
                                                   by nodes',
                                        centrobed.file = 'bed file of centromere
                                                          position'
){
  
  cent.v <- read.delim(file=centrobed.file, sep="\t", header=F, stringsAsFactors=F)[,1:3]
  cent.v <- as.numeric(cent.v[cent.v[,1]==chr,2:3])
  cent.v <- ceiling(cent.v/bin.len)
  
  # Mark nodes overlapping with centromere
  node.TF <- as.numeric(NETWRK$nodes$id)%in%(cent.v[1]:cent.v[2]) 
  NETWRK$nodes$shape[node.TF] <- "square"
  rm(node.TF)
  
  # Mark edges overlapping with centromere
  edge.len <- length(NETWRK$edges[,1])
  hits.mx <- WhichOverlap(start.subject=cent.v[1], 
                          end.subject=cent.v[2], 
                          space.subject=chr,
                          start.query=NETWRK$edges$from,
                          end.query=NETWRK$edges$to,
                          space.query=rep(x=chr, times=edge.len), 
                          maxgap=-1L, minoverlap=1L, type="any")
  
  edge.TF <- is.na(NETWRK$edges$label) & 1:edge.len%in%hits.mx[,"query"]
  width <- min(NETWRK$edges$width)
  NETWRK$edges$width[edge.TF] <- width*3L
 
  return(NETWRK)
  
}
################################################################################

# rm(list=ls()); gc()