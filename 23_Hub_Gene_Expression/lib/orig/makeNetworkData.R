################################################################################
# Network representation of chromosome substructures formed by contacts. Added 
# feature of marking areas overlapping with a feature provided as a bed file. 
# Option to output network representation of chr with and without contaacts or
# connections.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# PERSIST.MX 

# library(RColorBrewer)
# library(visNetwork)
# source(paste0(wk.dir, "/lib/modifyNetworkData.R"))
# library(compiler)
# library(GenomicRanges)
### FUNCTION ###################################################################
makeNetworkData <- function(
  # i-j contact matrix
  chr = "chr",
  contact.mx = PERSIST.MX$hits,
  cp.v = PERSIST.MX$ntis,
  bin.last = 1204,
  bin.len = 40000,
  edge.res = 50,
  # Multiplier to length of edges to make the difference in lengths more obvious
  edgelen.mult = 2L,
  #-------------------
  # Bed file format (chr, start, end). Ranges will not be reduced to make them
  # non-overlapping. Ideally, each range should corrrespond to a unique feature/marker.
  # Internally, row numbers are added as uniqueID of the feature (before removing
  # ranges with NA start and end coordinates to make sure that the row numbers
  # correspond to the original bed file provided.)
  mark.df = NULL, #mark.df[mark.df[,"chr"]==chr,], 
  bed.CS="1-based", # 1-based or 0-based
  header.bed = TRUE,
  olap.col = "#730707",
  # Add extra bins from mark.df?
  addNodes = TRUE,
  # Marking centromores
  centrobed.file = NULL #'bed file of centromere position; If NULL, do not mark'
){
  
  NETWRK <- list()
  
  if( nrow(contact.mx)==0 ){
    
    ubins.ij <- NULL
    print(paste0(chr, ": Network with no contacts/connections."), quote=FALSE)
    
  } else {
    
    # EDGES - contacts
    uniq.cp <- unique(cp.v)
    #col.v <- tail(brewer.pal(n=5, name="Reds"), n=length(uniq.cp))
    #names(col.v) <- sort(as.numeric(uniq.cp), decreasing=FALSE)
    col.v <- brewer.pal(n=5, name="Reds")
    names(col.v) <- 17:21
    NETWRK$edges <- data.frame(#id=1:length(contact.mx[,1]), 
                               label=cp.v, title=cp.v, from=contact.mx[,1], 
                               to=contact.mx[,2], length=edgelen.mult, width=20, 
                               color=NA, font="60px arial", arrows=NA,
                               stringsAsFactors=FALSE) 
    NETWRK$edges$color <- col.v[as.character(NETWRK$edges$label)]
    rm(contact.mx, cp.v, uniq.cp, col.v); gc()
    # Unique bins of persistent contacts
    ubins.ij <- sort(as.numeric(
      unique(c(NETWRK$edges$from, NETWRK$edges$to))
    )) 
    
    print(paste0(chr, ": Network with contacts/connections."), quote=FALSE)
    
  } 
  
  #---------------------------------------NODES
  
  # Add filler nodes based on edge.res. If, edge.res=5, there will be 5 bins
  # in between two non-contact nodes. 
  ubins.ext <- seq(from=1, to=bin.last, by=edge.res+1L)
  ubins <- sort( unique(c(1, ubins.ij, ubins.ext, bin.last)), decreasing=FALSE )
  ubins.len <- length(ubins)
  
  # Nodes final dataframe
  NETWRK$nodes <- data.frame(id=ubins, label=NA, title=as.character(ubins), 
                             shape="dot", size=12, color="black", 
                             font=NA, borderWidth=0, stringsAsFactors=FALSE)
  
  if( !is.null(ubins.ij) ){
    
    # Special format for persistent contact nodes
    ubins.ij.TF <- ubins%in%unique(c(1, ubins.ij, bin.last))
    NETWRK$nodes$color[ubins.ij.TF] <- adjustcolor("#d93229", alpha.f=0.5)
    rm(ubins.ij.TF); gc()
    
  }
  
  # First and last bin black
  node.last <- nrow(NETWRK$nodes)
  NETWRK$nodes$shape[c(1, node.last)] <- "circle"
  NETWRK$nodes$label[c(1, node.last)] <- "----"
  NETWRK$nodes$color[c(1, node.last)] <- "#ba9c4a"
  NETWRK$nodes$font[c(1, node.last)] <- "14px arial white"
  rm(node.last)
  
  #---------------------------------------
  # EDGES - connect consecutive bins
  
  # -1 to from and to bins to not overlap with nodes
  base <- data.frame(label=NA, title=NA, from=ubins[-(ubins.len)],
                     to=ubins[-1], length=(diff(ubins)-1)*edgelen.mult, 
                     width=15, color=adjustcolor("#aaadad", alpha.f=0.5),
                     font=NA, arrows=NA, stringsAsFactors=FALSE)
  
  if( !is.null(ubins.ij) ){
    NETWRK$edges <- rbind(NETWRK$edges, base)
  } else {
    NETWRK$edges <- base
  }

  NETWRK$edges$font[!is.na(NETWRK$edges$label)] <- "40px arial"
  NETWRK$edges <- NETWRK$edges[order(NETWRK$edges$from, NETWRK$edges$to, decreasing=FALSE),]
  NETWRK$edges$arrows[NETWRK$edges$to%in%c(ubins.ext, bin.last) & is.na(NETWRK$edges$label)] <- "to"
  
  #---------------------------------------
  
  if( !is.null(mark.df) ){
    NETWRK <- modifyNetworkData(
      chr=chr,
      NETWRK=NETWRK,
      bin.len=bin.len,
      mark.df=mark.df, 
      bed.CS=bed.CS,
      header.bed=header.bed,
      olap.col=olap.col,
      addNodes=addNodes,
      edgelen.mult=edgelen.mult
    )
  }
  
  # Mark centromere
  if( !is.null(centrobed.file) ){
    NETWRK <- markCentromereOnNetworkData(NETWRK=NETWRK, chr=chr, bin.len=bin.len,
                                          centrobed.file=centrobed.file)
  }
  
  return(NETWRK)
  
}
################################################################################
makeNetworkData <- cmpfun(makeNetworkData, options=list(suppressUndefined=TRUE))
################################################################################
