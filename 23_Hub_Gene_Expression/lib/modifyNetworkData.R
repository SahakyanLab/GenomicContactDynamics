################################################################################
# Make network data out of  contact data (PERSIST.MX)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# NETWRK
# library(compiler)
# library(GenomicRanges)
# source(paste0(wk.dir, "/lib/addNodesToNetwork.R"))
### FUNCTION ###################################################################
modifyNetworkData <- function(
  chr = "chr",
  NETWRK = NULL,
  bin.len = 40000,
  # Multiplier to length of edges to make the difference in lengths more obvious
  edgelen.mult = 2L,
  # Bed file format (chr, start, end, name, score). Ranges will not be reduced to make them
  # non-overlapping. Ideally, each range should correspond to a unique feature/marker.
  # Internally, row numbers are added as uniqueID of the feature (before removing
  # ranges with NA start and end coordinates to make sure that the row numbers
  # correspond to the original bed file provided.)
  mark.df = mark.df[mark.df[,1]==chr,], 
  bed.CS="1-based", # 1-based or 0-based
  header.bed = FALSE,
  olap.col = "#730707",
  # Add extra bins from mark.df?
  addNodes = TRUE,
  # Marking centromores
  centrobed.file = NULL, #'bed file of centromere position'
  colouring.style = "discrete" # "binary" | "discrete" | "gradient"
){
  
  if(bed.CS == "0-based"){
    mark.df[,2] <- mark.df[,2] + 1L
    print("Bed file converted to 1-based.", quote=FALSE)
  } else if (bed.CS == "1-based"){
    print("Bed file already one-based.", quote=FALSE)
  } else {
    stop("Invalid bed.CS input ('1-based' or '0-based' only).")
  }

  col.len <- ncol(mark.df)
  row.len <- nrow(mark.df)
  if(col.len >= 5){
    mark.df <- mark.df[,1:5]
  } else if(col.len == 4){
    mark.df$score.val <- 1
  } else if(col.len == 3){
    mark.df$name.feat <- paste0("R_", 1:row.len)
    mark.df$score.val <- 1
  } else {
    stop("modifyNetworkData(): Not enough columns for a bed file.")
  }

  mark.df[,2] <- ceiling(mark.df[,2] / bin.len)
  mark.df[,3] <- ceiling(mark.df[,3] / bin.len)
  colnames(mark.df) <- c("chr", "bin.from", "bin.to", "name.feat", "score.val")
  mark.df$uniqueID <- 1:row.len
  
  #-------------------
  
  # ADD NODES from mark.df and update edges to include new nodes
  if(addNodes){
    NETWRK <- addNodesToNetwork(mark.df=mark.df, NETWRK=NETWRK, edgelen.mult=edgelen.mult)
    NETWRK <- markCentromereOnNetworkData(NETWRK=NETWRK, chr=chr, bin.len=bin.len, centrobed.file=centrobed.file)
  }
  
  #-------------------
  
  # IDENTIFY NUMBER OF MARKERS OVERLAPPING (num.olap) 
  for(x in c("nodes", "edges")){
    
    if(x == "nodes"){
      # Overlaps of marks with nodes 
      start.q <- end.q <- NETWRK$nodes$id
      edge.cp.ind <- NULL
    } else {
      
      # Overlaps of marks with edges. Edges are denoted by the two nodes they should
      # connect so to avoid treating edges to have an overlap just because of their 
      # nodes, we do not include the nodes in counting the overlaps of edges with
      # the markfile (hence the +1 in start bin and -1 bin in end bin of edges). 
      # Do this only for edges with length greater than 1, otherwise, those edges
      # will have negative length. 
      start.q <- NETWRK$edges$from
      end.q <- NETWRK$edges$to
      adjustby1.TF <- (NETWRK$edges$to - NETWRK$edges$from) > 1
      start.q[adjustby1.TF] <- start.q[adjustby1.TF] + 1
      end.q[adjustby1.TF] <- end.q[adjustby1.TF] - 1
      rm(adjustby1.TF)
      # For edges connecting highlighted contacts, overlap counting will be separate
      # because we want to count overlaps WITHIN the region they span (the loop they form)
      edge.cp.ind <- which(!is.na(NETWRK$edges$label))
      
    }
    
    hits.df <- WhichOverlap(start.query=start.q, 
                            end.query=end.q, 
                            space.query=rep(x=chr, times=length(start.q)),
                            start.subject=as.numeric(mark.df$bin.from),
                            end.subject=as.numeric(mark.df$bin.to),
                            space.subject=mark.df$chr, 
                            maxgap=-1L, minoverlap=1L, type="any")
    colnames(hits.df) <- c("netwrk.ind", "marker.ind")
    
    NETWRK[[x]][,"num.olap"] <- 0
    NETWRK[[x]][,"marker.olap"] <- NA
    NETWRK[[x]][,"name.feat.olap"] <- NA
    NETWRK[[x]][,"score.olap.mean"] <- NA
    
    if( nrow(hits.df) == 0){
      print(paste0("No overlap with ", x, "!"), quote=FALSE)
      next
    }
    
    if(x == "edges" & length(edge.cp.ind) > 0){
      
      # Which of the markers are within the loops formed by the highlighted contacts?
      # Query and subject swapped
      edge.cp.hits.df <- WhichOverlap(start.query=as.numeric(mark.df$bin.from), 
                                      end.query=as.numeric(mark.df$bin.to), 
                                      space.query=mark.df$chr,
                                      start.subject=start.q[edge.cp.ind],
                                      end.subject=end.q[edge.cp.ind],
                                      space.subject=rep(x=chr, times=length(edge.cp.ind)), 
                                      maxgap=-1L, minoverlap=1L, type="within")
      # netwrk.ind here in reference to edge.cp.hits.df
      colnames(edge.cp.hits.df) <- c("marker.ind", "netwrk.ind")
      # netwrk.ind here in reference to hits.df
      edge.cp.hits.df[,"netwrk.ind"] <- edge.cp.ind[edge.cp.hits.df[,"netwrk.ind"]]
      # Replace edge.cp hits in hits.df with edge.cp.hits.df; note that query and
      # subject swapped. 
      # query in hits.df, netwrk.ind in edge.cp.hits.df
      # subject in hits.df, marker.ind in edge.cp.hits.df
      hits.df <- hits.df[!hits.df[,"netwrk.ind"] %in% edge.cp.ind,]
      rm(edge.cp.ind)
      hits.df <- rbind(hits.df, edge.cp.hits.df[,c("netwrk.ind","marker.ind")])
      hits.df <- hits.df[order(hits.df[,"netwrk.ind"]),]
      rm(edge.cp.hits.df)
      
    } 
    
    hits.df <- data.frame(hits.df, 
                          uniqueID=mark.df[,"uniqueID"][hits.df[,"marker.ind"]],
                          name.feat=mark.df[,"name.feat"][hits.df[,"marker.ind"]],
                          score.val=mark.df[,"score.val"][hits.df[,"marker.ind"]],
                          stringsAsFactors=FALSE)
    
    hits.df <- by(data=hits.df[,c("marker.ind", "uniqueID", "name.feat", "score.val")], 
                  INDICES=hits.df[,"netwrk.ind"], 
                  FUN=function(x){
                    
                    cbind.data.frame(
                      num.olap=length(unique(x[,"marker.ind"])),
                      marker.olap=paste(x[,"uniqueID"], collapse=";"),
                      name.feat.olap=paste(x[,"name.feat"], collapse=";"),
                      # na.rm=T so those with overlap but no score will give NaN
                      score.olap.mean=mean(x[,"score.val"], na.rm=T),
                      stringsAsFactors=FALSE
                    )
                    
                  })
    
    netwrk.ind <- as.numeric(names(hits.df))
    if(is.unsorted(netwrk.ind)){ stop("Indices not increasing.") }
    hits.df <- do.call("rbind.data.frame", hits.df)
    if( nrow(hits.df) != length(netwrk.ind) ){
      stop("hits.df rows not equal to number of nodes/edges with overlap.")
    }
    
    NETWRK[[x]][netwrk.ind,colnames(hits.df)] <- hits.df
    
    # Check if nodes/edges with num.olap==0 has no marker.olap (NA)
    if(! identical( NETWRK[[x]][,"num.olap"] == 0, is.na(NETWRK[[x]][,"marker.olap"]) )){
      stop("Inconsistency between num.olap and marker.olap.")
    }
    
    rm(hits.df, netwrk.ind, start.q, end.q); gc()
    
    print(paste0(x, " done!"), quote=FALSE)
    
  } # for loop end
  
  #-------------------
  
  nonij.TF <- is.na(NETWRK$edges$label)
  edgewithbin.TF <- (NETWRK$edges$to - NETWRK$edges$from) > 1
  edgeToColor.TF <- nonij.TF & edgewithbin.TF
  
  # Hover information
  hover.text <- cbind(nodes.title=NETWRK$nodes$title, 
                      nodes.name=NETWRK$nodes$name.feat.olap,
                      nodes.score=NETWRK$nodes$score.olap.mean)
  hover.text[is.na(hover.text)] <- ""
  NETWRK$nodes$title <- paste0(hover.text[,"nodes.title"], "<br>", 
                               "<p style='text-align:left;'>",
                               "Name: ", hover.text[,"nodes.name"], "<br>",
                               "Score: ", hover.text[,"nodes.score"], "</p>")
  hover.text <- cbind(edges.name=NETWRK$edges$name.feat.olap,
                      edges.score=NETWRK$edges$score.olap.mean)[edgeToColor.TF,]
  hover.text[is.na(hover.text)] <- ""
  NETWRK$edges$title[edgeToColor.TF] <- paste0("<p style='text-align:left;'>", 
                                               "Name: ", hover.text[,"edges.name"], "<br>",
                                               "Score: ", hover.text[,"edges.score"], "</p>")
  rm(hover.text)
  
  #-------------------Colour nodes and edges 
  
  if( sum(nonij.TF & edgewithbin.TF) > 0 ){
    NETWRK$edges[nonij.TF & !edgewithbin.TF, "color"] <- "black"
  }
  
  # With overlap but no name or no scores to take the mean of

  col.nme <- ifelse(colouring.style == "discrete", "name.feat.olap", "score.olap.mean")
 
  for( prt in c("nodes", "edges") ){
    
    fnx0 <- ifelse(colouring.style == "discrete", 'is.na', 'is.nan')
    eval(parse(text=paste0(
      'nan.TF <- ', fnx0, '(NETWRK$', prt, '$', col.nme, ') & (NETWRK$', prt, '$num.olap > 0)'
    )))
    if(prt == "edges"){ nan.TF <- nan.TF & edgeToColor.TF }
    if( sum(nan.TF) > 0 ){
      NETWRK[[prt]][nan.TF,"color"] <- "deeppink"
    }
    
  }
  
  # With overlaps and names/scores
  
  vals <- c(NETWRK$nodes[[col.nme]], NETWRK$edges[[col.nme]][edgeToColor.TF])
  fnx <- ifelse(colouring.style == "discrete", '!is.na', 'is.finite')
  eval(parse(text=paste0( 'nonijfinite.all.TF <- ', fnx, '(vals)' )))
  
  if( sum(nonijfinite.all.TF) > 0 ){
    
    nonijfinite.val <- vals[nonijfinite.all.TF]
    
    ## Discrete
    
    if(colouring.style == "discrete"){
      
      print("modifyNetworkData(): Discrete colouring style.")
      
      nonijfinite.val <- vals[nonijfinite.all.TF]
      col.hex <- rep(NA, length=length(nonijfinite.val))
      
      nonijfinite.val <- strsplit(x=nonijfinite.val, split=";")
      nonijfinite.val <- lapply(X=nonijfinite.val, FUN=unique)
      
      col.hex[lengths(nonijfinite.val) > 1] <- "blue" # Multiple categories
      
      grp.TF <- lengths(nonijfinite.val) == 1
      val.grps <- unique(unlist(nonijfinite.val[grp.TF]))
      
      col.pal <- setNames(colorRampPalette(brewer.pal(n=8, name="Dark2"))(length(val.grps)),
                          nm=val.grps)
      col.hex[grp.TF] <- col.pal[unlist(nonijfinite.val[grp.TF])]
      col.hex <- adjustcolor(col.hex, alpha.f=0.7)
      
    } else if( colouring.style %in% c("binary", "gradient") ){
      
      print("modifyNetworkData(): Binary/Gradient colouring style.")
      
      if(colouring.style == "binary"){ nonijfinite.val <- rep(1, length(nonijfinite.val)) }
    
      min.val <- min(nonijfinite.val)
      col.val <- nonijfinite.val - min.val # Convert to > 0 values
      col.val <- (col.val - min(col.val))
      denom <- diff(range(col.val))
      if(denom > 0){
        col.val <- col.val / denom
      }
      col.pal <- colorRamp(viridis(5))
      col.hex <- adjustcolor( rgb(col.pal(col.val), max=255), alpha.f=0.7 )
      
    } else {
      stop("modifyNetworkData(): Invalid colouring.style argument.")
    }
    
    for( prt in c("nodes", "edges") ){
      
      eval(parse(text=paste0( 'nonijfinite.TF <- ', fnx, '(NETWRK$', prt, '$', col.nme, ')' )))
      if(prt == "edges"){ nonijfinite.TF <- nonijfinite.TF & edgeToColor.TF }
      nonijfinite.len <- sum(nonijfinite.TF)
      if(nonijfinite.len == 0){
        next
      } else {
        
        NETWRK[[prt]][nonijfinite.TF,"color"] <- col.hex[1:nonijfinite.len]
        col.hex <- col.hex[-(1:nonijfinite.len)]
 
      }
      
    }
    
  } # if end
  
  return(NETWRK)
  
}
################################################################################
modifyNetworkData <- cmpfun(modifyNetworkData, options=list(suppressUndefined=TRUE))
################################################################################
