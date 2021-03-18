################################################################################
# Determine overlap of supplied coordinates with supplied bed files.
# In output, 1 means overlapping, 0 means not overlapping, NA means coordinate is
# missing.

# NOTE:
# GEN_WhichOverlap.R function used ignores strand information. Modify if needed.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(GenomicRanges)
# source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
### FUNCTION ###################################################################
whichOverlapToBed <- function(bed.dir, header, seqnames, start, end, strand){
                                 
  #-------------------Query coordinates
  coord.gr <- data.frame(seqnames=seqnames, start=start, end=end, 
                         strand=as.character(strand), stringsAsFactors=FALSE)
  
  #-------------------Bed files
  bedNme.v <- list.files(bed.dir, full.names=FALSE)
  bed.gr <- list()  
  bed.str.v <- c()
  for(bedNme in bedNme.v){
    
    bed.gr[[bedNme]] <- read.table(file=paste0(bed.dir, "/", bedNme), header=header, 
                                   sep="\t", stringsAsFactors=FALSE)
    if( ncol(bed.gr[[bedNme]]) >= 6 ){
      
      print(paste0(bedNme, " have strand information."))
      bed.str.v <- c( bed.str.v, bed.gr[[bedNme]][,6] )
      
    } else {
      
      print(paste0(bedNme, " has no strand information. Defaulting to '*'."))
      bed.str.v <- c(
        bed.str.v, rep( "*", times=length(bed.gr[[bedNme]][,1]) )
      )
      
    }
    
    # Take only chr, start and end
    bed.gr[[bedNme]] <- bed.gr[[bedNme]][,1:3]
    
  } # bedNme.v for loop end
  
  bed.gr <- do.call("rbind.data.frame", c(bed.gr, stringsAsFactors=FALSE))
  colnames(bed.gr) <- c("seqnames", "start", "end")
  rownames(bed.gr) <- NULL
  if( length(bed.str.v)!=length(bed.gr[,1]) ){
    stop("Unequal length of beds' columns and strand info.")
  }
  bed.gr <- data.frame(bed.gr, strand=as.character(bed.str.v), stringsAsFactors=FALSE)
  rm(bed.str.v); gc()
  
  #-------------------
  # Do checks for mimssing values 
  for( x in c("bed", "coord") ){
    
    # Check for NAs
    eval(parse(text=paste0(
      'NA.TF <- apply(X=', x, '.gr, MARGIN=2, FUN=is.na)'
    )))
    
    if( sum(NA.TF)>0 ){
      warning(paste0("Missing values in supplied ", x, "."))
    }
    
    # Start and end NAs identical?
    if( !identical(NA.TF[,2], NA.TF[,3]) ){
      
      stop(paste0(x, ": Missing values of start and end different."))
      
    } else {
      
      eval(parse(text=paste0(
        x, 'Incl.TF <- !NA.TF[,2]; rm(NA.TF)'
      )))
      
    }
    
    # Check strands
    eval(parse(text=paste0(
      'str <- ', x, '.gr$strand'
    )))
    if( any(!str%in%c("+", "-", "*")) ){
      stop(paste0(x, ": Invalid strand."))
    }
    rm(str)
  
  }
  
  #-------------------Determine coordinates that overlap
  out.v <- rep( 0, times=length(coord.gr[,1]) )
  out.v[!coordIncl.TF] <- NA
  
  bed.gr <- bed.gr[bedIncl.TF,]
  coord.gr <- coord.gr[coordIncl.TF,]
  
  olap <- WhichOverlap(start.query=coord.gr$start, end.query=coord.gr$end, 
                       space.query=coord.gr$seqnames,
                       start.subject=bed.gr$start, end.subject=bed.gr$end,
                       space.subject=bed.gr$seqnames, 
                       maxgap=-1L, minoverlap=1L, type="any")
  
  # Coordinates overlapping with beds
  out.v[ which(coordIncl.TF)[as.numeric(olap[,"query"])] ] <- 1
  
  return(out.v)
  
}
