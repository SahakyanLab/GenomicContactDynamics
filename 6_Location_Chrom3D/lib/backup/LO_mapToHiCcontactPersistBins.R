################################################################################
# Use LiftOver on Hi-C contact persist bins
# DEPENDENCIES
## lib = "/Users/ltamon/DPhil/lib"
## library(rtracklayer)
## library(ggplot2)
## source(paste0(lib, "/TrantoR_liftOver/GEN_liftOveR.R"))
## source(paste0(lib, "/TrantoR_liftOver/GEN_liftOverLoadChain.R"))
## source(paste0(lib, "/plotLengthDist.R"))
## source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
LO_mapToHiCcontactPersistBins <- function(out.dir = "/dir", 
                                          out.name = paste0(chr, "_Persist_", gcb),
                                          chr = "chr", 
                                          start = bin.start, end = bin.end,
                                          LOchain = "hg19ToHg38", rmchain = FALSE
                                          ){
  
  len <- length(end)
  LO.out <- liftOveR(conversion=LOchain, 
                     space=rep(x=chr, times=len), 
                     start=as.numeric(start),
                     end=as.numeric(end),
                     strand=rep(x="*", times=len),
                     getchain=FALSE,
                     rmchain=rmchain,
                     returnGRangesList=FALSE)
  LO.out <- LO.out[LO.out$seqnames==chr,]
  # Make sure ranges come from the same chr
  if(unique(LO.out$seqnames)!=chr){stop("Liftover ranges from different chr.")}
  # For regions broken into pieces, just get the coordinates spanning all pieces
  # Pieces may be not in order (lowest to highest)
  newRanges <- by(data=LO.out[,c("start", "end", "group")], INDICES=LO.out$group, FUN=function(x){
    xlen <- nrow(x)
    return( c(start=min(x[,"start"]), end=max(x[,"end"])) )
  })
  rm(LO.out); gc()
  bin.conv <- as.numeric(names(newRanges))
  newRanges <- do.call("rbind", newRanges)
  # Plot lengths of converted ranges   
  plotLengthDist(df=cbind.data.frame(variable="LO",
                                     value=log10(newRanges[,"end"]-newRanges[,"start"])),
                 vline.v=log10(40000),
                 col="deepskyblue3",
                 out.name=paste0(out.name, "_", length(bin.conv), "bins_", LOchain, "_LO"),
                 out.dir=out.dir,
                 label.x="log10 Length")
  # Order of ranges should correspond to order of bins
  if( identical(1:len, bin.conv) ){
    print(paste0("All bins of ", chr, " converted."))
    return(newRanges)
  } else {
    # For bins that cannot be converted, populate start and end with NA
    mx <- matrix(data=NA, nrow=len, ncol=2, dimnames=list(NULL, c("start", "end"))
                 )
    mx[bin.conv,] <- newRanges
    
    print(paste0(len-length(bin.conv), " bin/s of ", out.name, " not converted."), quote=FALSE)
    return(mx)
  }

} # Function end
