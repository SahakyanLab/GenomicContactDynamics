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
                                          chr = NULL, bin.len = NULL, bin.num = NULL,
                                          start = NULL , end = NULL,
                                          LOchain = NULL, rmchain = FALSE
                                          ){
  len <- length(end)
  # Note that a bin can correspond to multiple ranges, even from other chr
  LO.mx <- liftOveR(conversion=LOchain, 
                    space=rep(x=chr, times=len), 
                    start=as.numeric(start),
                    end=as.numeric(end),
                    strand=rep(x="*", times=len),
                    getchain=FALSE,
                    rmchain=rmchain,
                    returnGRangesList=FALSE)
  if(nrow(LO.mx)==0){ stop("No bins successfully converted.") }
  # Convert LO.mx$group to actual bin number
  LO.mx$group <- bin.num[LO.mx$group]
  save(LO.mx, file=paste0(out.dir, "/", out.name, "_", LOchain, "_RAW.RData"))
  # Remove bins with at least 1 converted range from other chr
  LO.mx <- LO.mx[!LO.mx$group%in%LO.mx[LO.mx$seqnames!=chr,"group"],]
  if(unique(LO.mx$seqnames)!=chr){stop("Liftover ranges from different chr.")}
  
  tot.len.range <- aggregate(x=LO.mx$width, by=list(LO.mx$group), FUN=sum)[,"x"]
  # Plot total length of converted ranges per bin   
  plotLengthDist(df=cbind.data.frame(variable="LO", value=tot.len.range/10000),
                 vline.v=bin.len/10000,
                 col.v="deepskyblue3",
                 out.name=paste0(out.name, "_", length(unique(LO.mx$group)), "bins_", LOchain, "_LO"),
                 out.dir=out.dir,
                 label.x="Length/10000")

  return(LO.mx)
} # Function end
