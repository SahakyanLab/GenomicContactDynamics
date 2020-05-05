###############################################################################
# reformat cmm file from Chrom3D to isolate x,y,z coordinates and radii of domains
# works regardless of ploidy
###############################################################################
###############################################################################
getcmmXYZR <- function(cmmPath = cmmfile){
  cmm <- readLines(cmmPath)
  # remove <marker_set name="chrom3d_model"> (first line)
  cmm <- cmm[-1] 
  # remove "</marker_set>" (last line)
  cmm <- cmm[-length(cmm)]
  # remove linkers
  cmm <- cmm[-(grep(pattern="<link id", 
                    x=cmm, ignore.case=TRUE))]
  DOMAIN.XYZR.MX <- do.call("rbind", 
                            strsplit(x=cmm, split='"'))[,c(18,4,6,8,10,20)]
  DOMAIN.XYZR.MX <- cbind.data.frame(DOMAIN.XYZR.MX, stringsAsFactors=FALSE)
  colnames(DOMAIN.XYZR.MX) <- c("chr", "x", "y", "z", "cmmRadius", "id")
  
  DOMAIN.XYZR.MX$chr <- as.character(DOMAIN.XYZR.MX$chr)
  DOMAIN.XYZR.MX$id <- as.character(DOMAIN.XYZR.MX$id)
  DOMAIN.XYZR.MX$x <- as.numeric(DOMAIN.XYZR.MX$x)
  DOMAIN.XYZR.MX$y <- as.numeric(DOMAIN.XYZR.MX$y)
  DOMAIN.XYZR.MX$z <- as.numeric(DOMAIN.XYZR.MX$z)
  DOMAIN.XYZR.MX$cmmRadius <- as.numeric(DOMAIN.XYZR.MX$cmmRadius)
  
  # check if there are domains duplicated
  if(sum(duplicated(DOMAIN.XYZR.MX$id))!=0){
    stop("Duplicated domain ids.")
  }
 return(DOMAIN.XYZR.MX)
}

# rm(list=ls())
