################################################################################
# For making metaplots
# Dependencies:
## library(ggplot2)
## source(paste0(lib, "/GG_bgr.R"))
################################################################################
################################################################################
convertToIncrements <- function(vec = FETACP.MX[1,],
                                refval = FETACP.MX[1,2201]){
  names(refval) <- NULL
  
  len <- length(vec)
  vec.inc <- sapply(X=1:len, simplify=TRUE, FUN=function(x){
    return(vec[x]/refval)
  })
  
  return(vec.inc)
}

