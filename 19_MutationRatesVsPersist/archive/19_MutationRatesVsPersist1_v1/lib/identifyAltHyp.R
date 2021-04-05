############################################################################### 
# Determine appropriate alternative hypothesis by comparing two statistic of a 
# distribution (e.g. mean or median).
###############################################################################
# "greater" means x > y
identifyAltHyp <- function( x=c(1,2), y=c(3,4) ){
  
  if( length(x)!=length(y) ){
    stop("Lengths of x and y vectors are different.")
  }
  
  alt.v <- x-y
  alt.v <- alt.v/abs(alt.v)
  
  if( !all(alt.v%in%c(-1, 1, NaN)) ){
    stop("Checkpoint 1.")
  }
  
  alt.v <- c(`-1`="less", `1`="greater","NaN"="two.sided")[as.character(alt.v)]
  names(alt.v) <- c("MEAN", "MEDIAN")
  return(alt.v)
  
}
