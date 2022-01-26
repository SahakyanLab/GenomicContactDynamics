################################################################################
# Function to split vector x into non-overlapping (BIN) and overlapping intervals.
# BIN function will split f such that each bin will contain roughly equal 
# percentage of datapoints. Use R base function cut() if you want to supply 
# breaks. I didn't incorporate this here because using cut() directly will allow
# control of type of interval (closed which side?). Output will be a list of chunks
# containing indices of values on x vector.
################################################################################
splitNumericVector <- function(x=c(1,2,4,6), 
                               # If action=BIN, d is the percentage of points per
                               # boxplot. If action=SLIDE, it is the radius of the
                               # sliding window.
                               d=0.5,
                               action="SLIDE", # BIN | SLIDE
                               # For action="SLIDE" only
                               numSLwind=1000, # Number of windows 
                               # Ideally, max and min values of x
                               boundsSL=c(0, max(numVec)) 
                              ){
  
  if( !is.numeric(x) ){
    stop("x should be numeric.")
  }
  
  if( any(!is.finite(x)) ){
    stop("x has non-finite values.")
  }
                      
  if( is.unsorted(x) ){
    warning("x is not sorted. This is needed for action=SLIDE.")
  }
  
  #-------------------
  if(action=="BIN"){
    
    lst <- ceiling(seq_along(x)/(length(x)*d/100)) 
    names(lst) <- 1:length(lst)
    
  } else if(action=="SLIDE"){
    
    # Midpoints of overlapping windows; ±d to boundsSL to make sure that
    # span of windows won't exceed limits set by boundsSL.
    mid.v <- seq(from=boundsSL[1]+d, to=boundsSL[2]-d, length.out=numSLwind)
    # Given mid.v and d, check if the windows will overlap
    if( !(mid.v[2]-mid.v[1]<2*d) ){
     warning("d (the radius of the overlapping window) is smaller than 
              difference between r values. Radial windows not overlapping.") 
    }
    
    lst <- sapply(X=mid.v, simplify=FALSE, FUN=function(mid){
      bound <- c(mid-d, mid+d)
      return( which(x>=bound[1] & x<=bound[2]) )
    })
    names(lst) <- 1:length(lst)
    lst[["midpoints"]] <- mid.v
    
  } else {
    stop('Action can only be "BIN" or "SLIDE".')
  }
  
  return(lst)

}

# rm(list=ls()); gc()

