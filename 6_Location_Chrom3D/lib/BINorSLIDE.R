################################################################################
# Function to get N points across the range of the radial distance such that
# windows formed by point±dr is overlapping (SLIDE) or non-overlapping (BIN)
################################################################################
BINorSLIDE <- function(numVec=DTA$radDist, 
                       action="SLIDE",
                       numPoints=1000, # For action="SLIDE" only
                       # Radius of bin or sliding window
                       dr=0.5){
  
  if(!is.numeric(numVec)){
    stop("Only numeric vector is accepted.")
  }
  
  #rLast <- ceiling(max(numVec))
  rLast <- max(numVec) 
  
  if(action=="BIN"){
    
    rFirst <- min(numVec, na.rm=TRUE) 
    # Points separated by 2*dr 
    rVal = unique( c( seq(from=rFirst, to=rLast, by=2*dr), rLast) )
    
  } else if(action=="SLIDE"){
    
    rVal <- seq(from=0, to=rLast, length.out=numPoints)
    # Check if given the rVal values and dr, the windows will overlap
    if( !(rVal[2]-rVal[1]<2*dr) ){
     warning("dr is smaller than difference between r values. 
             Radial windows not overlapping.") 
    }
    
    return(rVal)
    
  } else {
    
    stop('Action can only be "BIN" or "SLIDE".')
    
  }
  return(rVal)
}


