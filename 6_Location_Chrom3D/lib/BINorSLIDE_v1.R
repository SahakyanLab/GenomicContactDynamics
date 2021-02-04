################################################################################
BINorSLIDE <- function(numVec=DTA$radPos, 
                       action="SLIDE",
                       numPoints=950,
                       # radius of bin or sliding window
                       dr=0.5){
  
  if(!is.numeric(numVec)){
    stop("Only numeric vector is accepted.")
  }
  
  if(action=="BIN"){
    
    rFirst <- min(numVec, na.rm=TRUE) 
    rLast <- max(numVec, na.rm=FALSE) 
    rVal = unique( c( seq(from=rFirst, to=rLast, by=2*dr), rLast) )
    
  } else if(action=="SLIDE"){
    
    # ± dr makes sures that the first and last window are within the data range
    # consequently, the first and last r is dr distance away from min and max value, respectively
    # in addition, the distance between the values within rFirst and rLast can be very small 
    # compared to dr depending on the number of points wanted
    rFirst <- min(numVec) + dr
    rLast <- max(numVec) - dr
    rVal <- seq(from=rFirst, to=rLast, length.out=numPoints-2)
    if( !(rVal[2]-rVal[1]<2*dr) ){
     warning("dr is smaller than difference between r values. 
             Radial windows not overlapping.") 
    }
    
    return(rVal)
    
  } else {
    
    stop(' action can only be "BIN" or "SLIDE". ')
    
  }
  return(rVal)
}


