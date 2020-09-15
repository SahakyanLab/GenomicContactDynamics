################################################################################
# Custom evaluation functions for the permutation test. Returns a named list called
# eval.f.lst of all functions here. 
################################################################################
eval.f.lst <- list(
  numOlapA = function(A, B, ...) {
    num <- regioneR::numOverlaps(A=A, B=B, count.once=TRUE)
    #return( num/length(A)*100 )
    return(num)
  },
  numOlapAwithin = function(A, B, ...) {
    num <- length(
      subsetByOverlaps(x=A, ranges=B, maxgap=-1L, minoverlap=1L, type="within",
                       invert=FALSE)
    )
    #return( num/length(A)*100 )
    return(num)
  },
  numOlapB = function(A, B, ...) {
    num <- regioneR::numOverlaps(A=B, B=A, count.once=TRUE)
    #return( num/length(B)*100 )
    return(num)
  },
  numOlapBwithin = function(A, B, ...) {
    num <- length(
      subsetByOverlaps(x=B, ranges=A, maxgap=-1L, minoverlap=1L, type="within",
                       invert=FALSE)
    )
    #return( num/length(B)*100 )
    return(num)
  },
  comOlap = function(A, B, ...) {
    return(sum( width(commonRegions(A,B)) )
    )
  },
  #comOlapA = function(A, B, ...) {
  #  total <- sum( width(A) )
  #  common <- sum( width(commonRegions(A,B)) )
  #  return( common/total*100 )
  #},
  #comOlapB = function(A, B, ...) {
  #  total <- sum( width(B) )
  #  common <- sum( width(commonRegions(A,B)) )
  #  return( common/total*100 )
  #},
  meandist = function(A, B, ...) {
    return( -(regioneR::meanDistance(A,B)) )
  },
  meanLenOlapA = function(A, B, ...) {
    x <- mean(width(
      subsetByOverlaps(x=A, ranges=B, maxgap=-1L, minoverlap=1L, type="any", 
                       invert=FALSE)
    ))
    return(x)
  }, 
  meanLenOlapB = function(A, B, ...) {
    x <- mean(width(
      subsetByOverlaps(x=B, ranges=A, maxgap=-1L, minoverlap=1L, type="any", 
                       invert=FALSE)
    ))
    return(x)
  }
)

return(eval.f.lst)
