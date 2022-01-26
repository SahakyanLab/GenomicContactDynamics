################################################################################
# Function to categorise values into 3 categories top X% highest values (1), 
# bottom X% lowest values (-1) and the rest (0). Note that the X% cut-off value 
# was determined using only non-NA values from those supplied.
### FUNCTION ###################################################################
categoriseValues <- function(val.v='vector of values',
                             cutoff='percentage e.g. 5 for 5%'
                             ){
  
  if( !is.numeric(val.v) | !is.numeric(cutoff) ){
    stop("categoriseCIIscore(): Supplied values or cutoff not numeric.")
  }
  
  if( !identical( is.na(val.v), !is.finite(val.v) ) ){
    stop("categoriseCIIscore(): There are other non-finite values aside from NA.")
  }
  
  val.v <- as.numeric(val.v)
  #val.srt <- val.v[!is.na(val.v)]
  # na.last=NA to remove NA
  val.srt <- sort(val.v, decreasing=T, na.last=NA)
  len <- length(val.srt)
  
  # Determine boundaries of categories
  
  # Ideal length of high and low categories based on cutoff percentage
  cutval <- ceiling(len*cutoff/100)
  # High category minimum value
  up.b <- val.srt[cutval] 
  # Low category maximum value
  low.b <- val.srt[len-cutval+1] # len-cutval was used for the chr1 simulation data
  
  if(low.b>=up.b){
    stop("categoriseCIIscore(): Boundaries of high and low categories cross over or
         are equal.")
  }
  
  # Categorise based on boundary values
  
  group <- val.v
  group[val.v > (up.b)] <- 1
  group[(val.v <= (up.b)) & (val.v >= (low.b))] <- 0
  group[val.v < low.b] <- -1
  
  if( !identical( is.na(val.v), is.na(group) ) ){
    stop("categoriseCIIscore(): Not matching.")
  } 
  
  return(group)
  
}

# rm(list=ls()); gc()
