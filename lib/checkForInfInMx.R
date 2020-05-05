# Check for NA/NaN/Inf in matrices 

checkForInfInMx <- function(mx = ELMTISSDYN.MX.norm){
  v <- apply( X=mx, MARGIN=1, FUN=function(ntiscol){sum(is.finite(ntiscol))} )
  return(which(v==FALSE))
}