################################################################################
# Trapezoid rule for double integration given fix step size along each dimension
# Reference: http://www.ohiouniversityfaculty.com/youngt/IntNumMeth/lecture24.pdf
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(reshape2)
### FUNCTION ###################################################################
getVUS <- function(X, Y, step.v=c(0.01,0.01)){
  
  # Shift Y to have minimum at 0
  Y <- Y-min(Y)
  
  # Checks
  if( nrow(X)!=length(Y) ){
    stop("Lengths of argument vectors not equal")
  }
  
  # Value and weight matrix
  w.mx <- mx <- reshape2::acast(data=cbind.data.frame(X, Y), formula=Var2~Var1)
  w.dim <- dim(w.mx)
  # Interior weights
  w.mx[!is.na(w.mx)] <- 4
  # Edge and corner weights
  w.mx[,1] <- w.mx[,w.dim[2]]  <- c(1, rep(2, times=w.dim[1]-2), 1)
  w.mx[1,] <- w.mx[w.dim[1],]  <- c(1, rep(2, times=w.dim[2]-2), 1)
  
  # Calculate step size for each dimension
  #d <- dimnames(mx)
  #d <- lapply(X=d, FUN=function(x) unique(diff(as.numeric(x))) )
  
  # Calculate VUS estimate
  sm <- sum(w.mx*mx)
  VUS <- (prod(step.v))/4*sm
  
  return(VUS)
  
}
################################################################################
