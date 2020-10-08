################################################################################
# Trapezoid rule to integrate and get area/volume under curve/surface
### FUNCTION ###################################################################
getVUS <- function(X, Y){
  
  # Shift Y to have minimum at 0
  Y <- Y-min(Y)
  
  X <- as.matrix(X)
  grid.col <- ncol(X)
  grid.row <- nrow(X)

  # Checks
  if( nrow(X)!=length(Y) ){
    stop("Lengths of argument vectors not equal")
  }
  
  df <- cbind.data.frame(X,Y)
  # Order y 
  for(i in 1:grid.col){ df <- df[order(df[[i]]),] }
  
  X <- lapply(df[,-(grid.col+1)], FUN=function(col){
    diff(sort(unique(col)))
  })
  X <- apply(X=expand.grid(X), MARGIN=1, FUN=prod)
  
  # Calculate VUS per chunk
  chunk.len <- length(unique(df$Var1))
  nchunk <- grid.row/chunk.len-1
  chunk.id <- rep(1:nchunk, each=chunk.len)
  
  sumBase <- sapply(X=1:nchunk, simplify=FALSE, function(n){
    y <- df$Y[chunk.id==n]
    y <- y[-chunk.len]+y[-1] 
    return(y)
  })
  sumBase <- unlist(sumBase)
  
  if( length(X)!=length(sumBase) ){
    stop("Final X and sumBase lengths differ.")
  }

  VUS <- sum(0.5*X*sumBase)
  
  return(VUS)
  
}
################################################################################
