################################################################################
# Convert contact values to contact probability by dividing by max value per bin
### FUNCTION ###################################################################
convertToContactProb <- function(df = 'data.frame(x, y, value)', # upper triangle only
                                 metric,
                                 tot.bin = 'total bins in df; if bins not sequential
                                 i and j values in returned df correspond to index
                                 of actual bins when arranged based on this
                                 order(df$j, df$i, decreasing=FALSE)'
){
  
  nonNA.TF <- !is.na(df$value)
  
  # Transform CII continuous to positive values that can be transformed
  # to probabilities. Did eumiro's answer so distribution will only be translated:
  # https://stackoverflow.com/questions/3931419/turn-a-negative-number-into-a-positive-for-probability
  if( grepl(x=metric, pattern="CII.cont.") ){
    # +1 so probability for minimum value will be 1 (not 0)
    df$value[nonNA.TF] <- df$value[nonNA.TF]-min(df$value, na.rm=TRUE)+1
  }
  
  # Convert df to matrix
  MX <- matrix(data=NA, nrow=tot.bin, ncol=tot.bin)
  
  df <- df[order(df$i, df$j, decreasing=FALSE),]
  MX[ lower.tri(MX, diag=FALSE) ] <- df$value
  
  df <- df[order(df$j, df$i, decreasing=FALSE),]
  MX[ upper.tri(MX, diag=FALSE) ] <- df$value
  
  if( !isSymmetric(MX) ){
    stop(paste0(out.name, ": Matrix not symmetrical."))
  }
  
  # Rows that are not all NAs
  nonNArw.TF <- !apply( X=MX, MARGIN=1, FUN=function(rw) all(is.na(rw)) )
  
  # Maximum value per row (per bin)
  maxrw.v <- rep(NA, times=length(nonNArw.TF))
  maxrw.v[nonNArw.TF] <- apply(X=MX[nonNArw.TF,], MARGIN=1, FUN=max, na.rm=TRUE) 
  
  if( !identical(nonNArw.TF, !is.na(maxrw.v)) ){
    stop(paste0("convertToContactProb(): Checkpoint 1."))
  }
  
  MX <- MX/maxrw.v
  df$value <- MX[ upper.tri(MX, diag=FALSE) ]
  
  rm(nonNArw.TF, MX, tot.bin, maxrw.v)
  
  return(df)
  
}
################################################################################
# rm(list=ls()); gc()
