################################################################################
# Normalise contact value (i.e. complementarity values) based on distance between 
# contact bins (the farther the bins, that harder it should be to interact).
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
scaleContactByDist <- function(df.lst = "list of two map values to be compared",
                               plot.file = "complete path of plot"){
  
  subj.TF <- grepl(x=names(df.lst), pattern="CII.cont", fixed=T) 
  
  if( sum(subj.TF)!=1 ){
    
    print("scaleContactByDist(): No map needing scaling.")
    return(df.lst)
    
  } else {
    
    subj <- names(df.lst)[subj.TF]
    ref <- names(df.lst)[!subj.TF]
    
  }
  
  # Get multipliers from ref map
  
  val.TF <- !is.na(df.lst[[ref]]$value) & df.lst[[ref]]$value > 0
    
  ref.gapbin <- abs( df.lst[[ref]]$j-df.lst[[ref]]$i )[val.TF]
  ref.gapbin <- table(ref.gapbin)
  
  #plot(x=as.numeric(names(ref.gapbin)), y=as.numeric(ref.gapbin),
  #     cex=0.1)
  
  min.val <- min(ref.gapbin)
  max.val <- max(ref.gapbin)
  #ref.gapbin <- (ref.gapbin - min.val) / (max.val-min.val)
  ref.gapbin <- ref.gapbin / max.val
  ref.gapbin <- cbind(gapjMINUSi=as.numeric(names(ref.gapbin)),
                      multiplier=ref.gapbin)
  dimnames(ref.gapbin)[[1]] <- ref.gapbin[,"gapjMINUSi"]
  
  # Transform CII continuous to positive values that can be transformed
  # to probabilities. Did eumiro's answer so distribution will only be translated:
  # https://stackoverflow.com/questions/3931419/turn-a-negative-number-into-a-positive-for-probability
  
  nonNA.TF <- !is.na(df.lst[[subj]]$value)

  # +1 so probability for minimum value will be 1 (not 0)
  df.lst[[subj]]$value[nonNA.TF] <- df.lst[[subj]]$value[nonNA.TF] - 
                                    min(df.lst[[subj]]$value, na.rm=T) + 1
  
  # Scale subj map
  
  subj.gapbin.v <- abs( df.lst[[subj]]$j-df.lst[[subj]]$i )
  
  withMult.TF <- subj.gapbin.v%in%ref.gapbin[,"gapjMINUSi"]
  subj.mult.v <- rep(x=1, times=length(subj.gapbin.v))
  subj.mult.v[withMult.TF] <- ref.gapbin[ as.character(subj.gapbin.v[withMult.TF]), "multiplier" ]
  
  df.lst[[subj]]$value <- df.lst[[subj]]$value * subj.mult.v
  
  return(df.lst)
  
}

# rm(list=ls()); gc()