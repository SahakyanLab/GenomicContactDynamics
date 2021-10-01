################################################################################
# Scaling by "SD" | "MAD" | "MEANSD" | "MEDMAD"
### FUNCTION ###################################################################
scaling <- function(x = 'numeric vector', 
                    approach = '"SD" | "MAD" | "MEANSD" | "TRIMMEANSD" | "MEDMAD" | "BOXWHISK"'
                    ){
  
  if( any(!is.finite(x)) ){
    stop("scaling(): Non-finite values in supplied numeric vector.")
  }
  
  if( approach%in%c("BOXWHISK", "TRIMMEANSD") ){
    
    # Range of boxplot whiskers (default settings, whisker set at 1.5*IQR from Q1 & Q3)
    
    # NAs ignored by default
    bp.stat <- boxplot.stats(x=x, coef=1.5, do.conf=F, do.out=F)$stats
    if( length(bp.stat)!=5 ){
      stop("scaling(): Length of boxplot.stats output unexpected.")
    }
    
  }
  #-------------------
  
  if(approach=="SD"){
    
    scalingFactor <- sd(x=x, na.rm=T)
    x.s <- x / scalingFactor
    print(paste0("scaling(): SD approach applied."), quote=F)
    
  } else if(approach=="MAD"){
    
    scalingFactor <- mad(x=x, center=median(x, na.rm=T), constant=1, low=F, high=F)
    x.s <- x / scalingFactor
    print(paste0("scaling(): MAD approach applied."), quote=F)
    
  } else if(approach=="MEANSD"){
    
    scalingFactor <- c(MEAN=mean(x=x, na.rm=T), SD=sd(x=x, na.rm=T))
    x.s <- (x-scalingFactor["MEAN"]) / scalingFactor["SD"]
    print(paste0("scaling(): MEANSD approach applied."), quote=F)
    
  } else if(approach=="TRIMMEANSD"){
    
    ## Scale with trimmed mean and sd, [min value, Q3+1.5*IQR]
    #x.trim <- x[ x>=min(x, na.rm=T) & x<=bp.stat[5] ]
    # Scale with trimmed mean and sd, [Q1-1.5*IQR, Q3+1.5*IQR]
    x.trim <- x[ x>=bp.stat[1] & x<=bp.stat[5] ]
    scalingFactor <- c(MEAN=mean(x.trim, na.rm=T), SD=sd(x.trim, na.rm=T))
    x.s <- (x-scalingFactor["MEAN"]) / scalingFactor["SD"]
    print(paste0("scaling(): TRIMMEANSD approach applied."), quote=F)
    
  } else if(approach=="MEDMAD"){
    
    scalingFactor <- c(MED=median(x=x, na.rm=T), 
                       MAD=mad(x=x, center=median(x, na.rm=T), 
                               constant=1, low=F, high=F))
    x.s <- (x-scalingFactor["MED"]) / scalingFactor["MAD"]
    print(paste0("scaling(): MEDMAD approach applied."), quote=F)
    
  } else if(approach=="BOXWHISK"){
    
    scalingFactor <- abs(bp.stat[5]-bp.stat[1])
    x.s <- x / scalingFactor
    
    print(paste0("scaling(): BOXWHISK approach applied."), quote=F)
    
  } else {
    stop(paste0("scaling(): Invalid approach argument."))
  }
  
  return( list(scaled=x.s, scalingFactor=scalingFactor) )
  
}
################################################################################
