############################################################################### 
# Make boxplot
### FUNCTION ###################################################################
makebp <- function(df=x, calc, xlab, ylab, addjitter, plot.id){
  
  ind.v <- as.character(levels(x$ind))
 
  if( !is.factor(df$ind) ){
    stop("Boxplot making: x$ind not a factor as expected.")
  }
  if( !is.numeric(df$mutbin) ){
    stop("Boxplot making: x$mutbin not numeric as expected.")
  }
  if( !is.numeric(df[[calc]]) ){
    stop(paste0("Boxplot making: x$", calc, " not numeric as expected."))
  }
  
  df$ind <- as.character(df$ind)
  bin.TF <- is.finite(df[[calc]])
  
  binPerCp <- binmutPerCp <- rep(x=0, times=length(ind.v))
  names(binPerCp) <- names(binmutPerCp) <- ind.v
  
  tmp <- table(df$ind[bin.TF])
  binPerCp[names(tmp)] <- tmp
  
  tmp <- table(df$ind[df$mutbin==1 & bin.TF])
  binmutPerCp[names(tmp)] <- tmp
  rm(tmp)
  
  if( !identical(names(binPerCp), names(binmutPerCp)) ){
    stop("Boxplot making: names of binPerCp and binmutPerCp not identical.")
  }
  
  percbinmut <- binmutPerCp/binPerCp*100
  percbinmut[!is.finite(percbinmut)] <- 0
  
  binPerCp.str <- paste0("\n binPerCp1to21_", paste(x=binPerCp, collapse="_"))
  binmutPerCp <- paste0("\n binmutPerCp1to21_", paste(x=binmutPerCp, collapse="_"))
  percbinmut.str <- paste0("\n %binmutPerCp1to21_",
                           paste(x=round(percbinmut, digits=3), collapse="_"))
  
  plot.id <- paste0(plot.id, binPerCp.str, binmutPerCp, percbinmut.str)
  rm(bin.TF, binPerCp.str, binmutPerCp, percbinmut.str)
  
  if( !is.character(df$ind) ){
    stop("Boxplot making: x$ind for boxplot not a character as expected.")
  }
  df$ind <- factor(x=df$ind, levels=ind.v)
  
  if( !identical(as.character(1:21), levels(df$ind)) ){
    stop("Boxplot making: wrong levels of x$ind for boxplot")
  }
  # By default, ignore missing values in either the response or the group.
  eval(parse(text=paste0(
    'boxplot(', calc, '~ind,', ' outline=FALSE, data=df, 
    xlab=xlab, ylab=ylab, boxwex=0.6, cex.axis=1.2, col="#FDC776", cex.main=0.2,
    main=plot.id)'
  )))
  
  if(addjitter){
    # Add data points
    levelprop.v <- summary(df$ind)/nrow(df)
    for( i in 1:length(levels(df$ind)) ){
      # Take the x-axis indices and add a jitter, proportional to the N in each level
      jitt <- jitter(rep(i, length(df[[calc]])), amount=levelprop.v[i]/2)
      points(jitt, df[[calc]], cex=1, col=adjustcolor("black", alpha.f=0.01), pch=16) 
      rm(jitt)
    }
  }
  
  return(
    list(percbinmut=percbinmut,  binPerCp=binPerCp)
  )
  
}