################################################################################
# Contact probability vs. distance map
# https://github.com/open2c/open2c_examples/blob/master/contacts_vs_distance.ipynb
# https://github.com/deWitLab/GENOVA/blob/master/vignettes/GENOVA_vignette.pdf
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# lib = "/Users/ltamon/DPhil/lib"
# source(paste0(lib, "/splitNumericVector.R"))
### FUNCTION ###################################################################
contactprobVsDistance <- function(df = "upper triangle; data.frame(i,j,value,..)",
                                  bin.len = 'contact map resolution',
                                  scale.diag.ind = 'Scale values by dividing by
                                  average value of this diagonal index (j-i);
                                  ideally the diagonal index with the highest
                                  expected average value. If NA, no scaling',
                                  
                                  # Smoothening of curve
                                  smooth.method = "BIN means discrete windows;
                                  SLIDE means overlapping windows; NONE means
                                  no smoothening",
                                  
                                  # smooth.method=BIN arguments
                                  n.breaks='Number of x breaks for binning
                                  to smoothen curve',
                                  
                                  # smooth.method=SLIDE arguments
                                  d='radius of sliding window; ideally a multiple
                                  of bin.len',
                                  numSLwind='number of sliding windows between
                                  min and max diag.bp; set d and numSLwind such that
                                  windows are overlapping'
                                  ){
                                      
  rownames(df) <- NULL

  # Remove unwanted contacts
  df <- df[!is.na(df$value),]
  
  #if( any(df$i>=df$j) ){
  #  stop("contactprobVsDistance(): Some not upper triangle contacts.")
  #}
  
  # Diagonal index, diag=0 means at diagonal and diag=1 means 1 step from
  # diag=0
  df$diag.ind <- as.factor(as.character( abs(df$j-df$i) ))
  
  #-------------------
  # Calculate number of valid contacts and sum of values per diagonal index
  df <- aggregate.data.frame(x=df$value, by=list(df$diag.ind), FUN=function(x){
    c(Nij=sum(!is.na(x)), val.sum=sum(x, na.rm=T)) 
  })
  df <- cbind.data.frame(diag.ind=df$Group.1, 
                         diag.bp=bin.len*as.numeric(as.character(df$Group.1)), 
                         Nij=df$x[,"Nij"], 
                         val.sum=df$x[,"val.sum"], 
                         val.ave=df$x[,"val.sum"]/df$x[,"Nij"])
  df <- df[order(df$diag.bp),]
  
  # Scale plot to average value of diagonal index = scale.diag.ind
  if( !is.na(scale.diag.ind) ){
    
    df$val.ave.plot <- df$val.ave/df$val.ave[df$diag.ind==as.character(scale.diag.ind)]
    print("contactprobVsDistance(): Contact values scaled.")
    
  } else {
    
    df$val.ave.plot <- df$val.ave
    print("contactprobVsDistance(): Contact values not scaled.")
    
  }
  
  #plot(x=log10(df$diag.bp), y=log10(df$val.ave.plot), type="l", 
  #     main=paste0("\n contactprobscaledtoAveofdiag=", scale.diag.ind, 
  #                 "_n.breaks=", n.breaks), cex.main=0.5,
  #     xlab="s (in 10^s bp)", ylab="Normalised contact probability P(s)")
  
  #-------------------
  
  # Smoothen the curve by binning the genomic distance (diagonals combined into bins)
  
  if(smooth.method=="BIN"){
    
    # Discrete windows/bins
    breaks <- 10^seq(from=log10(min(df$diag.bp)), to=ceiling(log10(max(df$diag.bp))), 
                     length.out=n.breaks)
    cuts <- cut(x=df$diag.bp, breaks=c(0, breaks), include.lowest=T, right=T)
    dfsmooth <- aggregate.data.frame(x=df[,c("diag.bp", "val.ave.plot")], by=list(cuts), 
                                     FUN=mean, na.rm=T)
    
  } else if(smooth.method=="SLIDE"){
    
    # Overlapping windows/bins
    cuts <- splitNumericVector(x=df$diag.bp, d=d, action="SLIDE", numSLwind=numSLwind, 
                               boundsSL=c(min(df$diag.bp), max(df$diag.bp))
                               )
    dfsmooth <- sapply(X=1:length(cuts), simplify=F, FUN=function(i){
      ind <- cuts[[i]]
      c(diag.bp=mean(df[ind,"diag.bp"], na.rm=T), 
        val.ave.plot=mean(df[ind,"val.ave.plot"], na.rm=T))
    })
    dfsmooth <- as.data.frame(do.call("rbind", dfsmooth))
    
  } else if(smooth.method=="NONE"){
    dfsmooth <- NULL
  } else {
    stop("contactprobVsDistance(): Invalid smooth.method argument.")
  }
  
  #lines(x=log10(dfsmooth$diag.bp), y=log10(dfsmooth$val.ave.plot), col="red")
  
  return( list(bydiagonal=df, smooth=dfsmooth) )
  
}
################################################################################
#  rm(list=ls()); gc()



