################################################################################
# Contact probability vs. distance map
# https://github.com/open2c/open2c_examples/blob/master/contacts_vs_distance.ipynb
# https://github.com/deWitLab/GENOVA/blob/master/vignettes/GENOVA_vignette.pdf
### FUNCTION ###################################################################
contactprobVsDistance <- function(df="upper triangle; data.frame(i,j,value,..)",
                                  bin.len='contact map resolution',
                                  scale.diag.ind='Scale values by dividing by
                                  average value of this diagonal index (j-i);
                                  ideally the diagonal index with the highest
                                  expected average value',
                                  n.breaks='Number of x breaks for binning
                                  to smoothen curve'
                                  ){
                                      
  rownames(df) <- NULL
  
  # Remove unwanted contacts
  df <- df[!is.na(df$value),]
  
  if( any(df$i>=df$j) ){
    stop("contactprobVsDistance(): Some not upper triangle contacts.")
  }
  
  # Diagonal index, diag=0 means at diagonal and diag=1 means 1 step from
  # diag=0
  df$diag.ind <- as.factor(as.character(df$j-df$i))
  
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
  df$val.ave.scaled <- df$val.ave/df$val.ave[df$diag.ind==as.character(scale.diag.ind)]
  
  #plot(x=log10(df$diag.bp), y=log10(df$val.ave.scaled), type="l", 
  #     main=paste0("\n contactprobscaledtoAveofdiag=", scale.diag.ind, 
  #                 "_n.breaks=", n.breaks), cex.main=0.5,
  #     xlab="s (in 10^s bp)", ylab="Normalised contact probability P(s)")
  
  #-------------------
  # Smoothen the curve by binning the genomic distance (diagonals combined into bins)
  breaks <- 10^seq(from=log10(min(df$diag.bp)), to=ceiling(log10(max(df$diag.bp))), 
                   length.out=n.breaks)
  cuts <- cut(x=df$diag.bp, breaks=c(0, breaks), include.lowest=T, right=T)
  dfsmooth <- aggregate.data.frame(x=df[,c("diag.bp", "val.ave.scaled")], by=list(cuts), 
                                   FUN=mean, na.rm=T)
  #lines(x=log10(dfsmooth$diag.bp), y=log10(dfsmooth$val.ave.scaled), col="red")
  
  return( list(bydiagonal=df, smooth=dfsmooth) )
  
}
################################################################################
#  rm(list=ls()); gc()



