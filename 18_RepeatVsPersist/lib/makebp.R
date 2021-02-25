############################################################################### 
# Make boxplot
### FUNCTION ###################################################################
makebp <- function(df=x, x="ind", y=calc, xlab=paste0(wm, "-bin"), ylab=calc, 
                   addjitter=addjitter, plot.id=""
){
  
  # By default, ignore missing values in either the response or the group.
  eval(parse(text=paste0(
    'boxplot(', y, '~', x, ', outline=FALSE, data=df, xlab=xlab, 
    ylab=ylab, boxwex=0.6, cex.axis=1.2, col="#FDC776", cex.main=0.1,
    main=plot.id)'
  )))
  
  if(addjitter){
    # Add data points
    levelprop.v <- summary(df[[x]])/nrow(df)
    for( i in 1:length(levels(df[[x]])) ){
      # Take the x-axis indices and add a jitter, proportional to the N in each level
      jitt <- jitter(rep(i, length(df[[y]])), amount=levelprop.v[i]/2)
      points(jitt, df[[y]], cex=1, col=adjustcolor("black", alpha.f=0.01), pch=16) 
      rm(jitt)
    }
  }
  
}