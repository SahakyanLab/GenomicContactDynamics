################################################################################
# Plot LENCP.DF data 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(reshape2)
### FUNCTION ###################################################################
plotLENCPDF <- function(DF='LENCP.DF', repfree='T/F', outliers='T/F', col='"cyan"'
                        ){
  
  col.id <- strsplit(x=colnames(DF)[2], split="_", fixed=T)[[1]]
  col.id <- tail(x=col.id, n=1)
  
  col.v = c(full="#FDC776", exon="#e25d6c", intron="#3288BD", intronBYexon="#bb88dd")
  col.v <- c(col.v[col.id], "gray91")
  
  if(repfree==F){
    DF <- DF[ ,!grepl(x=colnames(DF), pattern="repfree") ]
  }
  
  df <- reshape2::melt(data=as.data.frame(DF), id="cp")
  
  #df$value <- df$value/(10^3)
  
  df$cp <- factor(x=as.character(df$cp), 
                  levels=sort(unique(df$cp)))
  df$variable <- factor(x=df$variable, 
                        levels=sort(unique(df$variable), decreasing=F))

  cp.len <- length(levels(df$cp))
  var.len <- length(levels(df$variable))
  
  ylim.v = NULL
  if( grepl(x=colnames(DF)[2], pattern="fr_") ){
    ylim.v = c(0, 0.5)
  } 
  
  boxplot(formula=value~variable*cp, data=df, outline=outliers, ylim=ylim.v,
          xaxt="n", col=col.v[1:var.len], cex.main=1, cex.axis=1.5, boxwex=0.4, ylab="", xlab="",
          main=paste0("colLeftToRightVariable=", 
                      paste(levels(df$variable), collapse="-"))
  )
  
  if(var.len==1){
    at.v <- 1:cp.len
  } else if(var.len==2){
    
    at.v <- seq(1.5, cp.len*var.len + 1L, var.len)
    # Partition each x value by vertical line
    for(i in seq(0.5, cp.len*var.len + 1L, var.len)){ 
      abline(v=i, lty=1, col="gray80")
    }
    
  } else {
    stop("plotLENCPDF: > 2 variables.")
  }
  
  axis(side=1, at=at.v, labels=levels(df$cp), tick=T, cex.axis=1.5)
  
}

# rm(list=ls()); gc()