################################################################################
# Make basic violin ggplot
# Adopted from https://www.datanovia.com/en/lessons/ggplot-violin-plot/
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(ggplot2)
# library(Hmisc) # for mean_sdl function
# source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makeViolinPlot <- function(df, x.nme, y.nme, sd.mult, col.nme = NULL, line.cols = NULL, 
                           fill.nme = NULL, fill.cols = NULL, fill.legend = '"none" or factor name', 
                           plot.title, ylim.val = NULL, showOutlier, addmean,
                           geom.viol.scale="count"){
  
  p <- ggplot(data=df, aes_string(x=x.nme, y=y.nme))
  
  violin.note <- "\n geom_violin(trim=T) so tails at range of data"

  # Combine with box plot to add median and quartiles
  # Change fill color by groups, remove legend
  p <- p + 
    geom_violin(aes_string(fill=fill.nme, col=col.nme), 
                lwd=1.5, scale=geom.viol.scale, trim=T) + 
    scale_y_continuous(limits=ylim.val) + 
    scale_fill_manual(values=fill.cols) +
    scale_color_manual(values=line.cols) + 
    labs(title=paste0(plot.title, violin.note)) + 
    guides(fill=fill.legend) + 
    bgr2
  
  if(!showOutlier){ 
    p <- p + geom_boxplot(width=0.1, outlier.shape=NA)
  } else {
    p <- p + geom_boxplot(width=0.1)
  }
    
  # Add mean points +/- SD
  # Use geom = "pointrange" or geom = "crossbar"
  #p <- p + stat_summary(
  #  data=df,
  #  fun.data="mean_sdl", fun.args=list(mult=sd.mult), geom="point", 
  #  color="black", size=2.5
  #)
 
  # Add mean
  if(addmean){
    df.mean <- stack(by(df[[y.nme]], INDICES=df[[x.nme]], FUN=mean, na.rm=F))
    p <- p + 
      geom_point(data=df.mean, aes(x=ind, y=values), size=5, col="black", shape=7)
  }

  return(p)
  
}

# rm(list=ls()); gc()