################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(ggplot2)
# library(scales)
################################################################################
if(metric=="Cp"){
  
  usedlvl.v <- as.character(1:21)
  coul <- colorRampPalette( rev(RColorBrewer::brewer.pal(11, "Spectral")) )(length(usedlvl.v))
  coul <- setNames(object=coul, nm=usedlvl.v)
  x <- factor(x=as.character(x), levels=usedlvl.v)
  col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=FALSE, na.value="white",
                                    guide=guide_legend(ncol=1), name=bquote(c[p]))
  
# Discretised complementarity
} else if( grepl(x=metric, pattern="CII.disc.", fixed=TRUE) ){ 
  
  x <- factor(x=as.character(x), levels=c("-1", "0", "1"))
  coul <- setNames(object=brewer.pal(n=11, name="RdYlBu")[c(11, 5, 1)], nm=c("-1", "0", "1"))
  coul[["0"]] <- "ivory" 
  col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=FALSE, na.value="white")
  
# Rest of continuous metrics
} else if( grepl(x=metric, pattern="Cs.raw|Cs.norm|CII.cont", fixed=FALSE) ){
  
  bp.stat <- boxplot.stats(x=x, coef=1.5)$stats
  lim.min <- bp.stat[1] # lower whisker
  lim.max <- bp.stat[5] # upper whisker
  
  col.lst[[1]] <- scale_fill_gradientn(colours=col.v, na.value="white",
                                       # scales::squish() for squishing out of bounds values into range
                                       oob=scales::squish, limits=c(lim.min, lim.max))
  
}

# # Test - IGNORE
# library(ggplot2)
# df <- data.frame(
#   x = rep(c(2, 5, 7, 9, 12), 2),
#   y = rep(c(1, 2), each = 5),
#   z = factor(rep(1:5, each = 2)),
#   w = rep(diff(c(0, 4, 6, 8, 10, 14)), 2)
# )
# ggplot(df, aes(x, y)) +
#   geom_raster(aes(fill = z), colour = "grey50")
# 
# df <- expand.grid(x = 0:5, y = 0:5)
# df$z <- factor(as.character(rep(1:9, each = 4)), levels=as.character(9:1)) #runif(nrow(df))
# # default is compatible with geom_tile()
# ggplot(df, aes(x, y, fill = z)) +
#   geom_raster() +
#   scale_y_continuous(trans="reverse")
