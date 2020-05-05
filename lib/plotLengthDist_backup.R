################################################################################
# Make density plot of length (or any vector of values) distribution
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(ggplot2)
# source(paste0(lib, "GG_bgr.R"))
################################################################################
plotLengthDist <- function(vec = "",
                           vline.v = c(0.5e6, 2e6),
                           col = "deepskyblue3",
                           out.name = "",
                           out.dir = "/dir",
                           label.x = bquote(bold("log"["10"]~"("~"L"~")"))
                           ){
  
  p <- ggplot(data=data.frame(values=vec), aes(x=values)) +
    geom_density(fill=col, colour=col) +
    # aes(y=..scaled..)
    # aes(y=..count..)
    # aes(values, stat(count))
    labs(title=out.name,
         x=label.x, 
         y=expression(bold("Density"))
    ) +
    bgr2 
  
  if( !is.null(vline.v) ){
    p <- p + geom_vline( data=data.frame(vline.v), linetype="dashed",
                         colour="black", size=0.7, aes(xintercept=vline.v))
  }
  
  ggsave(filename=paste0(out.dir, "/", out.name, ".pdf"),
         units="in", width=10, height=10, plot=p)
  
  return(p)
  
}
