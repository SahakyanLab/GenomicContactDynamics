################################################################################
# Make density plot of length (or any vector of values) distribution
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(ggplot2)
# source(paste0(lib, "GG_bgr.R"))
################################################################################
plotLengthDist <- function(df = data.frame(variable=NA, value=NA),
                           vline.v = c(0.5e6, 2e6),
                           col.v = NULL,
                           out.name = "",
                           out.dir = "/dir",
                           label.x = bquote(bold("log"["10"]~"("~"L"~")")),
                           addlabs = TRUE,
                           addlegend= TRUE
                           ){
  colnames(df) <- c("variable", "value")
  p <- ggplot(data=df, aes_string(x=colnames(df)[2])) +
    geom_density(aes_string(fill=colnames(df)[1], col=colnames(df)[1])) + #alpha=0.25) +
    # aes(y=..scaled..)
    # aes(y=..count..)
    # aes(values, stat(count))
    bgr2 
  
  if(addlabs){
    p <- p + labs(title=out.name,
                  x=label.x, 
                  y=expression(bold("Density"))
    ) 
  } else {
    p <- p + labs(title=NULL, x=NULL, y=NULL, fill=NULL)
  }
  
  if(addlegend==FALSE){
    p <- p + theme(legend.position="none")
  }
  
  if( !is.null(vline.v) ){
    p <- p + geom_vline( data=data.frame(vline.v), linetype="dashed",
                         colour="black", size=2, aes(xintercept=vline.v))
  }
  if( !is.null(col.v) ){
    p <- p + 
      scale_fill_manual(values=col.v) + 
      scale_color_manual(values=col.v) 
  }
  ggsave(filename=paste0(out.dir, "/", out.name, ".pdf"),
         units="in", width=10, height=10, plot=p)
  return(p)
}
