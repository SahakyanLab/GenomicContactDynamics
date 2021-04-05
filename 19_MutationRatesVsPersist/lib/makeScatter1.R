############################################################################### 
# Specific function to do Mann-Whitney test using x and df object.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(ggsci)
#library(ggplot2)
#source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makeScatter <- function(df, calc, yint=NULL, plot.id){
  
  if( !is.character(df$ind) ){
    stop("Scatterplot making: ind not a character as expected.")
  }
  df$ind <- factor(as.character(df$ind), 
                   levels=c(as.character(sort(as.numeric(unique(df$ind)))))
                   )
  
  df$pval <- c(`0`="ns", `1`="sig")[ as.character(as.numeric(df$pval<0.05)) ]
  df$pval <- factor(df$pval, levels=c("ns", "sig"))
  alpha.v <- c(0.3, 1)
  names(alpha.v) <- levels(df$pval)
  
  shape.v <- c(15, 16, 17, 18, 19)
  names(shape.v) <- levels(df$sigEpLim.id)
  
  p <- ggplot(data=df, aes(x=ind, y=values)) +
    geom_point(aes(colour=SIG.id, alpha=pval, shape=sigEpLim.id), size=3) +
    scale_alpha_manual(values=alpha.v) + 
    scale_shape_manual(values=shape.v) + 
    labs(x="Cp", y=calc, title=plot.id)
  + 
    bgr2 + 
    theme(plot.title=element_text(size=15), legend.text=element_text(size=10),
          legend.title=element_text(size=10))
    
  if( !is.null(yint) ){
    p <- p + geom_hline(yintercept=yint, linetype="dashed", colour="gray70", size=0.5) 
  }
    
  return(p)
  
}
