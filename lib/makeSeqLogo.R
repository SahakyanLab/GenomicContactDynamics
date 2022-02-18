################################################################################
# Function to make DNA sequence logos out of a list of character vectors
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(ggseqlogo) # geom_logo
# library(ggplot2)
# source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makeSeqLogo <- function(
  method = "prob", # "bit"
  list = "DNAseq",
  title = "myseqlogo",
  ncol = 2L
){
  
  
  if(method=="prob"){
    lim.v <- c(0,1)
  } else if(method=="bit"){
    lim.v <- c(0,2)
  } else {
    stop("Invalid method argument.")
  }
  
  p <- ggplot() + 
    geom_logo(list, method=method) +
    scale_y_continuous(limits=lim.v) + 
    labs(title=title,
         x=expression(bold("Position"))
    ) +
    bgr2 + 
    theme(strip.text.x=element_text(size=20, face="bold", vjust=0.5),
          strip.background = element_blank(),
          aspect.ratio=1) + 
    facet_wrap(~seq_group, ncol=ncol, scales='free_x') 
  
  return(p)
  
}
