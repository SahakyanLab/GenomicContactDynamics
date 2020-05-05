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
  list = "DNAseq",
  title = "myseqlogo",
  ncol = 2L
){
  
  p <- ggplot() + 
    geom_logo(list, method="prob") +
    labs(title=title,
         x=expression(bold("Position"))
    ) +
    bgr2 + 
    theme(strip.text.x=element_text(size=20, face="bold", vjust=2.75),
          strip.background = element_blank(),
          aspect.ratio=1) + 
    facet_wrap(~seq_group, ncol=ncol, scales='free_x')
  
  return(p)
  
}
