################################################################################
# Density plot of Cs (raw) per cell/tissue.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/13_Visualise"
    data.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/13_Visualise"
    data.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste("chr", c(1:22, "X"), sep="")
xmax=10
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/multiplot.R"))

densplot <- function( df=df, out.name=paste0(id, "_Cs_densplot") ){
  
  coul <- colorRampPalette( rev( brewer.pal(11, "Spectral") ) )(21)
  
  p <- ggplot(data=df, aes(x=value, colour=variable)) + 
    geom_density() + 
    scale_colour_manual(values=coul) +
    scale_y_continuous(limits=c(0,10), breaks=c(0,5,10)) +
    labs( title=out.name, colour="",
          y=expression( bold( "Density") ), 
          x=expression( bold( "c"["s"] ) ) ) + 
    guides(colour=guide_legend(byrow=TRUE, nrow=1)) +
    theme( legend.text=element_text(size=22, face="bold"),
           legend.title=element_text(size=25, face="bold"),
           legend.position="top",
           aspect.ratio=1) +
    bgr2 
  
  return(p)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(gcb in gcb.v){
  
  p.lst.raw <- list()
  p.lst.scld <- list()
  
  for(chr in chr.v){
    
    id <- paste0(chr, "_", gcb)
    
    # Load PERSIST.MX
    load(paste0(data.dir, "/", chr, "_Persist_", gcb, ".RData"))
    
    # Set to NULL values not needed
    PERSIST.MX[["valsum"]] <- NULL
    PERSIST.MX[["control"]] <- NULL
    
    # Density plot
    p.lst.raw[[chr]] <- densplot(df=melt(PERSIST.MX$hits[,-(1:2)]), 
                                 out.name=paste0(id, "_CsRaw_dp")
    )
    
    p.lst.raw[[paste0(chr,".a")]] <- p.lst.raw[[chr]] +
      scale_x_continuous(breaks=0:xmax) +
      coord_cartesian(xlim=c(0, xmax))
    
  } # chr.v for loop end
  
  p.arr <- ggarrange(plotlist=p.lst.raw, nrow=6, ncol=4, 
                     common.legend=TRUE)
  ggexport(p.arr, width=26, height=40,
           filename=paste0(out.dir, "/", gcb, "_raw_CsDens.pdf"))
  
} # gcb.v for loop end

# rm(list=ls())
