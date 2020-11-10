################################################################################
# Plot bincount per position
# deva, R/3.6.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/8_FeatureVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/10_ChromatinFeatures"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = out.dir = paste0(wk.dir, "/out_bincount")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
celltiss.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
               "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC",
               "hg19")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(reshape)
library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(celltiss in celltiss.v){
  
  # Load bincount list
  load(file=paste0(data.dir, "/chrALL_", gcb, "_", celltiss, "_bincountplot.RData"))
  
  # Comparing bincount.exp vs. bincount.act
  bincount.act <- 100*(bincount$exp-bincount$act)/bincount$exp
  bincount.act <- melt.array(bincount.act)
  maxval <- format(max(bincount.act$value), digits=4)
  colnames(bincount.act) <- c("cp", "pos", "value")
  exp <- paste( paste(names(bincount$exp), bincount$exp, sep="-"), collapse="," )
  
  # Plot
  cp.v.len <- length(unique(bincount.act$cp))
  coul <- colorRampPalette( rev( brewer.pal(11, "Spectral") ) )(cp.v.len)
  
  min.x <- min(bincount.act$pos)
  max.x <- max(bincount.act$pos)
  p <- ggplot(data=bincount.act, aes(x=pos, y=value)) +
    geom_point( size=4, aes(colour=factor(cp)) ) +
    scale_x_continuous(labels=as.vector(rbind("",seq(from=min.x, to=max.x, by=2)))[-1],
                       breaks=min.x:max.x,
                       limits=c(min.x-0.5, max.x+0.5)
    ) + 
    scale_colour_manual(values=coul) + 
    labs(title=paste0( "chrALL_", gcb, "_", celltiss, "_Npos=", 
                       unique(table(bincount.act$cp)), "_max=", maxval, "%\n", exp),
         x=expression(bold("Pos")), 
         y=expression(bold("% Bins missing")),
         colour=expression(bold("c"["p"]))
    ) + 
    guides(colour=guide_legend(ncol=1)) + 
    bgr2 +
    theme(plot.title=element_text(size=5),
          legend.text=element_text(size=20, face="bold"),
          legend.title=element_text(size=25, face="bold")
           #,
           #panel.grid.major.x=element_line(colour="gray90")
           #,
           #panel.grid.minor.y=element_line(colour="gray80")
           ) 
  ggsave(filename=paste0(out.dir, "/chrALL_", gcb, "_", celltiss, 
                         "_bincountplot.pdf"), 
         units="in", width=10, height=10, plot=p)
  
  rm(bincount.act, bincount, exp, maxval, coul, cp.v.len); gc()
  print(paste0(celltiss, " done!"), quote=FALSE)
  
} # celltiss.v for loop end

# rm(list=ls()); gc()