################################################################################
# Density plot of log10 distances (in bp) from reference points
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/5_Location_2D")
src.dir = paste0(wk.dir, "/out_distToRefPoints")
out.dir = paste0(wk.dir, "/out_distToRefPoints_plot")
### OTHER SETTINGS #############################################################
Cp.v = 1:21
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
chr.id = "chrALL" #"chr21and22"
options(digits=15) # Digits beyond 15 could be unreliable
# Should match ref.points order in src data (IJDISTTOREF.MX)
ref.points.name = c("chr.start.bp", "chr.end.bp", "centromere.midP.bp") 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(RColorBrewer)
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
#################################################################################
ref.points.id <- paste(ref.points.name, collapse="_")
df <- sapply(X=chr.v, simplify=F, FUN=function(chr){
  load(paste0(src.dir, "/", gcb, "_", chr, "_ijdisttoref_", ref.points.id, ".RData"))
  print(paste0(chr, ": Data obtained."))
  return(IJDISTTOREF.MX)
})
df <- as.data.frame(do.call("rbind", df))

df$MINIMUMof.chr.start.bpVSchr.end.bp <- pmin(abs(df$chr.start.bp), abs(df$chr.end.bp))
df$centromere.midP.bp <- abs(df$centromere.midP.bp)
df$chr.start.bp <- df$chr.end.bp <- NULL

df <- reshape2::melt(df, id="Cp")
df$Cp <- factor(as.character(df$Cp), levels=as.character(Cp.v))

coul <- colorRampPalette( rev( brewer.pal(11, "Spectral") ) )( length(levels(df$Cp)))

# Density plot
p <- ggplot(data=df, aes(x=log10(value), group=Cp)) +
  geom_density( position="identity", aes(colour=Cp) ) + 
  scale_colour_manual(values=coul) + 
  labs(title=paste0(gcb, "_", chr.id, "_AbsoluteDistFromRef_lastbincontactsExcluded"), 
       x="log10(distinbp)", colour="Cp") +
  guides(colour=guide_legend(ncol=1)) +
  bgr2 + 
  facet_grid(. ~ variable)
  

ggsave(filename=paste0(out.dir, "/", gcb, "_", chr.id, "_ijdisttoref_", ref.points.id, "_density.pdf"),
       units="in", width=10, height=10, plot=p)


# rm(list=ls()); gc()