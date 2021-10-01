################################################################################
# Density plot of correlation coefficients from tissues
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar = T)

# Expands warnings
#options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
COR.dir = paste0(wk.dir, "/out_compare_plot_scatter")
out.dir = paste0(wk.dir, "/out_compare_plot_density")
### OTHER SETTINGS #############################################################
Cs.form = "HiCNormCs"
gcb = "min2Mb"
chr = "chrALL"
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
              "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC"))
type.v = c("kmer", "align") 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(reshape2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.v.len <- length(ct.v)
ct.v <- rep(x=ct.v, times=length(type.v))
type.v <- rep(x=type.v, each=ct.v.len)
i.len <- length(type.v)

COR.df <- sapply(X=1:i.len, simplify=F, FUN=function(i){
  
  load(paste0(COR.dir, "/", chr, "_", gcb, "_", ct.v[[i]], "_", Cs.form, "_", type.v[[i]],
              "_cor.RData"))
  COR <- data.frame(ijset=rownames(COR), COR, type=type.v[[i]], stringsAsFactors=F,
                    row.names=NULL)
  return(COR)
  
})
COR.df <- do.call("rbind", COR.df)
COR.df <- reshape2::melt(data=COR.df, measure.vars=c("pearson", "spearman"),
                         na.rm=F, factorsAsStrings=T)
COR.df$ijset <- factor(x=as.character(COR.df$ijset), levels=c("All", "LR", "SR"))
COR.df$variable <- factor(x=as.character(COR.df$variable), levels=unique(COR.df$variable))

p <- ggplot(data=COR.df, aes(x=value)) + 
  geom_density(aes(fill=ijset, colour=ijset), alpha=0.7) +
  scale_fill_manual(values=c("#F39B7FFF", "#00A087FF", "#8491B4FF")) + 
  scale_colour_manual(values=c("#F39B7FFF", "#00A087FF", "#8491B4FF")) + 
  labs(title=paste0(gcb, "_", chr, "_", Cs.form, "_eachdensityhasdatafromalltissues")) + 
  bgr2 +
  facet_grid(type~variable) +
  theme(strip.background=element_rect(fill="white"))

ggsave(filename=paste0(out.dir, "/", gcb, "_", chr, "_", Cs.form, "_cor_dens.pdf"),
       width=15, height=15, unit="in")

# rm(list=ls()); gc()