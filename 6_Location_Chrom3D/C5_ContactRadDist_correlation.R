################################################################################
# Plot correlation coefficients from C4, combining data from two models
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
model.ids = c("IMR90_LMNB1_GSE49341_hg19", "H1-hESC_LMNB1_hg38")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/6_Location_Chrom3D")
src.dir = out.dir = paste0(wk.dir, "/z_ignore_git/out_ContactRadDist")
### OTHER SETTINGS #############################################################
cor.method = "pearson" #"spearman"
out.name = paste0("chrALL_min2Mb_", paste0(model.ids, collapse="-"), 
                  "_haploid_ContactRadDist_", cor.method)
ntis.corPlot = 1:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
DF <- sapply(model.ids, simplify=F, FUN=function(model.id){
  
  src.model.dir <- paste0(src.dir, "/", model.id)
  src.id <- paste0("chrALL_min2Mb_", model.id, "_haploid_ContactRadDist_", cor.method)
  load(paste0(src.model.dir, "/", src.id, "_PCoeff.RData"))
  df <- cbind(df, model.id=rep(model.id))
  return(df)
  
})

DF <- do.call(rbind, DF)
rownames(DF) <- NULL
DF$ntis <- factor(as.character(DF$ntis), levels=as.character(ntis.corPlot))
coef.col.ind <- grepl(colnames(DF), pattern="Coef", fixed=T)
colnames(DF)[coef.col.ind] <- "Coef"

ggplot(data=DF, aes(x=ntis, y=Coef)) +
  geom_point(size=6, aes(colour=model.id), shape=1, stroke=2, alpha=0.7) +
  scale_y_continuous(limits=c(0,1)) +
  scale_colour_manual(values=c("black", "darkred")) + 
  labs(title=out.name, x=expression("c"["p"]), 
       y=paste0(cor.method, " coefficient")) +
  bgr1 + 
  theme(legend.position="bottom")

ggsave(filename=paste0(out.dir, "/", out.name, "_PCoeff.pdf"),
       height=10, width=10, units="in")

# rm(list=ls()); gc()