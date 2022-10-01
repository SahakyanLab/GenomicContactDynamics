################################################################################
# Correlation of contact gap ( j-i not j-i-1 to avoid log10(0) ) and complementarity
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
    home.dir = "/project/sahakyanlab/ltamon" 
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/25_LengthDependence")
complperChr.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_constraints_GfreeSingleNorm/merged_final")
complallChr.dir = complperChr.dir
out.dir = paste0(wk.dir, "/out_gap_complementarity_COR")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c("ALL", 1:22, "X"))
compl.type = "kmer" # kmer | Gfree | sdDifference | align
bins.val = 30
cuts.val = 12
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggpubr)
library(hexbin)
library(Hmisc)
library(reshape2)
library(viridis)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeHexbinggplot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( compl.type %in% c("Gfree", "sdDifference") ){
  obj.id <- "kmer"; val.id <- compl.type
} else if( compl.type %in% c("kmer", "align") ){
  obj.id <- compl.type; val.id <- "C||"
} 

chr.v.len <- length(chr.v)
p.lst <- vector("list", length=chr.v.len)
names(p.lst) <- chr.v
COR.MX <- matrix(NA, nrow=chr.v.len, ncol=2, dimnames=list(chr.v, c("pearson", "spearman")))
out.id <- paste0(gcb, "_", compl.type, "_bins", bins.val, "_cuts", cuts.val)

for(chr in chr.v){
  
  COMPL.dir <- complperChr.dir
  if(chr == "chrALL"){ COMPL.dir <- complallChr.dir }
  
  load(paste0(COMPL.dir, "/", chr, "_", obj.id, "_", gcb, ".RData"))
  rownames(CII.MX) <- NULL
  
  is_nonNA_compl <- !is.na(CII.MX[,val.id])
  gaps <- CII.MX[,"j"] - CII.MX[,"i"]

  p.lst[[chr]] <- makeHexbinggplot(xvar=log10(gaps[is_nonNA_compl]), 
                                   yvar=CII.MX[is_nonNA_compl,val.id], 
                                   bins=bins.val, 
                                   cuts=cuts.val,
                                   xlab="contact gap (j-i-1)",
                                   ylab=paste0("c||_", compl.type),
                                   title=paste0(chr, "_", out.id),
                                   col=viridis(cuts.val))$hexplot
  
  
  chr.ind <- which(chr.v == chr)
  COR.MX[chr.ind, "pearson"] <- cor(x=gaps[is_nonNA_compl], y=CII.MX[is_nonNA_compl,val.id], method="pearson")
  COR.MX[chr.ind, "spearman"] <- cor(x=gaps[is_nonNA_compl], y=CII.MX[is_nonNA_compl,val.id], method="spearman")
  
}

p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=6)
ggexport(p.arr, height=40, width=60, units="in",
         filename=paste0(out.dir, "/", out.id, "_gap_complementarity_correlation.pdf"))
save(COR.MX, file=paste0(out.dir, "/", gcb, "_", compl.type, "_gap_complementarity_correlation.RData"))

p <- ggplot(data=melt(COR.MX), aes(value)) + 
  geom_density(aes(col=Var2, fill=Var2), alpha=0.5) + 
  labs(title=paste0("chrALL_", out.id, "_allValsInCORMX"), x="cor coefficient") + 
  bgr2
ggsave(filename=paste0(out.dir, "/", gcb, "_", compl.type, "_gap_complementarity_correlation_allValsInCORMX.pdf"),
       width=10, height=10, unit="in")

# rm(list=ls()); gc()