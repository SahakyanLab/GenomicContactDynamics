################################################################################
# Heatmap of radial distances of beads on each chromosome; useful for comparing
# models across cell lines
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    repeat.dir = "/Users/ltamon/Database/ucsc_tables/hsa_RepeatMasker"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    repeat.dir  = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_RepeatMasker"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
domXYZR.dir = paste0(wk.dir, "/out_AddXYZR")
out.dir = paste0(wk.dir, "/out_chrPos")
### OTHER SETTINGS #############################################################
out.name = paste0(model.id, "_haploid")
# Chromosome in the model
chr.v = paste("chr", c(1:22, "X"), sep="")
bin.len = 40000 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape)
library(ggplot2)
library(gplots)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Load DOMXYZR.DF
load(file=paste0(domXYZR.dir, "/", out.name, "_domXYZR.RData"))

DOMXYZR.DF <- DOMXYZR.DF[,c("radDist", "chr", "start", "end")]
# Because end position of a domain is same as the start coordinate
# of the previous domain 
DOMXYZR.DF$start <- DOMXYZR.DF$start+1L

# Length of largest chromosome
# 2.49e8 ~ 2.5e8
format(as.numeric(max(DOMXYZR.DF$end)), digits=3, scientific=TRUE)

round(max(DOMXYZR.DF$end)/1e6)

pos <- seq(from=1, to=2.5e8, by=bin.len)
pos[pos==0] <- 1

# Data for heatmap
HM.MX <- matrix(data=NA, ncol=length(pos), nrow=length(chr.v),
                dimnames=list(chr.v, pos))

rw.len <- nrow(DOMXYZR.DF)
for(i in 1:rw.len){
  rw <- DOMXYZR.DF[i,]
  HM.MX[chr.v==rw$chr, (pos>=rw$start & pos<=rw$end) ] <- rw$radDist
}
dimnames(HM.MX)[[2]] <- pos/1e6

df <- melt(HM.MX)
colnames(df) <- c("chr", "pos", "r")
df$chr <- factor(as.factor(df$chr), levels=chr.v)

# Heatmap 
ggplot(df, aes(x=pos, y=chr, fill=r)) + 
  geom_tile()+
  scale_fill_gradientn(colours=colorRampPalette( rev( brewer.pal(11, "Spectral") ) )(n = 11),
                       na.value="white", limits=c(0,6)) + 
  labs(title=paste0(out.name, "_bin", bin.len),
       x=expression(bold("Chromosome position (Mb)")),
       y=NULL)+
  bgr1 + 
  theme(legend.text=element_text(size=10))
  
ggsave(filename=paste0(out.dir, "/", out.name, "_bin", bin.len, "_chrPos_plot.pdf"), 
       units="in", width=10, height=10)

# rm(list=ls())
