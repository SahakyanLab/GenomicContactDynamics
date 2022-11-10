################################################################################
# Analysis plots comparing complementarity with various values/features
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
valmx.dir = paste0(wk.dir, "/out_complementarityVsContactValue") 
out.dir = paste0(wk.dir, "/out_complementarityVsContactValue_plot") 
state.group.file = paste0(wk.dir, "/Phylo-HMRF_state_group.csv")
group.col.file = paste0(wk.dir, "/Phylo-HMRF_state_group_col.csv")
### OTHER SETTINGS #############################################################
value.id = "genome_state_Phylo-HMRF_mapping_contact50K_norm"
state.group = "Grp3"
gcb = "min2Mb"
type = "kmer"
compl.type = "CII"
chr.v = paste0("chr", c(17, 22))
bin.len = 50000
gap.range.bins.closed = NULL # j - i - 1, No gap filtering if NULL
out.id = paste0("hg38.108_bin", bin.len, "bp_", value.id, "_gaprangebins_", 
                gap.range.bins.closed[1], "_", gap.range.bins.closed[2],
                "_state", state.group)  

coul.spec <- c(`GSM3685694_panTro5_norm`="cyan2", `GSM3685695_panPan2_norm`="blue",
           `GSM3685696_gorGor4_norm`="darkgreen") # hotpink for human
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(vioplot)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
df <- sapply(chr.v, simplify=F, FUN=function(chr){
  
  load(paste0(valmx.dir, "/", chr, "_", type, "_", gcb, "_", value.id, ".RData"))
  VAL.MX$Cp <- NULL
  
  if( !is.null(gap.range.bins.closed) ){
    print(paste0(chr, ": Gap Filtering"))
    gaps <- unname(VAL.MX$j - VAL.MX$i - 1)
    VAL.MX <- VAL.MX[gaps >= gap.range.bins.closed[1] & gaps <= gap.range.bins.closed[2],]
  }
  
  print(paste0(chr, ": Data obtained."))
  return(VAL.MX)
  
})
df <- do.call("rbind.data.frame", df)

#
grp.df <- read.csv(state.group.file, header=T, stringsAsFactors=F)
is.nonNAstate <- !is.na(df$PHYLOHMRFstate)
df$grp[is.nonNAstate] <- grp.df[[state.group]][ as.numeric(as.character(df$PHYLOHMRFstate))[is.nonNAstate] ]

# Reorder states - increasing median complementarity measure and order as levels
median.df <- aggregate(x=df$`C||`, by=list(as.character(df$grp)), FUN=median, na.rm=T)
median.df <- median.df[order(median.df$x, na.last=F), ]

df$grp <- factor(as.character(df$grp), levels=as.character(median.df$Group.1))
median.df$Group.1 <- factor(as.character(median.df$Group.1), levels=levels(df$grp))
colnames(median.df) <- c("grp", "value")

#
data.table::setnames(x=df, old="C||", new="CII")

col.df <- read.csv(group.col.file, header=T, stringsAsFactors=F)
col.df$Col <- paste0("#", col.df$Col)
col.df <- col.df[col.df$Grp %in% levels(df$grp),]
coul.grp <- setNames(col.df$Col, nm=col.df$Grp)

# add horizontal line
p.cii <- ggplot(data=df[!is.na(df$grp),], aes_string(x="grp", y=compl.type)) +
  geom_violin(scale="width", aes(fill=grp), col="white", alpha=0.5, trim=T) +
  geom_point(data=median.df, aes(x=grp, y=value), alpha=1, size=4, pch=3, stroke=1.1) + 
  stat_boxplot(geom="errorbar", width=1/length(levels(df$grp))) + 
  scale_fill_manual(values=unname(coul.grp[levels(df$grp)])) + 
  scale_colour_manual(values=unname(coul.grp[levels(df$grp)])) + 
  bgr2 + 
  theme(legend.position="none", 
        axis.text.x=element_text(size=7, angle=90),
        aspect.ratio=0.25) 
ggsave(filename=paste0(out.dir, "/", out.id, "_violin.pdf"), height=7, width=15)

norm.df <- reshape2::melt(df[!is.na(df$grp),-(1:6)])
colnames(norm.df) <- c("grp", "species", "value")
norm.df$species <- factor(norm.df$species, levels=unique(norm.df$species))

p.cf <- ggplot(data=norm.df, aes(x=value)) +
  geom_density(aes(y= ..scaled.., col=species)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", -(x))) + 
  scale_y_continuous(breaks=c(0,1)) + 
  scale_color_manual(values=coul.spec[levels(norm.df$species)]) + 
  bgr2 +
  theme(legend.position="bottom", legend.text=element_text(size=2), 
        axis.text.x=element_text(size=7, angle=90),
        aspect.ratio=1.5) +
  facet_grid(. ~ grp) +
  coord_flip() 
  
p.arr <- ggarrange(plotlist=list(a=p.cii, b=p.cf), nrow=2, ncol=1, legend=NULL)
ggexport(p.arr, width=15, height=10, units="in",
         filename=paste0(out.dir, "/", out.id, "_violin.pdf")) 

# rm(list=ls()); gc()


# med1.df <- aggregate(x=df$`GSM3685694_panTro5_norm`, by=list(df$grp), FUN=median, na.rm=T)
# med2.df <- aggregate(x=df$`GSM3685695_panPan2_norm`, by=list(df$grp), FUN=median, na.rm=T)
# med3.df <- aggregate(x=df$`GSM3685696_gorGor4_norm`, by=list(df$grp), FUN=median, na.rm=T)
# 
# #coul <- c("#313695", brewer.pal(n=11, name="RdYlBu")[c(9, 6, 4, 1)])
# #coul <- setNames(colorRampPalette(coul)(length(levels(df$grp))),
# #                 nm=as.character(levels(df$grp)))
# 
# pdf(file=paste0(out.dir, "/", out.id, "_violin.pdf"), height=10, width=15)
# 
# par(mfrow=c(3,1))
# 
# boxplot(formula=`C||`~grp, data=df, outline=T, main=out.id , col=coul[levels(df$grp)], cex.main=0.8)
# legend("topright", legend=names(coul[ median.df$Group.1[!is.na(median.df$x)] ]),
#        fill=coul[ median.df$Group.1[!is.na(median.df$x)] ], cex=0.5, bty="n", 
#        x=0, y=max(df$`C||`[!is.na(df$grp)], na.rm=T), horiz=F, ncol=6)
# 
# vioplot(formula=`C||`~grp, data=df, las = 1, plotCentre="line", cex.main=0.8, drop=T,
#         main=paste0(out.id, "_horizontalLineIsMedian_missingValueCategoryNoDataPoint"),
#         col=coul[ median.df$Group.1[!is.na(median.df$x)] ])
# #vioplot(formula=GSM3685694_panTro5_norm~grp, data=df, las = 1, plotCentre="line", cex.main=0.8, drop=T,
# #        main=paste0(out.id, "_horizontalLineIsMedian_missingValueCategoryNoDataPoint"),
# #        col=coul[ median.df$Group.1[!is.na(median.df$x)] ])
# plot(`x`~Group.1, data=med1.df, ylim=c(0,6), col = "gray10")
# points(`x`~Group.1, data=med2.df, ylim=c(0,6), col = "gray30", pch=2)
# points(`x`~Group.1, data=med3.df, ylim=c(0,6), col = "gray50", pch=15)
# 
# #legend("topright", legend=names(coul[ median.df$Group.1[!is.na(median.df$x)] ]),
# #       fill=coul[ median.df$Group.1[!is.na(median.df$x)] ], cex=0.5, bty="n", 
# #       x=0, y=max(df$`C||`[!is.na(df$grp)], na.rm=T), horiz=F, ncol=6)
# 
# dev.off()
# 
# #pdf(file=paste0(out.dir, "/", out.id, "_box.pdf"), height=10, width=15)
# 
# #dev.off()

# rm(list=ls()); gc()