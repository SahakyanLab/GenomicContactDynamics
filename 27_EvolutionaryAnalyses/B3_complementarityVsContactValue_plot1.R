################################################################################
# Analysis plots comparing complementarity with various values/features
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses")
valmx.dir = paste0(wk.dir, "/out_complementarityVsContactValue") 
out.dir = paste0(wk.dir, "/out_complementarityVsContactValue_plot") 
state.group.file = paste0(wk.dir, "/Phylo-HMRF_state_group.csv")
group.col.file = paste0(wk.dir, "/Phylo-HMRF_state_group_col.csv")
group.display.file = paste0(wk.dir, "/Phylo-HMRF_state_group_display.csv")
### OTHER SETTINGS #############################################################
value.id = "genome_state_Phylo-HMRF_contact50K_norm_KR"
state.group = "Grp0"
display.group.id = "Grp0_disp2"
gcb = "min2Mb"
group.order.basis = "kmer.CII" # Use "II" for "||"
chr.v = paste0("chr", c(1,4:22))
bin.len = 50000
gap.range.bins.closed = c(0,39) # j - i - 1, No gap filtering if NULL
out.id = paste0("hg38.108_bin", bin.len, "bp_", value.id, "_gaprangebins_", 
                gap.range.bins.closed[1], "_", gap.range.bins.closed[2],
                "_state", state.group, "_", display.group.id)  
coul.spec <- c(`GSM3685694_panTro5_mapping`="cyan2", 
               `GSM3685695_panPan2_mapping`="blue", 
               `GSM3685696_gorGor4_mapping`="darkgreen", 
               `hg38.downsample`="hotpink") # hotpink for human
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(vioplot)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(ggsci)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareManyDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
df <- sapply(chr.v, simplify=F, FUN=function(chr){
  
  load(paste0(valmx.dir, "/", chr, "_", gcb, "_", value.id, ".RData"))
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

# Change "||" to "II" in column names for ggplot
colnames(df) <- gsub(pattern="||", replacement="II", x=colnames(df), fixed=T)

# Identify PHYLO-HMRF state groupings 
grp.df <- read.csv(state.group.file, header=T, stringsAsFactors=F)
is.nonNAstate <- !is.na(df$PHYLOHMRFstate)
df$grp[is.nonNAstate] <- grp.df[[state.group]][ as.numeric(as.character(df$PHYLOHMRFstate))[is.nonNAstate] ]

# Get only groups you want to be displayed
disp.df <- read.csv(group.display.file, header=T, stringsAsFactors=F, na.strings="")
group.display <- na.omit(disp.df[[display.group.id]])
df <- df[df$grp %in% group.display,] 

# Reorder states - increasing median complementarity measure and order as levels
median.df <- aggregate(x=df[[group.order.basis]], by=list(as.character(df$grp)), FUN=median, na.rm=T)
median.df <- median.df[order(median.df$x, na.last=F), ]
group.present <- as.character(median.df$Group.1)

df$grp <- factor(as.character(df$grp), 
                 levels=c(setdiff(group.display, group.present), group.present))
median.df$Group.1 <- factor(group.present, levels=levels(df$grp))
colnames(median.df) <- c("grp", "value")

# Plot

col.df <- read.csv(group.col.file, header=T, stringsAsFactors=F)
col.df$Col <- paste0("#", col.df$Col)
col.df <- col.df[col.df$Grp %in% levels(df$grp),]
coul.grp <- setNames(col.df$Col, nm=col.df$Grp)

#

compl.nmes <- colnames(df)
compl.nmes <- compl.nmes[ grepl(pattern="kmer.|align.", compl.nmes) ]
P.LST <- list()
for(nme in compl.nmes){
  
  df.tmp <- na.omit(df[,c("grp", nme)])
  compareManyDist(xval=df.tmp[[nme]], grp=df.tmp$grp, 
                  alt="two.sided", out.dir=out.dir, out.name=paste0(out.id, "_", nme))
  
  P.LST[[nme]] <- ggplot(data=df.tmp, aes_string(x="grp", y=nme)) +
    geom_violin(scale="width", aes(fill=grp), col="white", alpha=0.5, trim=T, lwd=2) +
    #geom_point(data=median.df, aes(x=grp, y=value), alpha=1, size=4, pch=3, stroke=1.1) +
    stat_boxplot(geom="errorbar", width=0, lty="dashed") + 
    stat_summary(fun="median", pch=1, position=position_dodge(0.5)) +
    scale_x_discrete(limits=levels(df.tmp$grp)) +
    scale_fill_manual(values=unname(coul.grp[group.present])) + 
    scale_colour_manual(values=unname(coul.grp[group.present])) + 
    labs(x=NULL, title=out.id) + 
    bgr2 + 
    theme(legend.position="none", 
          axis.text.x=element_text(size=7, angle=90),
          aspect.ratio=0.25) 
  
  rm(df.tmp)
  
}

#

KR.nmes <- colnames(df)[grepl(colnames(df), pattern="_KR", fixed=T)] 

#KR.df <- df[,c(KR.nmes, "grp")]
#KR.df <- reshape2::melt(KR.df, id="grp")
#colnames(KR.df) <- c("grp", "species", "value")
#KR.df$species <- factor(KR.df$species, levels=unique(KR.df$species))
#coul.nme <- gsub(pattern="_KR", replacement="", x=levels(KR.df$species))
# P.LST[["KR_cf_"]] <- ggplot(data=KR.df, aes(x=grp, y=log10(value))) +
#   geom_violin(scale="width", aes(col=grp), fill=adjustcolor("white", alpha.f=0.1), 
#               alpha=0.5, trim=T, lwd=1) +
#   stat_boxplot(geom="errorbar", width=0, lty="dashed") + 
#   stat_summary(fun="median", pch=1, position=position_dodge(0.5)) +
#   scale_x_discrete(limits=levels(KR.df$grp)) +
#   scale_fill_manual(values=unname(coul.grp[group.present])) + 
#   scale_colour_manual(values=unname(coul.grp[group.present])) + 
#   labs(x=NULL) + 
#   bgr2 + 
#   theme(legend.position="none", 
#         axis.text.x=element_text(size=7, angle=90),
#         aspect.ratio=0.25) +
#   facet_grid(species~.) 

ylim.val <- range(log10(df[,KR.nmes]), na.rm=T)
ylim.val <- c(floor(ylim.val[[1]]), ceiling(ylim.val[[2]]))
for(nme in KR.nmes){

  df.tmp <- na.omit(df[,c("grp", nme)])
  compareManyDist(xval=df.tmp[[nme]], grp=df.tmp$grp,
                  alt="two.sided", out.dir=out.dir, out.name=paste0(out.id, "_", nme))
  df.tmp[[nme]] <- log10(df.tmp[[nme]])

  P.LST[[nme]] <- ggplot(data=df.tmp, aes_string(x="grp", y=nme)) +
    geom_violin(scale="width", aes(col=grp), fill="white", alpha=0.1, trim=T, lwd=2) +
    stat_boxplot(geom="errorbar", width=0, lty="dashed") +
    stat_summary(fun="median", pch=1, position=position_dodge(0.5)) +
    scale_x_discrete(limits=levels(df.tmp$grp)) +
    scale_y_continuous(limits=ylim.val) + 
    scale_fill_manual(values=unname(coul.grp[group.present])) +
    scale_colour_manual(values=unname(coul.grp[group.present])) +
    labs(x=NULL, title=paste0(out.id, "y value log10 transformed")) +
    bgr2 +
    theme(legend.position="none",
          axis.text.x=element_text(size=7, angle=90),
          aspect.ratio=0.25)

  rm(df.tmp)

}

#

norm.df <- df[ !is.na(df$grp),
               colnames(df)[grepl(colnames(df), pattern="_norm|grp")] 
             ]
norm.df <- reshape2::melt(norm.df, id="grp")
colnames(norm.df) <- c("grp", "species", "value")
norm.df$species <- factor(norm.df$species, levels=unique(norm.df$species))

coul.nme <- gsub(pattern="_norm", replacement="", x=levels(norm.df$species))
P.LST[["norm_cf"]] <- ggplot(data=norm.df, aes(x=value)) +
  geom_density(aes(y= ..scaled.., col=species)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", -(x))) + 
  scale_y_continuous(breaks=c(0,1)) + 
  scale_color_manual(values=unname(coul.spec[coul.nme]), drop=F) + 
  labs(title=out.id) + 
  bgr2 +
  theme(legend.position="bottom", 
        legend.text=element_text(size=2), legend.title=element_text(size=2), 
        axis.text.x=element_text(size=7, angle=90),
        aspect.ratio=1.5) +
  facet_grid(. ~ grp) +
  coord_flip() 

#

cp.nmes <- colnames(df)[grepl(colnames(df), pattern=".consCp", fixed=T)] 
for(nme in cp.nmes){
  
  df.tmp <- na.omit(df[,c("grp", nme)])
  compareManyDist(xval=df.tmp[[nme]], grp=df.tmp$grp, 
                  alt="two.sided", out.dir=out.dir, out.name=paste0(out.id, "_", nme))
  rm(df.tmp)
  
}

cp.df <- df[ !is.na(df$grp), 
             colnames(df)[grepl(colnames(df), pattern=".consCp|grp")]
]
cp.df <- reshape2::melt(cp.df, id="grp")

variable.len <- length(levels(cp.df$variable))
P.LST[["consCp"]] <- ggplot(data=cp.df[!is.na(cp.df$grp),], aes(x=grp, y=value)) +
  stat_boxplot(geom="errorbar", aes(lty=variable), col="gray30", 
               position=position_dodge(0.5), width=0, lwd=0.5) + 
  stat_summary(fun="median", aes(pch=variable), position=position_dodge(0.5)) +
  scale_x_discrete(limits=levels(df$grp)) +
  scale_y_continuous(labels=function(y) sprintf("%.1f", -(y))) + 
  #scale_colour_manual( values=pal_npg("nrc")(variable.len) ) +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) + 
  scale_shape_manual(values=rep(1, times=variable.len)) + 
  labs(y="Cp", x=NULL, title=out.id) + 
  bgr2 + 
  theme(legend.position="top", 
        legend.text=element_text(size=2), legend.title=element_text(size=2), 
        axis.text.x=element_text(size=7, angle=90),
        aspect.ratio=0.25) 

# Save plots

for( p.nme in names(P.LST) ){
  
  print(paste0("Plotting ", p.nme, "..."), quote=F)
  ggsave(paste0(out.dir, "/", out.id, "_", p.nme, ".pdf"), plot=P.LST[[p.nme]], 
         width=15, height=5, units="in")

}

#pdf(paste0(out.dir, "/", out.id, "_", p.nme, ".pdf"), width=15, height=5)
#cowplot::plot_grid(plotlist=P.LST[[p.nme]], ncol=1, align='v')
#dev.off()

#p.arr <- ggarrange(plotlist=list(p.cp, p.cii, p.cf), nrow=3, ncol=1, legend=NULL)
#ggexport(p.arr, width=15, height=20, units="in",
#         filename=paste0(out.dir, "/", out.id, "_Cp_CII.pdf")) 

#ggsave(paste0(out.dir, "/", out.id, "_Cf.pdf"), plot=p.cf, width=15, height=10, units="in")

# rm(list=ls()); gc()
