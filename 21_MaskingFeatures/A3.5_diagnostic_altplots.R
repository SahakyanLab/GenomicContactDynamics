################################################################################
# Alternative diagnostic plots
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/21_MaskingFeatures")
diag.dir = paste0(wk.dir, "/out_diagnostic_plots")
out.dir = paste0(wk.dir, "/out_diagnostic_altplots") 
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
chr.id = "chrALL"
bin.len = 40000
kmer.len = 7
affix.binkmer = "_hg19_rm"
yval.mult = 1000 # y values in violin plots
unmskdThresh = 0.5
Cps = 1:21
lty.val = c(UMlen="solid", KCNT="dashed")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
library(RColorBrewer)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(gcb, "_", chr.id, "_", "kmerlen", kmer.len, "_binlen", bin.len, 
                   affix.binkmer, "_maskedfrlessthan", unmskdThresh)
plot.suff <- paste0("onlyBins>=", unmskdThresh, "unmaskedseq_binswithmissingseqNotIncluded_forRegionwiseplotUniqueBinsNotTaken",
                    "_pointAtMedian_trim=TsoTailUpToRangeOfData_YvalsDividedBy", yval.mult)

# Get data
mx <- sapply(X=chr.v, simplify=F, FUN=function(chr){
  load(paste0(diag.dir, "/", chr, "_", out.name, "_binval_vsCp.RData"))
  return(bins.mx)
})
df <- as.data.frame(do.call("rbind", mx))
rm(mx)

# Plot

df$Cp <- factor(as.character(df$Cp), levels=as.character(Cps))
df <- reshape2::melt(df)

coul <- colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(length(Cps))
coul <- setNames(object=c("black", coul), nm=c("0", as.character(Cps)))
df$line.col.grp <- as.character(df$Cp)
df$line.col.grp[df$variable == "UMlen"] <- "0"
df$line.col.grp <- factor(df$line.col.grp, levels=names(coul))

p <- ggplot(data=df, aes(x=value)) +
  geom_density(aes(y= ..scaled.., col=line.col.grp, lty=variable), size=1.5) +
  scale_x_continuous(labels=function(x) x / yval.mult) + 
  scale_y_continuous(breaks=c(0,1)) + 
  scale_color_manual(values=coul[levels(df$line.col.grp)]) +
  scale_linetype_manual(values=lty.val[levels(df$variable)]) + 
  labs(title=paste0(out.name, "\n", plot.suff)) + 
  guides(col=guide_legend(nrow=2)) + 
  bgr2 +
  theme(legend.text=element_text(size=5), legend.title=element_text(size=5), 
        legend.key.size=unit(0.1, "cm"), legend.position="bottom", 
        axis.text.x=element_text(size=7, angle=90), plot.title=element_text(size=2),
        strip.text=element_blank(),
        panel.background=element_rect(colour="gray22", size=1, fill="gray70")) + 
  facet_wrap(.~Cp, ncol=6) 

ggsave(filename=paste0(out.dir, "/", out.name, "_KCNT_numUMChar_vsCp_altplot.pdf"), plot=p,
       width=10, height=10, units="in")


# rm(list=ls()); gc()