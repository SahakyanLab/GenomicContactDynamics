################################################################################
# Alternative plot for complementarity values
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

options(warnPartialMatchDollar = T)
options(warn = 1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
lib = paste0(home.dir, "/DPhil/lib")
constraints.id = "GfreeSingleNorm" #"hg19_rm_GfreeSingleNorm"
compl.dir = paste0(wk.dir, "/out_constraints_", constraints.id, "/merged_final")
out.dir = paste0(wk.dir, "/out_constraints_alternative_plot")
### OTHER SETTINGS #############################################################
chr = "chrALL" 
gcb = "min2Mb"
type = "kmer"
affix = ""
Cps = 1:21
ylim.val = list(CIIkmer=c(-2.5,0), CIIalign=NULL, # CIIkmer=c(-2.5,0)
                Gfreekmer=c(-1.2,-0.4), sdDifferencekmer=c(0,0.00030))

out.id = constraints.id # gap_effect

DATA.PATH = list(paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"),
                 paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"),
                 paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"),
                 paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
GAP.RNG = list(NULL, c(50,50), c(50,125), c(125,Inf)) # j - i - 1, closed range, set to NULL if not filtering
lty.val = rep("dashed", 4) #c("solid", "dashed", "dotted", "dotdash")
shp.val = 21:24 #c(15,1,2,4)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(Rmisc)
library(data.table)
library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
LENS <- c(length(DATA.PATH), length(GAP.RNG), length(lty.val), length(shp.val))
if( any(LENS != length(DATA.PATH)) ){
  rm(GAP.RNG, DATA.PATH)
  stop("Check lengths of objects, whose lengths should be equal.")
}

if( !is.null(GAP.RNG) ){
  GAP.RNG.ID <- lapply(X=GAP.RNG, FUN=function(gap.rng){
    if(is.null(gap.rng)){
      return("NULLgap")
    } else {
      return( paste0("gap", paste(gap.rng, collapse="To"), "bins") )
    }
  })
}

#

df <- sapply(X=1:length(DATA.PATH), simplify=F, FUN=function(ind){
  
  load(DATA.PATH[[ind]])
  gap.rng <- GAP.RNG[[ind]]
  if( !is.null(gap.rng) ){
    print(paste0(GAP.RNG.ID[[ind]], ": Filtering contacts based on this closed gap range."))
    gaps <- CII.MX[,"j"] - CII.MX[, "i"] - 1
    CII.MX <- CII.MX[ gaps >= gap.rng[[1]] & gaps <= gap.rng[[2]], ]
    rm(gaps)
  } else {
    print(paste0(GAP.RNG.ID[[ind]], ": No filtering of contacts based on gap."))
  }
  CII.MX <- data.frame(CII.MX, dta.id=GAP.RNG.ID[[ind]], check.names=F)
  return(CII.MX)
  
})

# Plot

df <- do.call("rbind", df)
setnames(x=df, old="C||", new="CII")

# Remove contacts with no Cp
df <- df[!is.na(df$Cp),]

df$Cp <- factor(as.character(df$Cp), levels=as.character(Cps))
df$dta.id <- factor(as.character(df$dta.id), levels=unlist(GAP.RNG.ID))

compl.types <- setdiff(colnames(df), c("i", "j", "Cp", "dta.id"))
out.name <- paste0(out.id, "_", chr, "_", type, "_", gcb, affix)
dta.ids <- levels(df$dta.id)
gaps.noNACp <- df$j - df$i - 1 

coul <- colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(length(Cps))
pd <- position_dodge(0.7)

for(c.type in compl.types){
  
  # Calculate actual range of gap values and add to title
  
  is.nonNA.c.type <- !is.na(df[[c.type]])
  df.tmp <- df[is.nonNA.c.type,]
  
  ACT.GAP.RNG.ID <- sapply(X=dta.ids, simplify=T, FUN=function(dta.id){
    is.dta <- df$dta.id == dta.id
    gap.rng.noNA <- range( gaps.noNACp[is.nonNA.c.type & is.dta] ) 
    return( paste0("gap", paste(gap.rng.noNA, collapse="To"), "bins") )
  })
 
  plot.title <- paste0(out.name, "\nactualgapAfterNAxyvaldropped_", paste(ACT.GAP.RNG.ID, collapse="_"))
  
  ## PLOTS
  
  # Plot 1 - Median + Boxplot whiskers
  
  p.bp <- ggplot(data=df.tmp, aes_string(x="Cp", y=c.type)) +
    scale_y_continuous(limits=ylim.val[[paste0(c.type, type)]]) +
    scale_fill_manual(values=coul) +
    scale_linetype_manual(values=lty.val) +
    scale_shape_manual(values=shp.val) +
    guides(fill="none") +
    bgr2 +
    theme(aspect.ratio=0.8, plot.title=element_text(size=5),
          legend.position="bottom", legend.text=element_text(size=5),
          legend.title=element_text(size=5))
    
  p <- p.bp + 
    stat_boxplot(geom="errorbar", aes(lty=dta.id), col="gray20", width=0, position=pd) + 
    stat_summary(aes(shape=dta.id, fill=Cp), fun="median", position=pd) +
    geom_vline(xintercept=seq(1.5, length(levels(df$Cp))-0.5, 1), col="gray90", lty="solid", lwd=0.5) + 
    labs(title=paste0(plot.title, "\ngapequalsjMINUSiMINUS1_pointsAtMedian_errorbarSameAsDefaultBoxplotWhiskers"))
    
  ggsave(filename=paste0(out.dir, "/", out.name, "_", c.type, "_medianPlusBPwhiskers.pdf"), 
                         width=10, height=8, plot=p)
  
  # Plot 2 - Median + Boxplot hinges
  
  p <- p.bp + 
    stat_summary(aes(shape=dta.id, fill=Cp, lty=dta.id), col="gray20", fun.y="median", position=pd,
                 fun.min=function(a){quantile(a, 0.25)}, fun.max=function(a){quantile(a, 0.75)}) +
    geom_vline(xintercept=seq(1.5, length(levels(df$Cp))-0.5, 1), colour="gray90", lty="solid", lwd=0.5) + 
    labs(title=paste0(plot.title, "\ngapequalsjMINUSiMINUS1_pointsAtMedian_errorbarSameAsDefaultBoxplotHinges")) 
  
  if(c.type == "CII"){ # Modify
    p <- p + 
    scale_y_continuous(limits=c(-1.9, -0.5), breaks=seq(-1.8, -0.6, by=0.2))
  } 
  
    ggsave(filename=paste0(out.dir, "/", out.name, "_", c.type, "_medianPlusBPhinges.pdf"), 
           width=10, height=8, plot=p)
  
  # Plot 3 - Mean + 95% CI
  
  df.summ.SE <- summarySE(df.tmp, measurevar=c.type, groupvars=c("Cp","dta.id"), na.rm=F, .drop=F)
  save(df.summ.SE, file=paste0(out.dir, "/", out.name, "_", c.type, "_summary_statistics.RData"))
  
  setnames(x=df.summ.SE, old=c.type, new="value")
  
  p <- ggplot(data=df.summ.SE, aes(x=Cp, y=value)) + 
    geom_errorbar(aes(ymin=value - ci, ymax=value + ci, lty=dta.id), 
                  width=0, linewidth=0.6, position=pd) + 
    stat_summary(aes(shape=dta.id, fill=Cp), fun="mean", position=pd) +
    geom_vline(xintercept=seq(1.5, length(levels(df$Cp))-0.5, 1), colour="gray90", lty="solid", lwd=0.5) + 
    labs(title=paste0(plot.title, "\ngapequalsjMINUSiMINUS1_pointsAtmean_errorbarAt95PercCI"),
         y=c.type) +
    scale_y_continuous(limits=ylim.val[[paste0(c.type, type)]]) +
    scale_fill_manual(values=coul) + 
    scale_linetype_manual(values=lty.val) + 
    scale_shape_manual(values=shp.val) + 
    guides(fill="none") + 
    bgr2 +
    theme(aspect.ratio=0.8, plot.title=element_text(size=5), 
          legend.position="bottom", legend.text=element_text(size=5),
          legend.title=element_text(size=5))
  
  ggsave(filename=paste0(out.dir, "/", out.name, "_", c.type, "_meanPlus95PercCI.pdf"), 
         width=10, height=8, plot=p)
  
  #
  
  rm(is.nonNA.c.type, ACT.GAP.RNG.ID, plot.title, p, df.tmp, df.summ.SE)
  gc()
  
  print(paste0(c.type, " done!"), quote=F)
  
}

# rm(list=ls()); gc()

