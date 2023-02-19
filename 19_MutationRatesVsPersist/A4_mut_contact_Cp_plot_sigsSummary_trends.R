################################################################################
# Barplot of trends of Cp means relative to mean of previous Cp, make plot
# combining sig and calc (final) and plot for each calc (for checking variations
# across calc). I think combining calc is okay because we expect trend to be 
# similar between calcs.
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
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/19_MutationRatesVsPersist")
src.dir = paste0(wk.dir, "/z_ignore_git/out_mut_contact_Cp_plot")
out.dir = paste0(wk.dir, "/z_ignore_git/out_mut_contact_Cp_plot_sigsSummary_trends")

mutsig.file = paste0(data.dir, "/signal_mutSig/out_samplesForSignature/donorlist_signatureExposureRAW.csv")
### OTHER SETTINGS #############################################################
Cp.v = 1:21

mut.data.id = "ijfnxmean_donor_centric_PCAWG_Hg19"
sig.filter.id = "sigEperclimits_nosampfilter_ijmut"

mut.calcs = c("numWTSEQ", "Tmut", "Tmutnorm", "Nmsite", "Nmsitenorm", "TmutDIVNmsite")
mut.types = "All" #c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mut.types = gsub(">", "To", mut.types, fixed=T)
mut.locs = "exon" #c("exon", "intron", "intergenic", "intron_intergenic", 
                  #  "exon_intron_intergenic_intergenic.excl")

# Get list of signatures
sig.df = read.csv(file=mutsig.file)[,-(1:7)]
mut.sigs = "RefSig.1" #c("RefSig.MMR1_RefSig.MMR2", colnames(sig.df))
rm(sig.df)

fdr.cutoff = 0.05
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))

library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
library(cowplot)
library(ggpubr)
library(ggsci)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
combi.df <- expand.grid(type=mut.types, loc=mut.locs, calc=mut.calcs, sig=mut.sigs, 
                        stringsAsFactors=F)
combi.len <- length(combi.df$type)

# Data for plotting

df <- foreach(c.ind=1:combi.len, .combine="rbind"
                   
  ) %do% {
  
  type <- combi.df$type[[c.ind]]
  loc <- combi.df$loc[[c.ind]]
  calc <- combi.df$calc[[c.ind]]
  sig <- combi.df$sig[[c.ind]]
  
  src.nme <- paste0(calc, "_", mut.data.id, "_", type, "_", sig, "_", loc, "_", sig.filter.id)
  
  # Load pairwisedifftest.RData
  load(paste0(src.dir, "/chrALL_", src.nme, "_pairwisedifftest.RData"))
  
  # Summary of trends looking at if value at a Cp greater/less/ns that value of Cp - 1
  
  trend.mx <- matrix(data=NA, nrow=length(Cp.v) + 1, ncol=4,
                     dimnames=list(c("0", Cp.v), c("Cpafter", "BH", "mean.val", "mean.CpMINUSPrevCp")))
  
  # BH adjusted p-values
  BH.df <- reshape2::melt(t(TEST$pmw$p.value)) # transpose so Var1 is Cpbefore, Var2 is Cpafter
  BH.df <- BH.df[!is.na(BH.df$value),]
  colnames(BH.df) <- c("Cpbefore", "Cpafter", "BH")
  # Get pairs of consecutive Cps
  BH.df <- BH.df[ (as.numeric(BH.df$Cpafter) - as.numeric(BH.df$Cpbefore)) == 1, ]
  
  if( !identical(as.character(BH.df$Cpafter), as.character(Cp.v)) ){
    rm(trend.mx)
    stop(paste0(src.nme, ": Error in getting BH adjusted p-valus data."))
  } 
  
  # Add other data to trend.mx
  
  Cpafter <- c("0", Cp.v)
  trend.mx[,"Cpafter"] <- as.numeric(Cpafter)
  trend.mx[-1,"BH"] <- BH.df$BH
  trend.mx[,"mean.val"] <- do.call("rbind", TEST$meanmed)[Cpafter,"MEAN"]
  trend.mx[-1,"mean.CpMINUSPrevCp"] <- unname(diff(trend.mx[,"mean.val"]))
  
  trend.df <- cbind.data.frame(calc=calc, sig=sig, type=type, loc=loc, trend.mx)
  
  return(trend.df)
  
}

# Alter df for plotting
 
df <- df[df$Cpafter != 0,]
df$Cpafter <- factor(as.character(df$Cpafter), levels=as.character(Cp.v))

# Add other pertinent columns

# gr (Cp mean greater than mean of previous Cp), le (less), eq (equal)
df$trend.CpMINUSPrevCp <- df$mean.CpMINUSPrevCp
df$trend.CpMINUSPrevCp[ df$mean.CpMINUSPrevCp > 0 ] <- "gr"
df$trend.CpMINUSPrevCp[ df$mean.CpMINUSPrevCp < 0 ] <- "le"
df$trend.CpMINUSPrevCp[ df$mean.CpMINUSPrevCp == 0 ] <- "eq"

#  
df$significance <- df$BH
df$significance[ df$BH < fdr.cutoff & !is.na(df$BH) ] <- "s"
df$significance[ df$BH >= fdr.cutoff & !is.na(df$BH) ] <- "ns"
df$trend.group <- paste0(df$trend.CpMINUSPrevCp, ".", df$significance)
df$trend.group <- factor(as.character(df$trend.group),  
                         levels=c("gr.s", "gr.ns", "eq.s", "eq.ns", "le.s", "le.ns"))
cols.trend.group <- setNames(object=c("#E64B35FF", "#F39B7FFF", "#7E6148FF", 
                                      "#B09C85FF", "#3C5488FF", "#4DBBD5FF"),
                             nm=levels(df$trend.group))

# Variables needed for plots

out.id.general <- paste0("chrALL_", mut.data.id, "_", sig.filter.id)
out.id.fin <- paste0(out.id.general, "_barplot")

pd <- position_dodge(0.1)
types <- unique(df$type)
type.len <- length(types)
loc.len <- length(unique(df$loc))
plot.title <- paste0(out.id.general, "_onlyShownAreTrend.groupPresent")

# Plot, type vs. loc panel, combine sig, calc

P.LST <- sapply(1:type.len, simplify=F, FUN=function(t.ind){
  
  typ <- types[[t.ind]]
  df.typ <- df[df$type == typ,]
  
  trend.group.present <- levels(df$trend.group)[levels(df$trend.group) %in% df$trend.group]
  
  p <- ggplot(df.typ, aes(Cpafter, fill=trend.group)) +
    geom_bar(aes(fill=trend.group)) +
    scale_fill_manual(limits=trend.group.present,
                      values=cols.trend.group[trend.group.present]) +
    labs(title=plot.title) + 
    bgr1 +
    facet_grid(.~loc)
  return(p)
  
})

# Plot, type vs. loc panel, combine sig, calc

P.LST <- sapply(1:type.len, simplify=F, FUN=function(t.ind){
  
  typ <- types[[t.ind]]
  df.typ <- df[df$type == typ,]
 
  trend.group.present <- levels(df$trend.group)[levels(df$trend.group) %in% df$trend.group]
 
  p <- ggplot(df.typ, aes(Cpafter, fill=trend.group)) +
    geom_bar(aes(fill=trend.group)) +
    scale_fill_manual(limits=trend.group.present,
                      values=cols.trend.group[trend.group.present]) +
    labs(title=plot.title) + 
    bgr1 +
    facet_grid(.~loc)
  return(p)
  
})

# Version with all details
p.arr <- ggarrange(plotlist=P.LST, nrow=type.len, common.legend=T, legend="top")
ggexport(p.arr, width=loc.len * 5, height=type.len * 5, 
         filename=paste0(out.dir, "/calcCombined_", out.id.fin, "_withLabels.pdf"))

# Final minimal version
P.LST <- lapply(P.LST, FUN=function(p){
  
  p <- p + 
    labs(title=NULL) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
          axis.title.y=element_blank(), legend.position="none",
          strip.text=element_blank()) 
  return(p)
  
})

p.arr <- plot_grid(plotlist=P.LST, nrow=type.len, align="hv", byrow=T)
save_plot(filename=paste0(out.dir, "/calcCombined_", out.id.fin, "_noLabels.pdf"), plot=p.arr,
          base_height=type.len * 5, base_width=loc.len * 5)

# Plot, type vs. loc panel per calc, combine sig

for( calc in unique(df$calc) ){
  
  is.calc <- df$calc == calc
  
  P.LST <- sapply(1:type.len, simplify=F, FUN=function(t.ind){
    
    typ <- types[[t.ind]]
    df.typ <- df[df$type == typ & is.calc,]
    
    trend.group.present <- levels(df$trend.group)[levels(df$trend.group) %in% df$trend.group]
    
    p <- ggplot(df.typ, aes(Cpafter, fill=trend.group)) +
      geom_bar(aes(fill=trend.group)) +
      scale_fill_manual(limits=trend.group.present,
                        values=cols.trend.group[trend.group.present]) +
      labs(title=plot.title) + 
      bgr1 +
      facet_grid(.~loc)
    return(p)
    
  })
  
  # Version with all details
  p.arr <- ggarrange(plotlist=P.LST, nrow=type.len, common.legend=T, legend="top")
  ggexport(p.arr, width=loc.len * 5, height=type.len * 5, 
           filename=paste0(out.dir, "/", calc, "_", out.id.fin, "_withLabels.pdf"))
  
}

# rm(list=ls()); gc()



