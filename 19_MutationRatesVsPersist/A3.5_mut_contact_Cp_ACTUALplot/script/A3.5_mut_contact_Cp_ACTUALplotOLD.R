################################################################################
# Per sig, calc, type, loc, generate df (with data from all chr), df.stat (with
# summary statistics of df including Cp=0 or all long-range contacts), calculate
# P-values comparing Cp distributions and calc vs. Cp, generate average calc
# values vs. Cp per signature (but mainly for nosampfilter). 
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
src.dir = paste0(wk.dir, "/out_mut_contact_Cp_plot_fast")
out.dir = paste0(wk.dir, "/out_mut_contact_Cp_ACTUALplot_fast")

mutsig.file = paste0(data.dir, "/signal_mutSig/out_samplesForSignature/donorlist_signatureExposureRAW.csv")
### OTHER SETTINGS #############################################################
mut.data.id = "ijfnxmean_donor_centric_PCAWG_Hg19"
sig.filter.id = "sigEperclimits_nosampfilter_ijmut" #"sigEperclimits_1000rawInf_ijmut"

nCPU = 3 # Number of combinations
#mut.calcs = c("numWTSEQ", "Tmut", "Nmsite", "Tmutnorm", "Nmsitenorm", "TmutDIVNmsite")
mut.calcs = c("Tmutnorm", "Nmsitenorm", "TmutDIVNmsite", "Tmut", "Nmsite")
mut.types = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mut.types = gsub(">", "To", mut.types, fixed=T)
mut.locs = c("exon", "intron", "intergenic", "intron_intergenic", "exon_intron_intergenic_intergenic.excl")

# Get list of signatures
sig.df = read.csv(file=mutsig.file)
mut.sigs = colnames(sig.df)[grepl("RefSig.", colnames(sig.df))]
mut.sigs = "RefSig.1" #c(mut.sigs, "RefSig.MMR1_RefSig.MMR2")
rm(sig.df)

# Cp MEAN or MED (median)?
average.fnx = "MED"
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

#### PARALLEL EXECUTION #########
DF.STAT <- foreach(itr=isplitVector(1:combi.len, chunks=nCPU), .inorder=F, .combine="rbind"
#DF.STAT <- foreach(c.ind=1:combi.len, .combine="rbind"
                   
) %op% {
  
  chunk <- sapply(X=itr, simplify=F, FUN=function(c.ind){
    
    type <- combi.df$type[[c.ind]]
    loc <- combi.df$loc[[c.ind]]
    calc <- combi.df$calc[[c.ind]]
    sig <- combi.df$sig[[c.ind]]
    
    # Per calc, sig, type, loc, combine data from all chr -> df
    src.nme <- paste0(calc, "_", mut.data.id, "_", type, "_", sig, "_", loc, "_", sig.filter.id)
    load(file=paste0(src.dir, "/chrALL_", src.nme, "_summ_stat.RData"))
    return(df.stat)
    
  })
  
  return( do.call("rbind", chunk) )
  
} # combi.len foreach loop end
### END OF PARALLEL EXECUTION ###

if(average.fnx == "MED"){
  DF.STAT$value <- DF.STAT$med
  DF.STAT$ci <- 0
}

## PLOTS

# Variables needed for plots

out.id.general <- paste0(mut.data.id, "_", sig.filter.id)

calc.type.combi.df <- expand.grid(calc=unique(DF.STAT$calc), type=unique(DF.STAT$type), stringsAsFactors=F)
ct.len <- length(calc.type.combi.df[,1])
calc.len <- length(unique(DF.STAT$calc))
type.len <- length(unique(DF.STAT$type))

#type.shapes <- setNames(object=c(4,15,2,18,0,16,17),
#                        nm=c("All", "CToA", "CToG", "CToT", "TToA", "TToC", "TToG"))
pd <- position_dodge(0.3)

# (1) Per sig panel (mut.types vs. calc) -> mainly for nosampfilter

for(sig in mut.sigs){
  
  is.sig <- DF.STAT$sig == sig
  out.id.sig <- paste0(sig, "_", out.id.general)
  
  P.LST <- sapply(1:ct.len, simplify=F, FUN=function(ct){
    
    calc <- calc.type.combi.df$calc[[ct]]
    type <- calc.type.combi.df$type[[ct]]
    is.ct <- DF.STAT$calc == calc & DF.STAT$type == type
    
    p <- ggplot(DF.STAT[is.sig & is.ct,], aes(x=Cp, y=value)) +
      geom_errorbar(aes(col=loc, ymin=value - ci, ymax=value + ci), width=0.4, linewidth=0.6, 
                    position=pd) +
      geom_point(aes(col=loc), size=3, position=pd, shape=1) + 
      scale_color_npg() + 
      labs(title=paste0(sig, "_", calc, "_", type), col=paste0(out.id.sig, "_loc")) + 
      bgr1 + 
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7))
    return(p)
    
  })
  
  out.id.fin <- paste0(out.id.sig, "_meanPlus95PercCI")
  if(average.fnx == "MED"){
    out.id.fin <- paste0(out.id.sig, "_medianNoCI")
  }
  
  # Version with all details

  p.arr <- ggarrange(plotlist=P.LST, ncol=calc.len, nrow=type.len, common.legend=T, legend="top")
  ggexport(p.arr, width=calc.len * 5, height=type.len * 5, 
           filename=paste0(out.dir, "/", out.id.fin, "_withLabels.pdf"))

  # Final minimal version
  P.LST <- lapply(P.LST, FUN=function(p){
    
    p <- p + 
      labs(title=NULL) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.title.y=element_blank(), legend.position="none")
    return(p)
    
  })
  
  p.arr <- plot_grid(plotlist=P.LST, ncol=calc.len, align="hv", byrow=T)
  save_plot(filename=paste0(out.dir, "/", out.id.fin, "_noLabels.pdf"), plot=p.arr,
            base_height=type.len * 5, base_width=calc.len * 5)
  
} # mut.sigs for loop end

# rm(list=ls()); gc()
