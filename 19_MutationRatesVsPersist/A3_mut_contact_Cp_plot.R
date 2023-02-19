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
src.dir = paste0(wk.dir, "/z_ignore_git/out_mut_contact_Cp_plotdata")
out.dir = paste0(wk.dir, "/z_ignore_git/out_mut_contact_Cp_plot")

mutsig.file = paste0(data.dir, "/signal_mutSig/out_samplesForSignature/donorlist_signatureExposureRAW.csv")
### OTHER SETTINGS #############################################################
Cp.v = 1:21

mut.data.id = "donor_centric_PCAWG_Hg19"
sig.filter.id = "sigEperclimits_nosampfilter_ijmut"

chrs = paste0("chr", c(21:22))
nCPU = 1
mut.calcs = c("numWTSEQ", "Tmut", "Tmutnorm", "Nmsite", "Nmsitenorm", "TmutDIVNmsite")
mut.types = "All" #c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mut.types = gsub(">", "To", mut.types, fixed=T)
mut.locs = "exon" #c("exon", "intron", "intergenic", "intron_intergenic", 
                  #  "exon_intron_intergenic_intergenic.excl")

# Get list of signatures
sig.df = read.csv(file=mutsig.file)[,-(1:7)]
mut.sigs = "RefSig.1" #c("RefSig.MMR1_RefSig.MMR2", colnames(sig.df))
rm(sig.df)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))

library(Rmisc)

library(car) # ANOVA for unbalanced dataset
source(paste0(lib, "/doVarTest.R")) # Update deva copy
source(paste0(lib, "/doCorTest.R")) # Update deva copy
source(paste0(lib, "/compareManyDist.R"))  # Update deva copy

library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
library(cowplot)
library(ggpubr)
library(ggsci)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrs.len <- length(chrs)

combi.df <- expand.grid(type=mut.types, loc=mut.locs, calc=mut.calcs, sig=mut.sigs, 
                        stringsAsFactors=F)
combi.len <- length(combi.df$type)

DF.STAT <- foreach(c.ind=1:combi.len, .combine="rbind"
                   
) %do% {
  
  type <- combi.df$type[[c.ind]]
  loc <- combi.df$loc[[c.ind]]
  calc <- combi.df$calc[[c.ind]]
  sig <- combi.df$sig[[c.ind]]
  
  # Per calc, sig, type, loc, combine data from all chr -> df
  
  src.nme <- paste0(calc, "_", mut.data.id, "_", type, "_", sig, "_", loc, "_", sig.filter.id)
  
  toExport <- c("src.dir", "chrs", "src.nme", "calc", "sig", "type", "loc")
  #### PARALLEL EXECUTION #########
  df <- foreach(itr=isplitVector(1:chrs.len, chunks=nCPU), 
                .inorder=F, .combine="rbind",
                .export=toExport, .noexport=ls()[!ls()%in%toExport]
                
  ) %op% {
    
    chunk <- sapply(itr, simplify=F, FUN=function(i){
      
      # Load IJ.MUT
      load(paste0(src.dir, "/", chrs[[i]], "_", src.nme, ".RData"))
      return(IJ.MUT)
      
    })
    return( do.call("rbind", chunk))
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  df <- as.data.frame(df)
  df$Cp <- factor(as.character(df$Cp), levels=c(as.character(Cp.v)))
  
  # df summary statistics
  df.stat.Cp1To21 <- summarySE(data=df, measurevar="value", groupvars="Cp",
                               na.rm=T, .drop=F)

  # Add summary statistics for Cp=0 (all long-range contacts)
  df.stat.Cp0 <- summarySE(data=df, measurevar="value", groupvars=NULL, na.rm=T, .drop=F)[,-1]
  df.stat.Cp0 <- cbind.data.frame(Cp=factor("0", level="0"), df.stat.Cp0)

  C0To21.meds <- c(median(df$value, na.rm=T), 
                   aggregate(df$value, by=list(df$Cp), median, na.rm=T)$x)
  
  # Final df.stat
  df.stat <- cbind(calc=calc, sig=sig, type=type, loc=loc, 
                   rbind(df.stat.Cp0, df.stat.Cp1To21),
                   med=C0To21.meds)
  
  save(df, file=paste0(out.dir, "/chrALL_", src.nme, ".RData"))
  save(df.stat, file=paste0(out.dir, "/chrALL_", src.nme, "_summ_stat.RData"))
  
  ## P-values
  
  # P-values, ANOVA/KW and correlation tests
  
  df <- na.omit(df)
  
  try(doVarTest( xval=df$value, grp=df$Cp, out.dir=out.dir, out.name=paste0("chrALL_", src.nme) ))
  
  try(doCorTest( xval=as.numeric(as.character(df$Cp)), yval=df$value, alt="two.sided",
                 exactpval=F, out.dir=out.dir, out.name=paste0("chrALL_", src.nme) ))
  
  # Add all values as Cp=0
  
  df.Cp0 <- df
  df.Cp0$Cp <- factor("0", levels="0")
  df <- rbind(df.Cp0, df)
  
  try(compareManyDist( xval=df$value, grp=df$Cp, alt="two.sided", out.dir=out.dir, 
                       out.name=paste0("chrALL_", src.nme) ))
  rm(df)
  
  return(df.stat)
  
} # combi.len foreach loop end

## PLOTS

backup <- DF.STAT
DF.STAT.dummy <- DF.STAT                 # REMOVE
DF.STAT.dummy$loc <- "intron"
DF.STAT.dummy$type <- "All" # REMOVE
DF.STAT <- rbind(DF.STAT, DF.STAT.dummy) # REMOVE

# Variables needed for plots

out.id.general <- paste0(mut.data.id, "_", sig.filter.id)

calc.type.combi.df <- expand.grid(calc=unique(DF.STAT$calc), type=unique(DF.STAT$type), stringsAsFactors=F)
ct.len <- length(calc.type.combi.df[,1])
calc.len <- length(unique(DF.STAT$calc))
type.len <- length(unique(DF.STAT$type))

type.shapes <- setNames(object=c(4,15,2,18,0,16,17),
                        nm=c("All", "CToA", "CToG", "CToT", "TToA", "TToC", "TToG"))
pd <- position_dodge(0.1)

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
      geom_point(aes(col=loc), size=2, position=pd, shape=type.shapes[[type]]) + 
      scale_color_npg() + 
      labs(title=paste0(sig, "_", calc, "_", type), col=paste0(out.id.sig, "_loc")) + 
      bgr1 
    return(p)
    
  })
  
  out.id.fin <- paste0(out.id.sig, "_meanPlus95PercCI")
  
  # Version with all details

  p.arr <- ggarrange(plotlist=P.LST, ncol=calc.len, common.legend=T, legend="top")
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