################################################################################
# Title
# deva, R/3.5.0-newgcc, gcc/4.9.2
# deva, R/3.6.0-newgcc, gcc/4.9.2
# Mac, R/3.5.2, R/3.6.1
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
src.dir = paste0(wk.dir, "/out_mut_contact_Cp_plotdata")
out.dir = paste0(wk.dir, "/out_mut_contact_Cp_plot")

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

library(ggsci)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
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
  
  # TO DO - Save df all chr combined
  # TO DO - Save df.stat all chr combined
  # TO DO - P-values
  
  return(df.stat)
  
} # combi.len foreach loop end

# Plot

backup <- DF.STAT
DF.STAT.dummy <- DF.STAT                 # REMOVE
DF.STAT.dummy$loc <- "intron"
DF.STAT.dummy$type <- "CToA" # REMOVE
DF.STAT <- rbind(DF.STAT, DF.STAT.dummy) # REMOVE

out.id.ci <- paste0(mut.data.id, "_", sig.filter.id)

#

calc.sig.combi.df <- expand.grid(calc=mut.calcs, sig=mut.sigs, stringsAsFactors=F)
cs.len <- length(calc.sig.combi.df[,1])  

# Per sig panel (calc vs. mut.types) -> mainly for nosampfilter

shapes <- setNames(object=c(4,15,2,18,0,16,17),
                   nm=c("All", "CToA", "CToG", "CToT", "TToA", "TToC", "TToG"))

calc.type.combi.df <- expand.grid(calc=mut.calcs, type=mut.types, stringsAsFactors=F)
ct.len <- length(calc.sig.combi.df[,1])

calc.len <- length(mut.calcs)
type.len <- length(mut.types)

for(sig in mut.sigs){
  
  is.sig <- DF.STAT$sig == sig
  out.name.sig <- paste0(sig, "_", out.id.ci)
  
  P.LST <- sapply(1:ct.len, simplify=F, FUN=function(ct){
    
    calc <- calc.type.combi.df$calc[[ct]]
    type <- calc.type.combi.df$type[[ct]]
    is.ct <- DF.STAT$calc == calc & DF.STAT$type == type
    
    p <- ggplot(DF.STAT[is.sig & is.ct,], aes(x=Cp, y=value)) +
      geom_errorbar(aes(col=loc, ymin=value - ci, ymax=value + ci), width=0.4, linewidth=0.6, 
                    position=pd) +
      geom_point(aes(col=loc), size=2, position=pd, shape=shapes[[type]]) + 
      scale_color_npg() + 
      bgr1 
      
    return(p)
    
  })
  
  p.arr <- ggarrange(plotlist=P.LST, nrow=calc.len, ncol=type.len, common.legend=T)
  ggexport(p.arr, height=calc.len * 5, width=type.len * 5,
           filename=paste0(out.dir, "/", out.name.sig, ".pdf"))
  
  
} # mut.sigs for loop end


  for(calc in mut.calcs){
    
    # Plot mean value at Cp (plus ci) vs. Cp
    # Per calc, sig, mut types separated because of ci bars
    
    out.name <- paste0(calc, "_", sig, "_", out.id.ci)
    is.row.cs <- DF.STAT$calc == calc & DF.STAT$sig == sig
    pd <- position_dodge(0.1)
    
    p <- ggplot(DF.STAT[is.row,], aes(x=Cp, y=value)) +
      geom_errorbar(aes(col=loc, ymin=value - ci, ymax=value + ci), width=0.4, linewidth=0.6, 
                    position=pd) +
      geom_point(aes(col=loc, shape=type), size=2, position=pd) + 
      scale_color_npg() + 
      bgr1 +
      labs(title=out.name) + 
      facet_wrap(.~type, nrow=2) +
      theme(strip.text=element_blank(), legend.position="bottom")
    
    ggsave(filename=paste0(out.dir, "/", out.name, "_meanPlus95PercCI.pdf"), 
           width=40, height=20, plot=p)
    
  } # mut.calcs for loop end
  
} # mut.sigs for loop end

for(cs in 1:cs.len){

  calc <- calc.sig.combi.df$calc[[cs]]
  sig <- calc.sig.combi.df$sig[[cs]]
  
  
}

# Per sig - Plot mean value at Cp (no ci) vs. Cp

# Plot mean value at Cp (plus ci) vs. Cp
# Per calc, sig, mut types separated because of ci bars

# rm(list=ls()); gc()