################################################################################
# 
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")

metric = "sumrep"
src.dir = paste0(wk.dir, "/z_ignore_git/out_minRepCounts/subfamALL_", metric, 
                 "_atleast2sumrep/pairwise")
out.dir = paste0(wk.dir, "/z_ignore_git/out_summary_trends")
element.file = paste0(wk.dir, "/Repeat_rankingbyAge/plot_GiorPubl372rankrepFamilies.csv")
### OTHER SETTINGS #############################################################
src.nme = paste0("chrALL_min2Mb_subfamALL_", metric, "Counts")
nCPU = 1 # elements
Cp.v = 1:21
# Cp MEAN or MED (median)?
average.fnx = "MED"
fdr.cutoff = 0.05

age.group.num = 1 # Choose number based on timepoints in the old density map
#age.numPer.group = 372 #7 / age.group.num # Should be integer

elm.file = paste0(wk.dir, "/z_ignore_git/out_combine/out_filterElmTRUE_repName.txt")
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
Cp.v.len <- length(Cp.v)
#elements <- read.csv(element.file, header=T, stringsAsFactors=F)$repName
elements <- readLines(elm.file) #c("AluY", "AluYa5", "L2a", "MIR", "MIR3", "MIRb", "MIRc")
elm.len <- length(elements)
age.numPer.group <- elm.len / age.group.num

# Data for plotting

#### PARALLEL EXECUTION #########
df <- foreach(itr=isplitVector(1:elm.len, chunks=nCPU), .inorder=T, .combine="rbind"

) %op% {
  
  chunk <- sapply(itr, simplify=F, FUN=function(elm.ind){
    
    elm <- elements[[elm.ind]]
    elm <- gsub(pattern="[^[:alnum:][:space:]]", replacement="", x=elm)
   
    # Load pairwisedifftest.RData
    #fle.path <- paste0(src.dir, "/", src.nme, "_", elm, "_", elm.ind, "_pairwisedifftest.RData")
    fle.path <- list.files(path=src.dir, pattern=paste0(src.nme, "_", elm, "_"), full.names=T)
    
    if( length(fle.path) == 1 ){
      load(fle.path)
    } else {
      print(paste0(elm.ind, " ", elm, ": pairwisedifftest.RData does not exist. Returning NULL."))
      rm(elm)
      return(NULL)
    }
    pval.mx <- TEST$pmw$p.value
    pval.mx[is.nan(pval.mx)] <- NA
    
    # Summary of trends looking at if value at a Cp greater/less/ns that value of Cp - 1

    trend.mx <- matrix(data=NA, nrow=length(Cp.v) + 1, ncol=4,
                       dimnames=list(c("0", Cp.v), c("Cpafter", "BH", "ave.val", "ave.CpMINUSPrevCp")))
    
    # BH adjusted p-values
    BH.df <- reshape2::melt(t(pval.mx)) # transpose so Var1 is Cpbefore, Var2 is Cpafter
    BH.df <- BH.df[!is.na(BH.df$value),]
    colnames(BH.df) <- c("Cpbefore", "Cpafter", "BH")
    # Get pairs of consecutive Cps
    BH.df <- BH.df[ (as.numeric(BH.df$Cpafter) - as.numeric(BH.df$Cpbefore)) == 1, ]
    
    # Fill out missing values
    BH.mx <- matrix(data=NA, nrow=Cp.v.len, ncol=3, dimnames=list(NULL, colnames(BH.df)))
    BH.mx[,"Cpafter"] <- Cp.v
    BH.mx[,"Cpbefore"] <- BH.mx[,"Cpafter"] - 1
    BH.mx[BH.df$Cpafter,"BH"] <- BH.df$BH
    BH.df <- as.data.frame(BH.mx)
      
    ave.df <- as.data.frame(do.call("rbind", TEST$meanmed))
    ave.mx <- matrix(data=NA, nrow=Cp.v.len, ncol=ncol(ave.df), dimnames=list(as.character(Cp.v), colnames(ave.df)))
    for( col.nme in colnames(ave.df) ){
      ave.mx[rownames(ave.df),col.nme] <- ave.df[[col.nme]]
    }
    ave.df <- as.data.frame(ave.df)
    
    rm(BH.mx, ave.mx)
    
    #
    if( !identical(as.character(BH.df$Cpafter), as.character(Cp.v)) ){
      rm(trend.mx)
      stop(paste0(elm, ": Error in getting BH adjusted p-valus data."))
    } 
    
    # Add other data to trend.mx
    
    Cpafter <- c("0", Cp.v)
    trend.mx[,"Cpafter"] <- as.numeric(Cpafter)
    trend.mx[-1,"BH"] <- BH.df$BH
    trend.mx[,"ave.val"] <- ave.df[Cpafter,average.fnx]
    trend.mx[-1,"ave.CpMINUSPrevCp"] <- unname(diff(trend.mx[,"ave.val"]))
    
    trend.df <- as.data.frame(trend.mx)
    trend.df <- cbind(type=ceiling(elm.ind / age.numPer.group), trend.df)
    
    print(elm.ind, quote=F)
    
    if(nrow(trend.df) != 22 ){
      stop(paste0(elm, ": trend.df not 22 rows."))
    }
    
    trend.df <- cbind(elm.ind=elm.ind, trend.df)
    return(trend.df)
    
  })
  
  return(do.call("rbind", chunk))
  
}
### END OF PARALLEL EXECUTION ###

## PLOT

# Alter df for plotting
 
df <- df[! df$Cpafter %in% c(0,1),]
df$Cpafter <- factor(as.character(df$Cpafter), levels=as.character(Cp.v))
df$BH[!is.finite(df$BH)] <- 1.5

if(! identical(df, na.omit(df)) ){
  warning("df has NA/s. Removed")
} else {
  df <- na.omit(df)
}

# Add other pertinent columns

# gr (Cp ave greater than ave of previous Cp), le (less), eq (equal)
df$trend.CpMINUSPrevCp <- df$ave.CpMINUSPrevCp
df$trend.CpMINUSPrevCp[ df$ave.CpMINUSPrevCp > 0 ] <- "gr"
df$trend.CpMINUSPrevCp[ df$ave.CpMINUSPrevCp < 0 ] <- "le"
df$trend.CpMINUSPrevCp[ df$ave.CpMINUSPrevCp == 0 ] <- "eq"

#  

df$significance <- df$BH
df$significance[ df$BH < fdr.cutoff & !is.na(df$BH) ] <- "s"
df$significance[ df$BH >= fdr.cutoff & !is.na(df$BH) ] <- "ns"
df$significance[ df$BH > 1 ] <- "nm" # not measured deal with this
df$trend.group <- paste0(df$trend.CpMINUSPrevCp, ".", df$significance)
df$trend.group <- factor(as.character(df$trend.group),  
                         levels=c("gr.s", "gr.ns", "eq.s", 
                                  "eq.ns", "le.s", "le.ns",
                                  "gr.nm", "eq.nm", "le.nm"))


age.id <- paste0("fdrcutoff", fdr.cutoff, "_Cp", average.fnx, "_age.group.num", age.group.num)
out.id.general <- paste0(src.nme, "_", age.id)
save(df, file=paste0(out.dir, "/", out.id.general, "plot.RData"))

# rm(list=ls()); gc()

## PLOT METRICS TOGETHER

# Variables needed for plots

out.id.fin <- paste0(out.id.general, "_barplot")

cols.trend.group <- setNames(object=c("#E64B35FF", "#F39B7FFF", "#7E6148FF", 
                                                 "#B09C85FF", "#3C5488FF", "#8491B4FF",
                                                 "#DC0000FF", "#7E6148FF", "#4DBBD5FF"),
                                                 nm=levels(df$trend.group))
pd <- position_dodge(0.1)

# Plot, type vs. loc panel, combine sig, calc

src.nmes <- c(paste0("chrALL_min2Mb_subfamALL_sumrepCounts"),
              paste0("chrALL_min2Mb_subfamALL_skewrepCounts"),
              paste0("chrALL_min2Mb_subfamALL_minrepCounts")
              )
#src.nmes <- c(#paste0("chrALL_min2Mb_subfamALL_sumrepCounts"),
#              paste0("chrALL_min2Mb_subfamALL_skewrepCounts")#,
#              #paste0("chrALL_min2Mb_subfamALL_minrepCounts")
#)

P.LST <- sapply(src.nmes, simplify=F, FUN=function(src.nme){
  
  load(paste0(out.dir, "/", src.nme, "_", age.id, "plot.RData"))
  
  trend.group.present <- levels(df$trend.group)[levels(df$trend.group) %in% df$trend.group]
  plot.title <- paste0(src.nme, "_onlyShownAreTrend.groupPresent")
  
  # p <- ggplot(df, aes(Cpafter)) +
  #   geom_bar(aes(fill=trend.group)) +
  #   #scale_y_continuous(breaks=seq(0,age.numPer.group,10)) + 
  #   scale_fill_manual(limits=trend.group.present,
  #                     values=cols.trend.group[trend.group.present]) +
  #   labs(title=plot.title) + 
  #   bgr1 +
  #   facet_grid(type~.) +
  #   theme(legend.title=element_text(size=5), legend.text=element_text(size=5),
  #         plot.title=element_text(size=5))
  
  #df$elm.ind <- as.character(df$elm.ind)
  df$elm <- elements[df$elm.ind]
  df$elm <- factor(df$elm, levels=rev(elements))
  
  p <- ggplot(df, aes(x=Cpafter, y=elm)) +
    geom_tile(aes(fill=trend.group), colour="white") + 
    #geom_bar(aes(fill=trend.group)) +
    #scale_y_continuous(breaks=seq(0,age.numPer.group,10)) + 
    scale_fill_manual(limits=trend.group.present,
                      values=cols.trend.group[trend.group.present]) +
    labs(title=plot.title) + 
    bgr1 +
    facet_grid(type~.) +
    theme(legend.title=element_text(size=5), legend.text=element_text(size=5),
          plot.title=element_text(size=5))
  
  return(p)
  
})


# Version with all details
row.num <- 1
col.num <- length(src.nmes)
p.arr <- ggarrange(plotlist=P.LST, nrow=1, ncol=3, common.legend=T, legend="top")
ggexport(p.arr, width=col.num * 5, height=row.num * age.group.num * 5, 
         filename=paste0(out.dir, "/calcCombined_", out.id.fin, "_withLabels.pdf"))

# Final minimal version
P.LST <- lapply(P.LST, FUN=function(p){
  
  p <- p + 
    labs(title=NULL) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.title.y=element_blank(), 
          legend.position="none", strip.text=element_blank()) 
  return(p)
  
})

p.arr <- plot_grid(plotlist=P.LST, nrow=1, ncol=3, align="hv", byrow=T)
save_plot(filename=paste0(out.dir, "/calcCombined_", out.id.fin, "_noLabels.pdf"), plot=p.arr,
          base_width=col.num * 5, base_height=row.num * age.group.num * 5)

# rm(list=ls()); gc()



