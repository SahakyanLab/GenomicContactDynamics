################################################################################
# Generate contact-wise plot of RT vs. Cp. Calculate p-values i.e. ANOVA/KRUSKAL, 
# correlation, pairwise comparison of Cp distributions plus each Cp distribution 
# vs. all long-range contacts. Generate plots and do p-value calculations per 
# rt.type (i.e. all, nontumor, tumor). 
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
src.dir = paste0(wk.dir, "/out_contact_Cp_plotdata")
out.dir = paste0(wk.dir, "/out_contact_Cp_plot")
### OTHER SETTINGS #############################################################
chrs = paste0("chr", c(21:22))
nCPU = 1
src.id = "ij.fnxmean_allrttypes_mean_normNot1_setpointcountGrEq3_40000bpHiCres"
Cp.v = 1:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))

library(reshape2)
library(gghalves)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))

source(paste0(lib, "/doVarTest.R"))
source(paste0(lib, "/doCorTest.R")) 
source(paste0(lib, "/compareManyDist.R"))  # Update deva copy
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrs.len <- length(chrs)
toExport <- c("chrs", "src.dir", "src.id")
#### PARALLEL EXECUTION #########
df <- foreach(itr=isplitVector(1:chrs.len, chunks=nCPU), .combine="rbind", .inorder=F,
              .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    chr <- chrs[[i]]
    load(paste0(src.dir, "/", chr, "_", src.id, ".RData"))
    
    return(IJ.RT.TYPES)
    
  }
  
}
### END OF PARALLEL EXECUTION ###

# Transform df for plotting

df <- as.data.frame(df)
df$Cp <- factor(as.character(df$Cp), levels=as.character(Cp.v))
df <- reshape2::melt(df, id="Cp")
colnames(df) <- c("Cp", "rt.type", "value")

# Plot

rt.type.cols <- c(all="honeydew3", nontumor="#3288BD", tumor="#e25d6c")
for( rt.type in unique(df$rt.type) ){
  
  df.tmp <- df[df$rt.type == rt.type, ]
  df.tmp <- na.omit(df.tmp)
  bp.stat <- boxplot.stats(x=df.tmp$value)$stats
  
  plot.title <- paste0(rt.type, "_", src.id, "_solidLineMedianALL_dashed1st3rdQuartileAll")
  out.name <- paste0(rt.type, "_", src.id)
  
  p <- ggplot(data=df.tmp, aes(x=Cp, y=value)) +
    geom_hline(yintercept=bp.stat[[3]], col="black") + 
    geom_hline(yintercept=bp.stat[c(2,4)], col="gray50", lty="dashed") + 
    geom_half_violin(side="r", fill=rt.type.cols[[rt.type]], scale="width", lwd=0.6, 
                     width=0.8, trim=T) +
    geom_half_boxplot(fill=rt.type.cols[[rt.type]], lwd=0.6, width=0.5, outlier.shape=1) + 
    scale_y_continuous(limits=c(-2,2)) + 
    labs(title=plot.title) + 
    bgr1 + 
    theme(plot.title=element_text(size=7)) 
  
  ggsave(filename=paste0(out.dir, "/", out.name, ".png"),
         width=10*300, height=10*300, plot=p, units="px")
  
  message(paste0(rt.type, " done!"))
  
  # P-values
  
  try(doVarTest(xval=df.tmp$value, grp=df.tmp$Cp, out.dir=out.dir, out.name=out.name))
  
  try(doCorTest(xval=as.numeric(as.character(df.tmp$Cp)), yval=df.tmp$value, alt="two.sided", 
                exactpval=F, out.dir=out.dir, out.name=out.name))
  
  # All Cp=0 (all contacts) to data for ALL vs. each Cp p-value calculation
  
  df.tmp$Cp <- as.character(df.tmp$Cp)
  df.tmp.all <- df.tmp
  df.tmp.all$Cp <- "0"
  df.tmp <- rbind(df.tmp, df.tmp.all)
  df.tmp$Cp <- factor(as.character(df.tmp$Cp), levels=as.character(c("0", Cp.v)))
  
  try(compareManyDist(xval=df.tmp$value, grp=df.tmp$Cp, alt="two.sided", out.dir=out.dir, out.name=out.name))
  
}

# rm(list=ls()); gc()