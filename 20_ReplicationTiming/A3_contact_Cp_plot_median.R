################################################################################
# Generate contact-wise plot of RT vs. Cp. Calculate p-values i.e. ANOVA/KRUSKAL, 
# correlation, pairwise comparison of Cp distributions plus each Cp distribution 
# vs. all long-range contacts (Cp=0). Generate plots and do p-value calculations 
# per rt.type (i.e. all, nontumor, tumor). 
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
src.dir = paste0(wk.dir, "/out_contact_Cp_plotdata")
out.dir = paste0(wk.dir, "/out_contact_Cp_plot")
### OTHER SETTINGS #############################################################
chrs = paste0("chr", c(1:22, "X"))
nCPU = 4
#src.id = "ij.fnxmean_allrttypes_mean_mincelllineWithData59_40000bpHiCres"
src.id = "ij.fnxmean_allrttypes_median_mincelllineWithData59_40000bpHiCres"
Cp.v = 1:21

#rt.type.cols = c(all="#c1cdc1", nontumor="#3288BD", tumor="#e25d6c")
rt.type.cols = c(all="#c1cdc1", nontumor="#4DBBD5FF", tumor="#F39B7FFF")

ci.plot.ylim = c(-0.5, 0.3)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))

# Summary statistics
library(Rmisc)

library(reshape2)
library(gghalves)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))

library(colorspace) # darken colour
rt.type.cols.dark = rgb(hex2RGB(rt.type.cols)@coords * 0.7)
names(rt.type.cols.dark) <- names(rt.type.cols)

library(car) # ANOVA for unbalanced dataset
source(paste0(lib, "/doVarTest_new.R")) # Update deva copy
source(paste0(lib, "/doCorTest.R")) # Update deva copy
source(paste0(lib, "/compareManyDist.R"))  # Update deva copy
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrs.len <- length(chrs)
toExport <- c("chrs", "src.dir", "src.id")

# Combine contact-wise RT data from all chr

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

# Melt df, Convert Cp to factors

df <- as.data.frame(df)
df$Cp <- factor(as.character(df$Cp), levels=as.character(Cp.v))
df <- reshape2::melt(df, id="Cp")
colnames(df) <- c("Cp", "rt.type", "value")

# Remove NAs

df <- na.omit(df)

# Y limits for all distribution plots
y.rng <- range(df$value)
ylim.dist <- round(y.rng, digits=1)
ylim.dist[[1]] <- ylim.dist[[1]] - 0.2
ylim.dist[[2]] <- ylim.dist[[2]] + 0.2

y.label <- paste0("Data_range: ", paste(y.rng, collapse="To"), 
                  "_Plot_ylim: ", paste(ylim.dist, collapse="To"))

# df summary statistics

df.summ.SE <- summarySE(df, measurevar="value", groupvars=c("Cp","rt.type"), na.rm=T, .drop=F)
df.summ.SE.Cp0 <- cbind(Cp="0", summarySE(df, measurevar="value", groupvars="rt.type", na.rm=T, .drop=F))
df.summ.SE <- rbind(df.summ.SE.Cp0, df.summ.SE)

SUMM.STAT <- list(summ.SE=df.summ.SE)
save(SUMM.STAT, file=paste0(out.dir, "/", src.id, "_summary_statistics.RData"))
rm(SUMM.STAT)

# Plot mean + 95% CI

df.summ.SE$Cp <- factor(as.character(df.summ.SE$Cp), levels=as.character(c("0", Cp.v)))
df.summ.SE$rt.type <- factor(as.character(df.summ.SE$rt.type), levels=names(rt.type.cols))

pd <- position_dodge(0.3)
p <- ggplot(data=df.summ.SE, aes(x=Cp, y=value)) +
  geom_errorbar(aes(ymin=value - ci, ymax=value + ci, col=rt.type), 
                width=1, linewidth=0.6, position=pd) + 
  stat_summary(aes(col=rt.type), fun="mean", size=0.35, position=pd) +
  scale_y_continuous(limits=ci.plot.ylim) + 
  scale_colour_manual(values=rt.type.cols.dark[levels(df.summ.SE$rt.type)]) +
  labs(title=paste0(src.id, "_meanPlus95PercCI")) + 
  bgr1

ggsave(filename=paste0(out.dir, "/", src.id, "_meanPlus95PercCI.pdf"), width=10, height=10, 
       plot=p)


# Plot mean + 95% CI (but translated to set Cp=0 means to 0)

# Translate Cp-0 means to 0, translating other Cp values accordingly
df.summ.SE$value.tr <- NA
for( rt.type in levels(df.summ.SE$rt.type) ){
  
  Cp0.mean.val <- df.summ.SE$value[df.summ.SE$Cp == "0" & df.summ.SE$rt.type == rt.type]
  df.summ.SE$value.tr[df.summ.SE$rt.type == rt.type] <- df.summ.SE$value[df.summ.SE$rt.type == rt.type] - Cp0.mean.val
  
}

p <- ggplot(data=df.summ.SE, aes(x=Cp, y=value.tr)) +
  #geom_point(data=aggs, aes(x=Group.2, y=x, aes(col="Group.1"))) + 
  geom_errorbar(aes(ymin=value.tr - ci, ymax=value.tr + ci, col=rt.type), 
                width=1, linewidth=0.6, position=pd) + 
  # stat_summary causes this warning Warning: Removed 66 rows containing missing values (`geom_segment()`).
  # but I've checked that the mean values are right with aggregate() and geom_point()
  stat_summary(aes(col=rt.type), fun="mean", size=0.35, position=pd) +
  scale_y_continuous(limits=ci.plot.ylim) + 
  scale_colour_manual(values=rt.type.cols.dark[levels(df.summ.SE$rt.type)]) +
  labs(title=paste0(src.id, "_meanPlus95PercCI_meanTranslatedtoSetCp0MeanAt0")) + 
  bgr1

ggsave(filename=paste0(out.dir, "/", src.id, "_meanPlus95PercCI_Cp0MeanAt0.pdf"), width=10, height=10, 
       plot=p)

# Boxplots of rt values vs. Cp generated per rt.type

for( rt.type in unique(df$rt.type) ){
  
  df.tmp <- df[df$rt.type == rt.type, ]
  df.tmp <- na.omit(df.tmp)
  bp.stat <- boxplot.stats(x=df.tmp$value)$stats
  
  plot.title <- paste0(rt.type, "_", src.id, "_solidLineMedianALL_dashed1st3rdQuartileAll")
  out.name <- paste0(rt.type, "_", src.id)
  
  p <- ggplot(data=df.tmp, aes(x=Cp, y=value)) +
    geom_hline(yintercept=bp.stat[[3]], col="black") + 
    geom_hline(yintercept=bp.stat[c(2,4)], col="gray50", lty="dashed") + 
    geom_half_violin(side="r", fill=rt.type.cols[[rt.type]], scale="width", lwd=0.6, width=0.8, trim=T) +
    geom_half_boxplot(fill=rt.type.cols.dark[[rt.type]], lwd=0.6, width=0.5, outlier.shape=1) + 
    scale_y_continuous(limits=ylim.dist) + 
    labs(title=paste0(plot.title, "_", y.label)) + 
    bgr1 + 
    theme(plot.title=element_text(size=4)) 
  
  ggsave(filename=paste0(out.dir, "/", out.name, ".png"),
         width=10*300, height=10*300, plot=p, units="px")
  
  rm(p)

  # # P-values, ANOVA/KW and correlation tests
  # 
  # try(doVarTest(xval=df.tmp$value, grp=df.tmp$Cp, out.dir=out.dir, out.name=out.name))
  # 
  # try(doCorTest(xval=as.numeric(as.character(df.tmp$Cp)), yval=df.tmp$value, alt="two.sided",
  #               exactpval=F, out.dir=out.dir, out.name=out.name))
  # 
  # # Treat ALL long-range contact data distribution as extra Cp i.e. Cp = 0 then do
  # # Cp pairwise comparisons
  # 
  # # Add Cp=0
  # df.tmp$Cp <- as.character(df.tmp$Cp)
  # df.tmp.all <- df.tmp
  # df.tmp.all$Cp <- "0"
  # df.tmp <- rbind(df.tmp.all, df.tmp)
  # rm(df.tmp.all)
  # df.tmp$Cp <- factor(as.character(df.tmp$Cp), levels=as.character(c("0", Cp.v)))
  # 
  # try(compareManyDist( xval=df.tmp$value, grp=df.tmp$Cp, alt="two.sided", out.dir=out.dir, 
  #                      out.name=paste0(out.name, "_Cp0To21compare") ))
  
  rm(df.tmp)
  
  message(paste0(rt.type, " done!"))
  
}

## Other p-values

# 1. P-values, pairwise comparisons per rt.type, per Cp 1 to 21

try(compareManyDist( xval=df$value, grp=paste0(df$Cp, df$rt.type), alt="two.sided", out.dir=out.dir, 
                     out.name=paste0(src.id, "_Cp1To21compare") ))

try(compareManyDist( xval=df$value, grp=df$rt.type, alt="two.sided", out.dir=out.dir, 
                     out.name=paste0(src.id, "_Cp0CompareRt.types") ))

# 2. P-values, two-way (Cp and rt.type) ANOVA, Kruskal-Wallis
# Refer to http://www.sthda.com/english/wiki/two-way-anova-test-in-r#compute-two-way-anova-test-in-r-for-unbalanced-designs

TEST <- list()
TEST[["ano2"]] <- aov(value ~ Cp * rt.type, data=df)
TEST[["ano2_unbTyp3"]] <- Anova(TEST[["ano2"]], type="III")
save(TEST, file=paste0(out.dir, "/", src.id, "_ano2_varbasedtest.RData"))

# 3. P-values described below

for( rt.type in unique(df$rt.type) ){
  
  df.tmp <- df[df$rt.type == rt.type, ]
  df.tmp <- na.omit(df.tmp)
  out.name <- paste0(rt.type, "_", src.id)
  
  # P-values, ANOVA/KW and correlation tests
  
  try(doVarTest(xval=df.tmp$value, grp=df.tmp$Cp, out.dir=out.dir, out.name=out.name, lightenAOV=F))
  
  try(doCorTest(xval=as.numeric(as.character(df.tmp$Cp)), yval=df.tmp$value, alt="two.sided",
                exactpval=F, out.dir=out.dir, out.name=out.name))
  
  # Treat ALL long-range contact data distribution as extra Cp i.e. Cp = 0 then do
  # Cp pairwise comparisons
  
  # Add Cp=0
  df.tmp$Cp <- as.character(df.tmp$Cp)
  df.tmp.all <- df.tmp
  df.tmp.all$Cp <- "0"
  df.tmp <- rbind(df.tmp.all, df.tmp)
  rm(df.tmp.all)
  df.tmp$Cp <- factor(as.character(df.tmp$Cp), levels=as.character(c("0", Cp.v)))
  
  try(compareManyDist( xval=df.tmp$value, grp=df.tmp$Cp, alt="two.sided", out.dir=out.dir, 
                       out.name=paste0(out.name, "_Cp0To21compare") ))
  
  rm(df.tmp)
  
  message(paste0(rt.type, " p-value calc done!"))
  
}

# rm(list=ls()); gc()