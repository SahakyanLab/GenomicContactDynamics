################################################################################
# Plot contact recombination rates vs. Cp. Calculate all p-values needed.
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
    home.dir = "/project/sahakyanlab/ltamon" 
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/24_Recombination")
ijrates.dir = paste0(wk.dir, "/out_ijConsensusValue/tmp")
out.dir = paste0(wk.dir, "/out_ijConsensusValue_plot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(21:22))
Cp.v = 1:21
consensus.method = "MEAN.MEAN" #"MEDIAN.MEDIAN" 
src.id = "min.countPerBin3_ij_recomRates_2011_01_phaseII_B37_Myers"
col.hex = "#c1cdc1"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(scales) # scientific()
library(Rmisc)  # summarySE()

library(ggplot2)
library(gghalves)
source(paste0(lib, "/GG_bgr.R"))
library(colorspace)

library(car) # ANOVA for unbalanced dataset
source(paste0(lib, "/doVarTest.R")) # Update deva copy
source(paste0(lib, "/doCorTest.R")) # Update deva copy
source(paste0(lib, "/compareManyDist.R"))  # Update deva copy
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Plot data

BOX.STAT <- matrix(data=NA, nrow=length(Cp.v) + 1, ncol=5 + 3, 
                   dimnames=list(c("0", Cp.v), c("low.whisk", "low.hinge", "MEDIAN", "upp.hinge", "upp.whisk", 
                                         "MEAN", "MIN", "MAX") )) 

BOX.DF <- list() 
# Count contacts with value per Cp
num.ij.Cp <- setNames(rep(NA, times=length(Cp.v)), nm=Cp.v)
  
for(Cp in Cp.v){
  
  IJ.RATES.CP <- sapply(X=chr.v, simplify=T, FUN=function(chr){
    load(paste0(ijrates.dir, "/mthd_", consensus.method, "_Cp", Cp, "_", chr, "_", gcb, "_", src.id, ".RData"))
    return(IJ.RATES)
  })
  IJ.RATES.CP <- unname(unlist(IJ.RATES.CP))
  
  # Match order with dimnames
  BOX.STAT[as.character(Cp),] <- c( boxplot.stats(x=IJ.RATES.CP, coef=1.5, do.conf=F, do.out=F)$stats, # "NAs and NaNs are allowed and omitted"
                                    mean(IJ.RATES.CP, na.rm=T), min(IJ.RATES.CP, na.rm=T), max(IJ.RATES.CP, na.rm=T) ) 
  
  BOX.DF[[as.character(Cp)]] <- IJ.RATES.CP
  
  num.ij.Cp[as.character(Cp)] <- sum(!is.na(IJ.RATES.CP), na.rm=T)
    
  rm(IJ.RATES.CP)
  
  print(paste0(Cp, ": Chr data collected."))
     
}

BOX.DF <- stack(BOX.DF)
BOX.DF$ind <- factor(as.character(BOX.DF$ind), levels=as.character(Cp.v))
setnames(BOX.DF, old="ind", new="Cp")
print(dim(BOX.DF))

# Remove NAs

BOX.DF <- na.omit(BOX.DF)

#

out.name <- paste0("mthd_", consensus.method, "_", gcb, "_allChr_ijratesVsCp_", src.id)

# Boxplot statistics

# Stats for Cp=0 (all long-range contacts)
BOX.STAT["0",] <- c( boxplot.stats(x=BOX.DF$values, coef=1.5, do.conf=F, do.out=F)$stats, # "NAs and NaNs are allowed and omitted"
                     mean(BOX.DF$values, na.rm=T), min(BOX.DF$values, na.rm=T), max(BOX.DF$values, na.rm=T) ) 

save(BOX.STAT, file=paste0(out.dir, "/", out.name, "_boxplot_stats.RData"))
#rm(BOX.STAT)

# Summary statistics

SUMM.STAT <- data.frame(Cp="0", summarySE(BOX.DF, measurevar="values", na.rm=T, .drop=F)[,-1])
SUMM.STAT <- rbind(SUMM.STAT,
                   summarySE(BOX.DF, measurevar="values", groupvars="Cp", na.rm=T, .drop=F))
save(SUMM.STAT, file=paste0(out.dir, "/", out.name, "_summstat.RData"))

# Boxplor no outlier

BOX.STAT <- as.data.frame(BOX.STAT, stringsAsFactors=F)
BOX.STAT$Cp <- factor(rownames(BOX.STAT), levels=rownames(BOX.STAT))

bp.stat <- boxplot.stats(x=BOX.DF$values)$stats
col.dark <- rgb(hex2RGB(col.hex)@coords * 0.7)

#ggplot(BOX.STAT) +
p <- ggplot(BOX.STAT[BOX.STAT$Cp != "0",]) +
  geom_hline(yintercept=bp.stat[[3]], col="black") + 
  geom_hline(yintercept=bp.stat[c(2,4)], col="gray50", lty="dashed") +
  geom_errorbar(aes(x=Cp, ymin=low.whisk, ymax=upp.whisk), width=0.3, linewidth=0.6) + 
  geom_boxplot(stat="identity", fill=col.dark, lwd=0.6, width=0.5, 
               aes(x=Cp, ymin=low.whisk, ymax=upp.whisk, lower=low.hinge, middle=MEDIAN, upper=upp.hinge)) +
  labs(y="values", title=paste0(out.name, "_NoOutlier")) + 
  bgr1 + 
  theme(plot.title=element_text(size=7)) 

ggsave(filename=paste0(out.dir, "/", out.name, "_NoOutlierbp.png"), 
       width=10*300, height=10*300, plot=p, units="px")

# Boxplot no outlier (truncated), half boxplot half violin

p <- ggplot(data=BOX.DF, aes(x=Cp, y=values)) +
  geom_hline(yintercept=bp.stat[[3]], col="black") + 
  geom_hline(yintercept=bp.stat[c(2,4)], col="gray50", lty="dashed") + 
  geom_half_violin(side="r", fill=col.hex, scale="width", lwd=0.6, width=0.8, trim=T) +
  geom_half_boxplot(fill=col.dark, lwd=0.6, width=0.5, outlier.shape=1) + 
  #geom_errorbar(data=BOX.STAT[BOX.STAT$Cp != "0",], inherit.aes=F,
  #              aes(x=Cp, ymin=low.whisk, ymax=upp.whisk), width=0.3, linewidth=0.6) + 
  #geom_boxplot(fill=col.dark, lwd=0.6, width=0.5, outlier.shape=NA) + 
  labs(title=paste0(out.name, "_withOutlier")) + 
  bgr1 + 
  theme(plot.title=element_text(size=7)) +
  coord_cartesian(ylim=c(0,2.5))

ggsave(filename=paste0(out.dir, "/", out.name, "_halfbp_halfviol_truncbp.png"), 
       width=10*300, height=10*300, plot=p, units="px")

# Full boxplot with outlier

p <- ggplot(data=BOX.DF, aes(x=Cp, y=values)) +
  geom_hline(yintercept=bp.stat[[3]], col="black") + 
  geom_hline(yintercept=bp.stat[c(2,4)], col="gray50", lty="dashed") + 
  geom_half_violin(side="r", fill=col.hex, scale="width", lwd=0.6, width=0.8, trim=T) +
  geom_half_boxplot(fill=col.dark, lwd=0.6, width=0.5, outlier.shape=1) + 
  #geom_errorbar(data=BOX.STAT[BOX.STAT$Cp != "0",], inherit.aes=F,
  #              aes(x=Cp, ymin=low.whisk, ymax=upp.whisk), width=0.3, linewidth=0.6) + 
  #geom_boxplot(fill=col.dark, lwd=0.6, width=0.5, outlier.shape=NA) + 
  labs(title=paste0(out.name, "_withOutlier")) + 
  bgr1 + 
  theme(plot.title=element_text(size=7)) 

ggsave(filename=paste0(out.dir, "/", out.name, "_halfbp_halfviol_fullbp.png"), 
       width=10*300, height=10*300, plot=p, units="px")

# Plot mean + 95% CI

SUMM.STAT$Cp <- factor(as.character(SUMM.STAT$Cp), levels=as.character(c("0", Cp.v)))
SUMM.STAT$dummyleg <- "all"

col.dark1 <- rgb(hex2RGB(col.hex)@coords * 0.7)

pd <- position_dodge(0.3)
p <- ggplot(data=SUMM.STAT, aes(x=Cp, y=values)) +
  geom_errorbar(aes(col=dummyleg, ymin=values - ci, ymax=values + ci), width=0.4, linewidth=0.6, 
                position=pd) +
  stat_summary(aes(col=dummyleg), fun="mean", size=0.35, position=pd) +
  scale_colour_manual(values=col.dark1) +
  labs(title=paste0(out.name, "_meanPlus95PercCI")) + 
  bgr1

ggsave(filename=paste0(out.dir, "/", out.name, "_meanPlus95PercCI.pdf"), 
       width=10, height=10, plot=p)

## P-values

# P-values, ANOVA/KW and correlation tests

try(doVarTest(xval=BOX.DF$values, grp=BOX.DF$Cp, out.dir=out.dir, out.name=out.name))

try(doCorTest(xval=as.numeric(as.character(BOX.DF$Cp)), yval=BOX.DF$values, alt="two.sided",
              exactpval=F, out.dir=out.dir, out.name=out.name))

# Add all values as Cp=0

BOX.DF.Cp0 <- BOX.DF
BOX.DF.Cp0$Cp <- factor("0", levels="0")
BOX.DF <- rbind(BOX.DF.Cp0, BOX.DF)

try(compareManyDist( xval=BOX.DF$values, grp=BOX.DF$Cp, alt="two.sided", out.dir=out.dir, 
                     out.name=paste0(out.name, "_Cp0To21compare") ))

# rm(list=ls()); gc()