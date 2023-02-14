################################################################################
# Plot contact recombination rates vs. Cp. Calculate all p-values needed.
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
    home.dir = "/project/sahakyanlab/ltamon" 
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/24_Recombination")
ijrates.dir = paste0(wk.dir, "/z_ignore_git/out_ijConsensusValue/tmp")
out.dir = paste0(wk.dir, "/out_ijConsensusValue_plot1")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(21:22))
Cp.v = 1:21
consensus.method = "MEDIAN.MEDIAN" 
recom.id = "recomRates_2011_01_phaseII_B37_Myers" #"min.countPerBin3_ij_recomRates_2011_01_phaseII_B37_Myers" 
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
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Plot data

BOX.STAT <- matrix(data=NA, nrow=length(Cp.v), ncol=5 + 3, 
                   dimnames=list(Cp.v, c("low.whisk", "low.hinge", "MEDIAN", "upp.hinge", "upp.whisk", 
                                         "MEAN", "MIN", "MAX") )) 

BOX.DF <- list() 
# Count contacts with value per Cp
num.ij.Cp <- setNames(rep(NA, times=length(Cp.v)), nm=Cp.v)
  
for(Cp in Cp.v){
  
  IJ.RATES.CP <- sapply(X=chr.v, simplify=T, FUN=function(chr){
    load(paste0(ijrates.dir, "/mthd_", consensus.method, "_Cp", Cp, "_", chr, "_", gcb, "_ij_", recom.id, ".RData"))
    return(IJ.RATES)
  })
  IJ.RATES.CP <- unname(unlist(IJ.RATES.CP))
  
  # Match order with dimnames
  BOX.STAT[Cp,] <- c( boxplot.stats(x=IJ.RATES.CP, coef=1.5, do.conf=F, do.out=F)$stats, # "NAs and NaNs are allowed and omitted"
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

# 

out.name <- paste0("mthd_", consensus.method, "_", gcb, "_allChr_ijratesVsCp_", recom.id)

# BOXPLOT statistics

save(BOX.STAT, file=paste0(out.dir, "/", out.name, "_boxplot_stats.RData"))
rm(BOX.STAT)

# Summary statistics

SUMM.STAT <- data.frame(Cp="0", summarySE(BOX.DF, measurevar="values", na.rm=T, .drop=F)[,-1])
SUMM.STAT <- rbind(SUMM.STAT,
                   summarySE(BOX.DF, measurevar="values", groupvars="Cp", na.rm=T, .drop=F))
save(SUMM.STAT, file=paste0(out.dir, "/", out.name, "_summstat.RData"))

# Boxplot

bp.stat <- boxplot.stats(x=BOX.DF$values)$stats

p <- ggplot(data=BOX.DF, aes(x=Cp, y=values)) +
  geom_hline(yintercept=bp.stat[[3]], col="black") + 
  geom_hline(yintercept=bp.stat[c(2,4)], col="gray50", lty="dashed") + 
  geom_half_violin(side="r", fill="#c1cdc1", scale="width", lwd=0.6, 
                   width=0.8, trim=T) +
  geom_half_boxplot(fill=rgb(hex2RGB("#c1cdc1")@coords * 0.8), 
                    lwd=0.6, width=0.5, outlier.shape=1) + 
  labs(title=paste0(out.name, "_withOutlier")) + 
  bgr1 + 
  theme(plot.title=element_text(size=7)) 

ggsave(filename=paste0(out.dir, "/", out.name, ".png"),
       width=10*300, height=10*300, plot=p, units="px")


#num.ij.Cp <- paste(paste("Cp", names(num.ij.Cp), scientific(num.ij.Cp, digits=3), sep="_"), collapse=";")

png(filename=paste0(out.dir, "/", out.name, "_bp.png"), height=300*10, width=300*10, res=300)

plot(x=BOX.DF$Cp, y=BOX.DF$values, panel.first=c(abline(h=bp.stat[[3]], col="black", lty=1, lwd=2)))

plot(BOX.DF$Cp, BOX.DF$values, type = "n")
grid()
points(x, y, col = 'red', type = 'o', lwd = 3, pch = 15)

plot(0,type='n',axes=FALSE,ann=FALSE)
abline(h=bp.stat[[3]], col="black", lty=1, lwd=2)
abline(h=bp.stat[c(2,4)], col="gray50", lty=2, lwd=2)
boxplot(formula=values~Cp, data=BOX.DF, outline=F, xlab="Cp", ylab=consensus.method, 
        main=paste0(out.name, "_noOutlier"), col="#c1cdc1", cex.main=0.25, add=T)



dev.off()

# Mean + ci

# rm(list=ls()); gc()