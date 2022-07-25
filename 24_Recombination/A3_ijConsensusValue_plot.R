################################################################################
# Plot contact recombination rates vs. Cp
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
chr.v = paste0("chr", c(1:22, "X"))
Cp.v = 1:21
consensus.method = "MEDIAN.MEDIAN" 
recom.id = "recomRates_2011_01_phaseII_B37_Myers"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(scales) # scientific()
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
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

out.name <- paste0("mthd_", consensus.method, "_", gcb, "_allChr_ijratesVsCp_bp")
save(BOX.STAT, file=paste0(out.dir, "/", out.name, "_stat_", recom.id, ".RData"))
rm(BOX.STAT)

BOX.DF <- stack(BOX.DF)
BOX.DF$ind <- factor(as.character(BOX.DF$ind), levels=as.character(Cp.v))

print(dim(BOX.DF))

num.ij.Cp <- paste(paste("Cp", names(num.ij.Cp), scientific(num.ij.Cp, digits=3), sep="_"), collapse=";")

png(filename=paste0(out.dir, "/", out.name, "_", recom.id, ".png"), height=300*10, width=300*10, res=300)

boxplot(formula=values~ind, data=BOX.DF, outline=F, xlab="Cp", ylab=consensus.method, 
        main=paste0(out.name, "\n", num.ij.Cp), col="#FDC776", cex.main=0.25)

dev.off()

# rm(list=ls()); gc()