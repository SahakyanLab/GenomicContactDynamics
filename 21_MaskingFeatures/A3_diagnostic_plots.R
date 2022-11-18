################################################################################
# Diagnostic plots for assessing masking:
# (Region-wise, unique of bins not taken) : Unmasked length vs. Cp, Total k-mer
# count vs. Cp
# (Contact-wise): sd(unmasked lengths of pairs) vs. Cp
# Calculate pairwise significance
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
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/21_MaskingFeatures")
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
unmskdThresh = 0.5
binkmerMasked.dir = paste0(wk.dir, "/out_binkmer_maskingThreshold/unmskdThresh", unmskdThresh)
out.dir = paste0(wk.dir, "/out_diagnostic_plots")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
chr.id = "chrALL"
bin.len = 40000
kmer.len = 7
affix.binkmer = "_hg19_rm"
Cps = 1:21
yval.mult = 1
plotOnly = T
getPval = F
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareManyDist.R"))
### FUNCTION ###################################################################
plot.viol <- function(df, x.nme, y.nme, ylim.val, plot.title, out.file, ylab){
  
  p <- ggplot(data=df, aes_string(x=x.nme, y=y.nme)) +
    geom_violin(scale="width", alpha=0.5, trim=T, position="dodge", col="gray70", fill="gray50") +
    stat_boxplot(geom="errorbar", width=0.5 / length(levels(df$Cp)), lty="dashed") +
    stat_summary(fun="median") + 
    scale_x_discrete(expand=c(0.05,0)) + 
    scale_y_continuous(limits=ylim.val) + 
    labs(title=plot.title, y=ylab) + 
    bgr2 + 
    theme(aspect.ratio=0.8,
          plot.title=element_text(size=7),
          axis.title.y=element_text(size=10)) 
  
  ggsave(filename=out.file, plot=p, width=10, height=8)
  
  return(p)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(gcb, "_", chr.id, "_", "kmerlen", kmer.len, "_binlen", bin.len, 
                   affix.binkmer, "_maskedfrlessthan", unmskdThresh)

if(!plotOnly){
  
  for(chr in chr.v){
    
    load(file=paste0(binkmerMasked.dir, "/", chr, "_BinKmer", kmer.len, affix.binkmer, ".RData"))
    UMlen <- BINKMER.MX[,4]
    UMlen[ is.na(BINKMER.MX[,5]) ] <- NA
    KCNT <- rowSums(BINKMER.MX[,-(1:4)], na.rm=F)
    rm(BINKMER.MX)
    
    if( !identical(is.na(UMlen), is.na(KCNT)) ){
      rm(UMlen, KCNT)
      stop(paste0(chr, ": Bins considered not identical."))
      break
    } 
    
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
    rm(PERSIST.MX)
    
    bins.mx <- rbind( ij.mx[,c("i", "Cp")], ij.mx[,c("j", "Cp")] ) 
    bins.mx <- cbind( Cp=bins.mx[,"Cp"], KCNT=KCNT[bins.mx[,"i"]], UMlen=UMlen[bins.mx[,"i"]] )
    save(bins.mx, file=paste0(out.dir, "/", chr, "_", out.name, "_binval_vsCp.RData"))
    
    abs.diffs <- abs(apply(X=cbind(UMlen[ij.mx[,"i"]], UMlen[ij.mx[,"j"]]), MARGIN=1, FUN=diff, na.rm=F))
    Sds <- apply(X=cbind(UMlen[ij.mx[,"i"]], UMlen[ij.mx[,"j"]]), MARGIN=1, FUN=sd, na.rm=F)
    ij.mx <- cbind(Cp=ij.mx[,"Cp"], abs.diff.val=abs.diffs, sd.val=Sds)
    save(ij.mx, file=paste0(out.dir, "/", chr, "_", out.name, "_ijval_vsCp.RData"))
    
    rm(UMlen, KCNT, bins.mx, abs.diffs, Sds, ij.mx)
    
  }
  
} 

plot.suff <- paste0("onlyBins>=", unmskdThresh, "unmaskedseq_binswithmissingseqNotIncluded_forRegionwiseplotUniqueBinsNotTaken",
                    "_pointAtMedian_trim=TsoTailUpToRangeOfData_YvalsDividedBy", yval.mult,
                    "_horizontalLineIsDefaultstat_boxplotGeomErrorbar")
plot.title <- paste0(out.name, "\n", plot.suff)

# Bin/Region-wise plot

BINVAL <- list()
for(chr in chr.v){
  
  load(file=paste0(out.dir, "/", chr, "_", out.name, "_binval_vsCp.RData"))
  #file.remove()
  BINVAL[[chr]] <- bins.mx
  rm(bins.mx)
  
}
BINVAL <- do.call("rbind.data.frame", BINVAL)
BINVAL$Cp <- factor(as.character(BINVAL$Cp), levels=as.character(Cps))
BINVAL$KCNT <- BINVAL$KCNT / yval.mult
BINVAL$UMlen <- BINVAL$UMlen / yval.mult
if( identical(is.na(BINVAL$KCNT), is.na(BINVAL$UMlen)) ){
  BINVAL <- na.omit(BINVAL)
} else {
  rm(BINVAL)
  stop("Checkpoint 1.")
}

ylim.val <- range(c(BINVAL$KCNT, BINVAL$UMlen), na.rm=T)

# Plot k-mer count

#pdf(file=paste0(out.dir, "/", out.name, "_KCNT_vsCp_bins.pdf"), width=10, height=8)
#boxplot(formula=KCNT~Cp, data=BINVAL, outline=F, xlab="Cp", ylab="sum(kmer count, na.rm=F)", 
#        cex.main=0.5, main=paste0(out.name, "_", plot.suff), ylim=ylim.val)
#dev.off()

plot.id <- "_KCNT_vsCp_bins"
plot.viol(df=BINVAL, x.nme="Cp", y.nme="KCNT", ylim.val=ylim.val, 
          ylab="sum(kmer count, na.rm=F)", plot.title=plot.title, 
          out.file=paste0(out.dir, "/", out.name, plot.id, ".pdf"))

if(getPval){
  compareManyDist(xval=BINVAL$KCNT, grp=BINVAL$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, plot.id))
}

# Plot unmasked length

plot.id <- "_numUMChar_vsCp_bins"
plot.viol(df=BINVAL, x.nme="Cp", y.nme="UMlen", ylim.val=ylim.val, 
          ylab="numUMChar, unmasked length", plot.title=plot.title, 
          out.file=paste0(out.dir, "/", out.name, plot.id, ".pdf"))

if(getPval){
  compareManyDist(xval=BINVAL$UMlen, grp=BINVAL$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, plot.id))
}

rm(BINVAL, plot.id, ylim.val)
gc()

# Contact-wise plot

IJVAL <- list()
for(chr in chr.v){
  
  load(file=paste0(out.dir, "/", chr, "_", out.name, "_ijval_vsCp.RData"))
  #file.remove()
  IJVAL[[chr]] <- ij.mx
  rm(ij.mx)
  
}

# 
IJVAL <- do.call("rbind.data.frame", IJVAL)
IJVAL$Cp <- factor(as.character(IJVAL$Cp), levels=as.character(Cps))
IJVAL$abs.diff.val <- IJVAL$abs.diff.val / yval.mult
IJVAL$sd.val <- IJVAL$sd.val / yval.mult
if( identical(is.na(IJVAL$KCNT), is.na(IJVAL$UMlen)) ){
  IJVAL <- na.omit(IJVAL)
} else {
  rm(IJVAL)
  stop("Checkpoint 2.")
}

ylim.val <- range(c(IJVAL$abs.diff.val, IJVAL$sd.val), na.rm=T)

# Plot absolute difference of unmasked lengths of two contacting loci

plot.id <- "_absdiffnumUMChar_vsCp_ij"
plot.viol(df=IJVAL, x.nme="Cp", y.nme="abs.diff.val", ylim.val=ylim.val, 
          ylab="abs(diff(unmasked lengths of contact pairs, na.rm=F))", plot.title=plot.title, 
          out.file=paste0(out.dir, "/", out.name, plot.id, ".pdf"))

if(getPval){
  compareManyDist(xval=IJVAL$abs.diff.val, grp=IJVAL$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, plot.id))
}

# Plot sd of unmasked lengths of two contacting loci

plot.id <- "_sdnumUMChar_vsCp_ij"
plot.viol(df=IJVAL, x.nme="Cp", y.nme="sd.val", ylim.val=ylim.val, 
          ylab="sd(unmasked lengths of contact pairs, na.rm=F)", plot.title=plot.title, 
          out.file=paste0(out.dir, "/", out.name, plot.id, ".pdf"))

if(getPval){
  compareManyDist(xval=IJVAL$sd.val, grp=IJVAL$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, plot.id))
}

# rm(list=ls()); gc()