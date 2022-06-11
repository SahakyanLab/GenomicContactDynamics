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
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/21_MaskingFeatures")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/12_MaskingFeatures")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binkmerMasked.dir = paste0(wk.dir, "/out_binkmer_maskingThreshold")
out.dir = paste0(wk.dir, "/out_diagnostic_plots")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
chr.id = "chrALL"
bin.len = 40000
kmer.len = 7
affix.binkmer = "_hg19_rm"
generateSRC = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
source(paste0(lib, "/compareManyDist.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(generateSRC){
  
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
    save(bins.mx, file=paste0(out.dir, "/", gcb, "_", chr, "_KCNT_numUMChar_vsCp_bins.RData"))
    
    Sds <- apply(X=cbind(UMlen[ij.mx[,"i"]], UMlen[ij.mx[,"j"]]), MARGIN=1, FUN=sd, na.rm=F)
    ij.mx <- cbind(Cp=ij.mx[,"Cp"], sd.val=Sds)
    save(ij.mx, file=paste0(out.dir, "/", gcb, "_", chr, "_sdnumUMChar_vsCp_ij.RData"))
    
    rm(UMlen, KCNT, bins.mx, Sds, ij.mx)
    
  }
  
} else {
  
  out.name <- paste0(gcb, "_", chr.id, "_", "kmerlen", kmer.len, "_binlen", bin.len, "_maskedfrlessthan0.5")
  suff.name <- "onlyBins>=0.5unmaskedseq_binswithmissingseqNotIncluded"
  
  # Bin/Region-wise plot
  
  CNTLEN <- list()
  for(chr in chr.v){
    
    load(file=paste0(out.dir, "/", gcb, "_", chr, "_KCNT_numUMChar_vsCp_bins.RData"))
    #file.remove()
    CNTLEN[[chr]] <- bins.mx
    rm(bins.mx)
    
  }
  CNTLEN <- do.call("rbind.data.frame", CNTLEN)
  CNTLEN$Cp <- factor(as.character(CNTLEN$Cp), levels=as.character(1:21))
  CNTLEN <- na.omit(CNTLEN)
  
  pdf(file=paste0(out.dir, "/", out.name, "_KCNT_vsCp_bins.pdf"), width=10, height=8)
  boxplot(formula=KCNT~Cp, data=CNTLEN, outline=F, xlab="Cp", ylab="sum(kmer count, na.rm=F)", 
          cex.main=0.5, main=paste0(out.name, "_", suff.name, "_uniqueBinsNotTaken"))
  dev.off()
  
  compareManyDist(xval=CNTLEN$KCNT, grp=CNTLEN$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, "_KCNT_vsCp_bins"))  
  
  #
  pdf(file=paste0(out.dir, "/", out.name, "_numUMChar_vsCp_bins.pdf"), width=10, height=8)
  boxplot(formula=UMlen~Cp, data=CNTLEN, outline=F, xlab="Cp", ylab="numUMChar, unmasked length", 
          cex.main=0.5, main=paste0(out.name, "_", suff.name))
  dev.off()
  
  compareManyDist(xval=CNTLEN$UMlen, grp=CNTLEN$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, "_numUMChar_vsCp_bins"))
  
  # Contact-wise plot
  
  SD <- list()
  for(chr in chr.v){
    
    load(file=paste0(out.dir, "/", gcb, "_", chr, "_sdnumUMChar_vsCp_ij.RData"))
    #file.remove()
    SD[[chr]] <- ij.mx
    rm(ij.mx)
    
  }
  
  SD <- do.call("rbind.data.frame", SD)
  SD$Cp <- factor(as.character(SD$Cp), levels=as.character(1:21))
  SD <- na.omit(SD)
  
  pdf(file=paste0(out.dir, "/", out.name, "_sdnumUMChar_vsCp_ij.pdf"), width=10, height=8)
  boxplot(formula=sd.val~Cp, data=SD, outline=F, xlab="Cp", 
          ylab="sd(unmasked lengths of contact pairs, na.rm=F)", 
          cex.main=0.5, main=paste0(out.name, "_", suff.name))
  dev.off()
  
  compareManyDist(xval=SD$sd.val, grp=SD$Cp, alt="two.sided", out.dir=out.dir, 
                  out.name=paste0(out.name, "_sdnumUMChar_vsCp_ij"))
  
}

# rm(list=ls()); gc()