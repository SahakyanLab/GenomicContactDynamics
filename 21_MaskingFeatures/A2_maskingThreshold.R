################################################################################
# Select the valid chromosome bins based on the fraction of bases masked. 
# Populate with NAs those bins with unmasked fraction less than threshold plus
# bins with any missing sequence in unmasked version of the genome.
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
data.dir = paste0(home.dir, "/Database")
binkmer.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_allbins")
binkmerMasked.dir = paste0(wk.dir, "/out_binkmer_maskedKmerCounts_genomeMasked")
out.dir = paste0(wk.dir, "/out_binkmer_maskingThreshold")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
kmer.len = 7
unmskdThresh = 0.5
affix.binkmer = "_hg19_rm"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
req.bp <- bin.len * unmskdThresh

for(chr in chr.v){
  
  load(file=paste0(binkmer.dir, "/", chr, "_BinKmer", kmer.len, ".RData"))
  UM.MX <- BINKMER.MX
  rm(BINKMER.MX)
  
  load(file=paste0(binkmerMasked.dir, "/", chr, "_BinKmer", kmer.len, affix.binkmer, ".RData"))
  
  if( identical(UM.MX[,1:3], BINKMER.MX[,1:3]) ){
    
    NA.TF <- is.na(UM.MX[,5]) | (BINKMER.MX[,"numUMChar"] < req.bp)
    BINKMER.MX[NA.TF,-(1:4)] <- NA
    save(BINKMER.MX, file=paste0(out.dir, "/", chr, "_BinKmer", kmer.len, affix.binkmer, ".RData"))
    print(chr)
    rm(BINKMER.MX, UM.MX)
    
  } else {
    rm(BINKMER.MX, UM.MX); stop(paste0(chr, ": Wrong order of bins."))
  }
  
}

# rm(list=ls()); gc()