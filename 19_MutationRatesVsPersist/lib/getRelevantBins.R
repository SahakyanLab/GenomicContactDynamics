################################################################################
# Identify bins to be included; bins should not have any missing base, should
# not be the last bin of chr and should have at least one site that overlaps
# with given location (i.e. exon, intron, intergenic, intron + intergenic)
### FUNCTION ###################################################################
getRelevantbins <- function(basecont.dir = 'directory of BINKMER.MX', 
                            basecont.affix = 'affix of BINKMER.MX indicating location',
                            gcb = 'gcb of BINKMER.MX',
                            chr = 'chr of BINKMER.MX'){
  
  load(file=paste0(basecont.dir, "/", chr, "_BinKmer1", basecont.affix, ".RData"))
  a <- rowSums(x=BINKMER.MX[,c("A","C","G","T")], na.rm=F)
  
  # Drop bins not overlapping with given location or those with missing sequence
  drop.TF <- (a==0 & !is.na(a)) | is.na(a)
  
  # By default, exclude last bin in analyses because of <Hi-C bin length
  drop.TF[length(drop.TF)] <- T
  
  return(cbind(bins=BINKMER.MX[,"bins"], relevant=as.numeric(!drop.TF), numWTSEQ=a))
  
}
################################################################################

# rm(list=ls()); gc()