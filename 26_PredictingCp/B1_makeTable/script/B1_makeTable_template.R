################################################################################
# Make ML tabular dataset
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
    CII.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/z_ignore_git/out_constraints/merged_final")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    CII.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_constraints_GfreeSingleNorm/merged_final")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir  = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/26_PredictingCp")
anv.dir = paste0(wk.dir, "/out_getANV")
persist.dir  = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binkmer3.dir = paste0(persist.dir, "/binkmer3_allbins")
binkmer1.dir = paste0(persist.dir, "/out_binBaseContent/maskfile0") 
out.dir = paste0(wk.dir, "/out_makeTable")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrCHRREPLACE"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(Biostrings)
library(compiler)
source(paste0(lib, "/makeKmerStrandInvar.R"))
source(paste0(wk.dir, "/lib/feat_kmer.R"))
source(paste0(wk.dir, "/lib/feat_complement.R"))
source(paste0(wk.dir, "/lib/feat_anv.R"))
source(paste0(wk.dir, "/lib/feat_makeTablePerChr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
TBL <- feat_makeTablePerChr(anv.dir=anv.dir, binkmer3.dir=binkmer3.dir, 
                            binkmer1.dir=binkmer1.dir, CII.dir=CII.dir, 
                            persist.dir=persist.dir, chr=chr, gcb=gcb)

# Remove rows with NA for smaller CSV file
TBL <- TBL[! (is.na(TBL$grp.compl.ij_CIIkmer) | is.na(TBL$grp.compl.ij_CIIalign)),]
if( any(is.na(TBL)) ){
  rm(TBL); stop("Missing values in TBL.")
} else {
  write.csv(TBL, file=paste0(out.dir, "/", gcb, "_", chr, "_MLtbl.csv"), row.names=F)
}

# rm(list=ls()); gc()
