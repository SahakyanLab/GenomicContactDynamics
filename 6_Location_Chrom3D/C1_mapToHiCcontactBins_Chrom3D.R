################################################################################
# Use mapToHiCcontactPersistBins() to map DOMXYZR.DF to HiC contact bins.
# deva, R/3.5.0-newgcc, gcc/4.9.2
# Mac, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)
version

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    persist.dir = "/home/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
# DOMXYZR.df
feat.dir = paste0(wk.dir, "/out_AddXYZR")
out.dir = paste0(wk.dir, "/out_mapToHiCcontactBins_Chrom3D/", model.id)
### OTHER SETTINGS #############################################################
ploidy = "haploid"
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 1L #~5G
HiC.res = 4e4L
LOchain = NULL # NULL | "hg19ToHg38"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(itertools)
library(doParallel)
library(rtracklayer)
library(ggplot2)
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOveR.R"))
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOverLoadChain.R"))
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/loadRData.R"))
source(paste0(lib, "/plotLengthDist.R"))
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/mapToHiCcontactPersistBins.R"))
source(paste0(wk.dir, "/lib/LO_mapToHiCcontactPersistBins.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
LO.TF <- ifelse(is.null(LOchain), FALSE, TRUE)
featurefile <- paste0(feat.dir, "/", model.id, "_", ploidy, "_domXYZR.RData")
out.name <- paste0(model.id, "_", ploidy, "_Chrom3DVsPersist") 

mapToHiCcontactPersistBins(
  # PERSIST.MX object directory
  persist.dir=persist.dir,
  
  ## Feature directory
  feat.dir=feat.dir,
  # Separate files per chromosome?
  featSepByChr=FALSE,
  #---- if featSepByChr = TRUE
  # Program will look for each file on feat.dir using chr number eg. *chrX* as pattern
  # Beware of duplicate files for each chr
  
  #---- if featSepByChr = FALSE
  # Path to file
  featurefile=featurefile,
  # Combine output of chromosomes?
  combineOut=TRUE,
  
  feat.header=NULL,
  # Refer to featurefile, column numbers or names 
  start.coord="start",
  end.coord="end",
  # Table coordinate system
  # zero-based (0-start, half-open) or one-based (1-start, fully-closed (1-based))
  # Coordinates stored in UCSC Genome Browser tables are zero-based while those on 
  # web interfaces are one-based
  feat.coordSys="zero-based",
  chr.col="chr",
  # Remove irrelevant columns
  # Recommended that the featurefile should have a column uniquely identifying each entry/row
  drop.col=NULL,
  # Column names after dropping chr.col and drop.col (follow order in original table)
  # If NULL, keep original column names
  coluNames=NULL,
  
  # Output directory
  output.dir=out.dir,
  # Output identifier
  out.name=out.name,
  
  # Hi-C resolution
  bin.len=HiC.res,
  gcb=gcb,
  chr.v=chr.v,
  nCPU=nCPU,
  
  # Overlap parameters
  min.olap=1L,
  max.gap=-1L,
  type.olap="any",
  
  # LiftOver
  doLiftOver=LO.TF,
  LOchain=LOchain
)

# rm(list=ls())


