################################################################################
# Map recombination rate data to Hi-C contact persist bins
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
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
feat.dir = paste0(data.dir, "/recomRates_2011-01_phaseII_B37_Myers/BED")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/24_Recombination")
out.dir = paste0(wk.dir, "/out_mapToHiCcontactPersistBins")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
library(itertools)
library(doParallel)
library(GenomicRanges)
library(compiler)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/loadRData.R"))
source(paste0(lib, "/mapToHiCcontactPersistBins.R")) 
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
mapToHiCcontactPersistBins(
  
  # PERSIST.MX object directory
  persist.dir = persist.dir,
  
  ## Feature directory
  feat.dir = feat.dir,
  # Separate files per chromosome?
  featSepByChr = T,
  #---- if featSepByChr = TRUE
  # Program will look for each file on feat.dir using chr number and additional affixes 
  # eg. *pattern.pref + chrX + pattern.suff* as pattern
  # Beware of duplicate files for each chr
  pattern.pref = "GRCh37_",
  pattern.suff = ".bed",
  
  #---- if featSepByChr = FALSE
  # Path to file
  featurefile = NULL, #paste0(feat.dir, "/hg19_TOP2B_MCF7_GSE66753_atleast_two_peaks.bed"),
  # Combine output of chromosomes?
  combineOut = F,
  
  feat.header = F,
  # Column numbers or names
  start.coord = 2,
  end.coord = 3,
  # Table coordinate system
  # zero-based (0-start, half-open) or one-based (1-start, fully-closed (1-based))
  # Coordinates stored in UCSC Genome Browser tables are zero-based while those on 
  # web interfaces are one-based
  feat.coordSys = "zero-based",
  chr.col = 1,
  # Remove irrelevant columns
  # Recommended that the featurefile should have a column uniquely identifying each entry/row
  drop.col = NULL, #c("species", "species.count", "score", "r", "g", "b"), #NULL, 
  # Column names after dropping chr.col and drop.col (follow order in original table)
  # If NULL, keep original column names
  coluNames = NULL,
  
  # Output directory
  out.dir = out.dir,
  # Output identifier
  out.name = "recomRates_2011_01_phaseII_B37_Myers",
  
  # HiC resolution
  bin.len = 40000L,
  # 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap between contacting bins
  gcb = "min2Mb", 
  chr.v = paste("chr", c(22:1, "X"), sep=""),
  nCPU = 2L,
  
  # Overlap parameters
  min.olap = 1L,
  max.gap = -1L,
  type.olap = "any",
  
  # LiftOver
  doLiftOver = F,
  LOchain = NULL
)

# rm(list=ls()); gc()