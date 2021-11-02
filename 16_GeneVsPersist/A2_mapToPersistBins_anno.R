################################################################################
# Map features to HiC contact bins 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
    data.dir = "/Users/ltamon/Database"
    feat.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
feat.dir = paste0(data.dir, "/ucsc_tables/hsa_geneAnno")
out.dir = paste0(wk.dir, "/out_mapToPersistBins_anno")
### OTHER SETTINGS #############################################################
#featurefile.v = paste(feat.dir, c("/hg19anno_ALL", "/hg19anno_NM", "/hg19anno_NR"), sep="")
featurefile.v = paste(feat.dir, c("/hg19annoLTr_ALL", "/hg19annoLTr_NM", "/hg19annoLTr_NR"), sep="")
# Output identifier; should correspond with featurefile.v
#out.name.v = c("ALL_hg19_40KbBinAnno", "NM_hg19_40KbBinAnno", "NR_hg19_40KbBinAnno")
out.name.v = c("LTr_ALL_hg19_40KbBinAnno", "LTr_NM_hg19_40KbBinAnno", "LTr_NR_hg19_40KbBinAnno")
# HiC res
bin.len = 40000L
# Path to file
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb = "min2Mb"
chr.v = paste("chr", c(22:1, "X"), sep="")
nCPU = 2L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) 
library(compiler)
library(GenomicRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/loadRData.R"))
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/mapToHiCcontactPersistBins.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
len <- length(featurefile.v)
for(i in 1:len){
  
  mapToHiCcontactPersistBins(
    
    # PERSIST.MX object directory
    persist.dir=persist.dir, 
    
    ## Feature directory
    feat.dir=feat.dir,
    # Separate files per chromosome?
    featSepByChr = FALSE,
    #---- if featSepByChr = TRUE
    # Program will look for each file on feat.dir using chr number eg. *chrX* as pattern
    # Beware of duplicate files for each chr
    
    #---- if featSepByChr = FALSE
    # Path to file
    featurefile=featurefile.v[i],
    # Combine output of chromosomes?
    combineOut = FALSE,
    
    feat.header = TRUE,
    # Column numbers or names
    start.coord = "txStart",
    end.coord = "txEnd",
    # Table coordinate system
    # zero-based (0-start, half-open) or one-based (1-start, fully-closed (1-based))
    # Coordinates stored in UCSC Genome Browser tables are zero-based while those on 
    # web interfaces are one-based
    feat.coordSys = "zero-based",
    chr.col = "chrom",
    # Remove irrelevant columns
    # Recommended that the featurefile should have a column uniquely identifying each entry/row
    drop.col = c("bin", "name", "strand", "cdsStart", "cdsEnd", "exonCount", "exonStarts", 
                "exonEnds", "score", "cdsStartStat", "cdsEndStat", "exonFrames"), 
    # Column names after dropping chr.col and drop.col (follow order in original table)
    # If NULL, keep original column names
    coluNames = NULL,
    
    # Output directory
    out.dir=out.dir,
    # Output identifier
    out.name=out.name.v[i],
    
    # HiC resolution
    bin.len=bin.len,
    # 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap between contacting bins
    gcb=gcb, 
    chr.v=chr.v,
    nCPU=nCPU,
    
    # Overlap parameters
    min.olap = 1L,
    max.gap = -1L,
    type.olap = "any"
  )
  
} # featurefile.v for loop end
  
# rm(list=ls()); gc()

  



