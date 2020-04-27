################################################################################
# Merge IJ.FEAT.MX of all chromosomes.
# deva, R/3.6.0-newgcc, gcc/4.9.2
# athena, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    persist.dir = "/home/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
data.dir = out.dir = paste0(wk.dir, "/out_ContactRadDist/", model.id)
### OTHER SETTINGS #############################################################
ploidy = "haploid"
gcb.v = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 1L
################################################################################
# LIBRARIES & DEPEND9898870ANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel) 
library(itertools) 
source(paste0(lib, "/loadRData.R"))
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(model.id, "_", ploidy) 
for(gcb in gcb.v){
  
  chr.v.len <- length(chr.v)
  toExport <- c("out.dir", "gcb", "chr.v", "chr.v.len")
  
  #### PARALLEL EXECUTION #########
  
  IJ.FEAT.MX <- foreach( itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                         .combine="rbind", .inorder=FALSE,
                         .export=toExport, 
                         .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    IJ.FEAT.MX.chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      chr <- chr.v[i]
      mx <- loadRData(file=paste0(data.dir, "/", chr, "_", gcb, "_", 
                                  out.name, "_ContactRadDist.RData"))
      # Get only included contacts
      test1 <- !is.na(mx[,"iRadDist"])
      test2 <- !is.na(mx[,"jRadDist"])
      mx <- mx[test1&test2,]
    })
    do.call("rbind", IJ.FEAT.MX.chunk)
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  save(IJ.FEAT.MX, file=paste0(out.dir, "/chrALL_", gcb, "_", 
                               out.name, "_ContactRadDist.RData"))
  
}

# rm(list=ls())




