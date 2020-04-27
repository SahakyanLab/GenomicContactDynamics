################################################################################
# Make a matrix of the contacts per chr with their corresponding radial positions.
# Output IJ.FEAT.MX is cell/tissue-specific. 
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
feat.bin.dir = paste0(wk.dir, "/out_mapToHiCcontactBins_Chrom3D/", model.id)
out.dir = paste0(wk.dir, "/out_ContactRadDist/", model.id)
### OTHER SETTINGS #############################################################
celltiss = "FC" # "FC" | "ESC"
ploidy = "haploid"
gcb = "min2Mb"
chr.v = paste("chr", c(22:1, "X"), sep="")
# Only consider contacts with bins mapping to one domain; <= x or within [x,y]
setCountPerBin = c(1,1)
recordDropout <- TRUE
featureToAdd = c("iRadDist", "jRadDist")
nCPU = 5L #~40G
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel) 
library(itertools) 
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/makeContactVsFeatureTable.R"))
source(paste0(lib, "/HiC21celltissueCode.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(model.id, " ", celltiss, "..."), quote=FALSE)
out.name <- paste0(model.id, "_", ploidy) 
chr.v.len <- length(chr.v)
for(c in 1:chr.v.len){
  
  chr <- chr.v[c]
  
  # Load FEATURE.BIN.MX
  # Better to load this repeatedly and subset based on chr
  load( file=paste0(feat.bin.dir, "/chrALL_", gcb, "_", out.name, 
                    "_Chrom3DVsPersist.RData") )
  
  # Filter FEATURE.BIN.MX based on chromosome
  feature.bin.mx <- FEATURE.BIN.MX[FEATURE.BIN.MX$chr==chr, 
                                   c("bin", "radDist", "countPerBin")]
  rm(FEATURE.BIN.MX); gc()
  
  # Load PERSIST.MX
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  # Remove unnecessary data
  PERSIST.MX$valsum <- NULL
  PERSIST.MX$control <- NULL
  
  tot.ij <- length(PERSIST.MX$ntis)
  
  # Filter PERSIST.MX based on cell/tissue
  incl.TF <- PERSIST.MX$hits[[celltiss]]>0L
  # Number of contacts for the given cell/tissue and bins present even after liftOver
  ct.ij <- sum(incl.TF)
  PERSIST.MX$hits <- PERSIST.MX$hits[incl.TF,c("i", "j")]
  PERSIST.MX$ntis <- PERSIST.MX$ntis[incl.TF]
  rm(incl.TF); gc()

  # Initialize table where to add the features
  IJ.FEAT.MX <- makeContactVsFeatureTable(setCountPerBin=setCountPerBin,
                                          featureBinObj=feature.bin.mx,
                                          persistObj=PERSIST.MX,
                                          featureToAdd=featureToAdd)
  rm(PERSIST.MX); gc()
  
  incl.ij.ind <- which(IJ.FEAT.MX[,"include"]==1)
  # Remove include column
  IJ.FEAT.MX <- IJ.FEAT.MX[,-1]
  incl.ij.ind.len <- length(incl.ij.ind)
  IJ.FEAT.MX.len <- nrow(IJ.FEAT.MX)
  drop.ij.perc=paste0( round((1-(incl.ij.ind.len/IJ.FEAT.MX.len))*100, 
                             digits=4), "%" )
    if(recordDropout==TRUE){
      towrite <- paste(chr, gcb, tot.ij, ct.ij, incl.ij.ind.len, drop.ij.perc, 
                       paste(setCountPerBin, collapse="-"), sep=",")
      if(c==1){
        towrite <- c("chr,gcb,tot.ij,ct.ij, incl.ij, drop.ij.perc(reftoincl.ij), setCountPerBin[x,y]", towrite)
      }
    write(towrite, append=TRUE,
          file=paste0(out.dir, "/", gcb, "_", out.name, 
                      "_dropOutPercent_ContactRadDist"))
  }
  rownames(feature.bin.mx) <- feature.bin.mx$bin
  toExport <- c("feature.bin.mx", "IJ.FEAT.MX", "incl.ij.ind")
  
  #### PARALLEL EXECUTION #########
  
  add <- foreach( itr=isplitVector(1:incl.ij.ind.len, chunks=nCPU), 
                  .combine="rbind", .inorder=TRUE,
                  .export=toExport, 
                  .noexport=ls()[!ls()%in%toExport] 
  ) %op% {
    ir <- feature.bin.mx[as.character(IJ.FEAT.MX[incl.ij.ind[itr],"i"]),
                         "radDist"]
    ir <- as.numeric(unlist(lapply(X=strsplit(x=ir, split=";"), FUN=unique)))
    jr <- feature.bin.mx[as.character(IJ.FEAT.MX[incl.ij.ind[itr],"j"]),
                         "radDist"]
    jr <- as.numeric(unlist(lapply(X=strsplit(x=jr, split=";"), FUN=unique)))
    return(cbind(ir, jr))
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  # Add features to the initialized table
  IJ.FEAT.MX[incl.ij.ind, featureToAdd] <- add
  save(IJ.FEAT.MX, file=paste0(out.dir, "/", chr, "_", gcb, "_", 
                               out.name, "_ContactRadDist.RData"))
  
  rm(feature.bin.mx); gc()
  
  print(paste0(gcb, chr, ":", "ijRadDist done!"), quote=FALSE)
  
} # chr.v for loop end

# rm(list=ls())




