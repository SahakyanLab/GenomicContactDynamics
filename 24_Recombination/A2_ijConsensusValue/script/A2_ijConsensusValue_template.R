################################################################################
# Process recombination rate data contact wise. Filter bins (rate set to NA) based 
# on count of data points it contains with the min.countPerBin argument using 
# FEATURE.BIN.MX countPerBin column. 
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
    home.dir = "/project/sahakyanlab/ltamon" 
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/24_Recombination")
featbin.dir = paste0(wk.dir, "/out_mapToHiCcontactPersistBins")
out.dir = paste0(wk.dir, "/out_ijConsensusValue")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrarr1.repl" 
Cp = arr2.repl
nCPU = 1
chunk.size = 100 # Number of contacts to be processed at the same time per CPU
consensus.method = "MEDIAN.MEDIAN" #"MEAN.MEAN" #"MEDIAN.MEDIAN" 
min.countPerBin = 3 # 3 to be consistent with RT data
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/ijconsensusFUN.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(featbin.dir, "/", chr, "_", gcb, "_recomRates_2011_01_phaseII_B37_Myers.RData"))
if( any(duplicated(FEATURE.BIN.MX[,"bin"])) ){
  stop(paste0(chr, ": Duplicated bin in FEATURE.BIN.MX."))
}

# Set rate to NA of bins with < min.countPerBin count of rates
FEATURE.BIN.MX[ FEATURE.BIN.MX[,"countPerBin"] < min.countPerBin ,"V5"] <- as.character(NA)
rates <- setNames(FEATURE.BIN.MX[,"V5"], nm=FEATURE.BIN.MX[,"bin"])
rm(FEATURE.BIN.MX)

load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
ij.mx <- ij.mx[ ij.mx[,"Cp"] == Cp,-3, drop=F ]
rm(PERSIST.MX)

ij.len <- length(ij.mx[,1])
toExport <- c("rates", "ij.mx")
#### PARALLEL EXECUTION #########
IJ.RATES <- foreach(itr=isplitVector(1:ij.len, chunkSize=chunk.size), .inorder=T, .combine="c",
                    .export=toExport, .noexport=ls()[!ls()%in%toExport]
                    
) %op% {
  
  I.RATES <- strsplit(x=rates[ as.character(ij.mx[itr,"i"]) ], split=";")
  J.RATES <- strsplit(x=rates[ as.character(ij.mx[itr,"j"]) ], split=";")
  
  I.RATES <- lapply(I.RATES, FUN=as.numeric)
  J.RATES <- lapply(J.RATES, FUN=as.numeric)

  # Get consensus value per contact
  
  ij.rates <- sapply(X=1:length(I.RATES), simplify=T, FUN=function(ind){
    
    if( all( is.na(I.RATES[[ind]]) ) | all( is.na(J.RATES[[ind]]) ) ){
      return(NA)
    } else {
      return( ijconsensusFUN(i=I.RATES[[ind]], j=J.RATES[[ind]], METHOD=consensus.method) )
    }
  
  })
  
  return(ij.rates)
  
}
### END OF PARALLEL EXECUTION ###

if( length(IJ.RATES) != ij.len ){
  stop(paste0(chr, ": Length of rates and contacts different."))
}

OUT.dir <- paste0(out.dir, "/tmp")
if( !dir.exists(OUT.dir) ){ dir.create(OUT.dir, showWarnings=T) }
save(IJ.RATES, file=paste0(OUT.dir, "/mthd_", consensus.method, "_Cp", Cp, "_", chr, "_", gcb, 
                           "_min.countPerBin", min.countPerBin, "_ij_recomRates_2011_01_phaseII_B37_Myers.RData") )

# rm(list=ls()); gc()