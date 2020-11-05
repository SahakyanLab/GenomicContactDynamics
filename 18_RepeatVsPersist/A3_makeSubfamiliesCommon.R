################################################################################
# Determine proper binning of heatmap
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "fam" # "fam" | "subfam"
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
out.dir = paste0(wk.dir, "/out_makeSubfamiliesCommon/", rep.group)
### OTHER SETTINGS #############################################################
gcb.v = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 1
out.name = "GiorPubl372NotIncl"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(rep.group=="fam"){
  agerank <- read.csv(file=paste0(agerank.dir, "/repfam.csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,"repFamily"]
} else if(rep.group=="subfam"){
  agerank <- read.csv(file=paste0(agerank.dir, "/plot_GiorPubl372rankrepFamilies.csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,"repName"]
} else {
  stop("Invalid rep.group.")
}

chr.v.len <- length(chr.v)

for(gcb in gcb.v){
  
  print(paste0(gcb, "..."), quote=FALSE)
  
  toExport <- c("chr.v", "minelm.dir", "agerank", "gcb")
  #### FOREACH EXECUTION #########
  
  out <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                 .combine="rbind", .inorder=TRUE,
                 .export=toExport, 
                 .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      chr <- chr.v[i]
      load(paste0(minelm.dir,"/", chr, "_MinElm_", gcb, ".RData"))
      agerankNotIncl <- setdiff(agerank, colnames(MINELM.MX))
      rm(MINELM.MX); gc()
      agerankNotIncl.len <- length(agerankNotIncl)
      if(is.null(agerankNotIncl)){
        agerankNotIncl <- "none"
      }
      c( chr, agerankNotIncl.len, paste(agerankNotIncl, collapse=";") )
    })
    
    do.call("rbind", chunk)
    
  }
  
  colnames(out) <- c("chr", "number", "agerankNotincl")
  write.csv(x=out, file=paste0(out.dir, "/", gcb, "_", out.name, ".csv"),
            row.names=FALSE)
  
  ### END OF FOREACH EXECUTION ###
  
}

# rm(list=ls())
