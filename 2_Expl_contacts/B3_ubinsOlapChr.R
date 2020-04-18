################################################################################
# Determine the percentage of chr 40-kb bins that participate in Hi-C long-range
# contact; do per chromosome per cell/tissue
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_ubinsOlapChr")
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste0("chr", c(1:22, "X"), sep="")
nCPU = 2L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
          "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
ct.v.len <- length(ct.v)
chr.v.len <- length(chr.v)
chrLen.df <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))
rownames(chrLen.df) <- chrLen.df$chromosome
chrLen.v <- chrLen.df[chr.v, "bins.40kb"]
names(chrLen.v) <- chr.v
rm(chrLen.df)

for(gcb in gcb.v){
  toExport <- c("ct.v", "ct.v.len", "chr.v", "persist.dir", "gcb", "chrLen.v")
  #### PARALLEL EXECUTION #########
  OLAPCHR.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                       .inorder=TRUE, .combine="rbind",
                       .export=toExport, 
                       .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      chr <- chr.v[i]
      chrLen <- chrLen.v[chr]
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      perc.v <- rep(NA, times=ct.v.len)
      names(perc.v) <- ct.v
      for(ct in ct.v){
        log <- PERSIST.MX$hits[,ct]!=0L
        ubins <- unique( c(unique(PERSIST.MX$hits$i[log]), 
                           unique(PERSIST.MX$hits$j[log])) )
        # Percentage of chr ubins forming LR contact
        perc.v[ct] <- length(intersect(1:chrLen, ubins))/chrLen*100
      }
      return(perc.v)
    })
    do.call("rbind", chunk)
  }
  ### END OF PARALLEL EXECUTION ###
  totbins <- sum(chrLen.v)
  chrALL <- apply(OLAPCHR.MX, MARGIN=2, FUN=function(ctperc){
    sum(ctperc/100*chrLen.v)/totbins*100
  })
  OLAPCHR.MX <- rbind(OLAPCHR.MX, chrALL=chrALL); rm(chrALL, totbins)
  OLAPCHR.MX <- cbind(nbins=c(chrLen.v, sum(chrLen.v)), OLAPCHR.MX)
  dimnames(OLAPCHR.MX) <- list(c(chr.v, "chrALL"), c("nbins", ct.v))
  save(OLAPCHR.MX, file=paste0(out.dir, "/chrALL_", gcb, "_ubinsOlapChr.RData"))
  rm(OLAPCHR.MX); gc()
  print(paste0(gcb, " done!"), quote=FALSE)
}

# rm(list=ls())



 

