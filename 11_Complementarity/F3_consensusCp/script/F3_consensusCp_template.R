################################################################################
# Calculate consensus Cp based on overlapping Cp values with each hg38 to
# hg19 converted contact save as semi-colon separated string from F2_assignCp.R. 
# Code expects to use consensus functions giving numeric values. Output maintains
# order of the source of Cp strings, which is the same order as CII.MX hg38 
# contacts.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) 
options(warn=1)

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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
Cp.dir = paste0(wk.dir, "/out_assignCp") # 1-based
out.dir = paste0(wk.dir, "/out_consensusCp")
### OTHER SETTINGS #############################################################
src.id = "hg38ToHg19_LOwidth.min.bp30000_kmer_min2Mb_cons.methd2_finalStringOfCpsCp"
out.id = "hg38ToHg19_LOwidth.min.bp30000_kmer_min2Mb_consensusCp"
chrs = "chrCHRREPLACE"

# Tabulate works only for positive-integer-valued vector, so x should be the same
MODE <- function(x){ 
  tab.x <- tabulate(x) 
  mean(which(tab.x == max(tab.x)))
}
cons.funx.obj <- list("mean"=mean, "median"=median, "max"=max, "min"=min,
                      "MODE"=MODE)
nCPU = 2
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cons.funx.len <- length(cons.funx.obj)

for(chr in chrs){
  
  load(paste0(Cp.dir, "/", chr, "_", src.id, ".RData"))
  
  is.nonNA.Cp <- !is.na(consensusCp$value)
  nonNA.Cp.lst <- strsplit(x=consensusCp$value[is.nonNA.Cp], split=";")
  rm(consensusCp)
  nonNA.Cp.lst <- lapply(nonNA.Cp.lst, FUN=as.integer)
  
  print("Cp values:", quote=F)
  print( sort(unique(lengths(nonNA.Cp.lst))), quote=F ) 
  
  #### PARALLEL EXECUTION #########
  toExport <- c("cons.funx.obj", "nonNA.Cp.lst")
  
  nonNA.consCp <- foreach(itr=isplitVector(1:cons.funx.len, chunks=nCPU), .inorder=T, 
                          .combine="c", .export=toExport, .noexport=ls()[!ls()%in%toExport]
                         
  ) %op% {
    
    chunk <- sapply(itr, simplify=F, FUN=function(i){
      
      cons.funx <- unname( cons.funx.obj[[i]] )
      return( unlist(lapply(nonNA.Cp.lst, FUN=cons.funx)) )
      
    })
    
    return(chunk)
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  nonNA.consCp <- do.call("cbind", nonNA.consCp)
  consensusCp <- matrix(data=NA, nrow=length(is.nonNA.Cp), ncol=cons.funx.len,
                        dimnames=list( NULL, paste0(names(cons.funx.obj), ".consCp") )
                        )
  consensusCp[is.nonNA.Cp,] <- nonNA.consCp
  save(consensusCp, file=paste0(out.dir, "/", chr, "_", out.id, ".RData"))

  rm(is.nonNA.Cp, nonNA.Cp.lst, nonNA.consCp, consensusCp)
  
  print(paste0(chr, ": done!"), quote=F)
  
}

# rm(list=ls()); gc()
