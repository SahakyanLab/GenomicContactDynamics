################################################################################
# Generate .RData and .csv per chromosome containing complementarity values
# from all 3 methods. Format intended for sharing data.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar = TRUE)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
    CII.dir = paste0(wk.dir, "/z_ignore_git/out_constraints/merged_final")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints"
    CII.dir = paste0(wk.dir, "/out_constraints/merged_final")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
out.dir = paste0(wk.dir, "/out_complementarity_dataset")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 1
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
chr.v.len <- length(chr.v)

toExport <- c("CII.dir", "out.dir", "gcb", "chr.v")
#### PARALLEL EXECUTION #########

foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), .inorder=F,
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    chr <- chr.v[i]
    
    load(paste0(CII.dir, "/", chr, "_align_", gcb, ".RData"))
    align <- CII.MX
    
    load(paste0(CII.dir, "/", chr, "_kmer_", gcb, ".RData"))
  
    if( !identical(align[,c("i", "j", "Cp")], CII.MX[,c("i", "j", "Cp")]) ){
      stop("Contacts in CII.MX kmer and align in different order.")
    }
    
    CII.MX <- cbind(CII.MX, align=align[,"C||"])
    rm(align); gc()
    
    # Fix column name
    clnme <- colnames(CII.MX)
    if( !"C||"%in%clnme ){
      stop("C|| column not present.")
    }
    clnme[clnme=="C||"] <- "kmer"
    colnames(CII.MX) <- clnme
    
    save(x=CII.MX, file=paste0(out.dir, "/", chr, "_CII_", gcb, ".RData"))
    write.csv(x=CII.MX, file=paste0(out.dir, "/", chr, "_CII_", gcb, ".csv"),
              row.names=T)
    
  }
  
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()