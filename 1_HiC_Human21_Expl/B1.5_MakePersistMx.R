################################################################################
# Make PERSIST.MX using HiCNorm matrices and make sure that contacts present
# per tissue and corresponding Cp is identical with PERSIST.MX from raw values.
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
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rawpmx.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts")
out.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
### OTHER SETTINGS #############################################################
norm = "HiCNorm_primary_cohort"
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 1L
gcb = "min05Mb"
min.dist = 2000000
hic.resol = 40000
top.perc = 100
min.tiss = 1
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/TrantoRextr/GEN_getMELTMXpersist.R"))  
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, " ", norm), quote=FALSE)

melt.dir <- paste0(melt.dir, "/", norm)
chr.v.len <- length(chr.v)

toExport <- c("rawpmx.dir", "melt.dir", "out.dir", "chr.v", "min.dist", 
              "hic.resol", "top.perc", "min.tiss", "gcb")
#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
        .inorder=TRUE, .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    chr <- chr.v[i]
    load(paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))
    PERSIST.MX <- getMELTMXpersist(MELT.MX=MELT.MX, min.dist=min.dist,
                                   hic.resol=hic.resol, top.perc=top.perc,
                                   min.tiss=min.tiss)
    rm(MELT.MX); gc()
    save(PERSIST.MX, file=paste0(out.dir, "/", chr, "_Persist_", gcb, ".RData"))

    newPMX <- PERSIST.MX
    rm(PERSIST.MX); gc()
    ct.v <- setdiff(colnames(newPMX$hits), c("i","j"))
    if( (length(ct.v)!=ncol(newPMX$hits)-2) ){ stop(paste0(chr, ": Checkpoint")) }
    #-------------------Compare with PERSIST.MX made from raw Cs
    load(paste0(rawpmx.dir, "/", chr, "_Persist_", gcb, ".RData"))
    
    if( !identical(PERSIST.MX$control, newPMX$control) ){
      print(paste0(chr, ": Control not consistent."), quote=FALSE)
    } 
  
    # Identical Cps?
    if( !identical(PERSIST.MX$ntis, newPMX$ntis) ){
      print(paste0(chr, ": Cp not consistent."), quote=FALSE)
      next
    } 
    
    # Identical ij?
    if( !identical(PERSIST.MX$hits[,c("i","j")], newPMX$hits[,c("i","j")]) ){
      print(paste0(chr, ": ij not consistent."), quote=FALSE)
      next
    } 
    
    PERSIST.MX$hits <- data.matrix(PERSIST.MX$hits[,ct.v])
    PERSIST.MX$hits[PERSIST.MX$hits!=0] <- 1
    newPMX$hits <- data.matrix(newPMX$hits[,ct.v])
    newPMX$hits[newPMX$hits!=0] <- 1
    # Identical set of contacts per tissue?
    if( !identical(PERSIST.MX$hits, newPMX$hits) ){
      print(paste0(chr, ": Set of contacts per tissue not consistent."), quote=FALSE)
      next
    }
    rm(PERSIST.MX, newPMX, ct.v); gc()
    print(paste0(chr, " done!"), quote=FALSE)
  }
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()