################################################################################
# For each contact in a chromosome, get mean, median, sd of Hi-C values across 
# tissues as well as the number of tissues having that contact. Do this using 
# raw and the types of normalised matrices. Contacts consired are only those
# present in at least one tissue (those in MELT.MX$upper.tri). These set of 
# contacts the same for RAW_primary_cohort and HiCNorm_primary_cohort but not for
# HiCNorm_QQ_primary_cohort. The code also compares LR contacts before and
# after normalisation (within and across datasets). It checks whether the set of
# normalised LR contacts is the same with the raw, if the Cp of contacts is 
# consistent and if each contact is expressed in the same set of tissues.
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/1_HiC_Human21_Expl"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/1_HiC_Human21_Expl"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts")
out.dir = paste0(wk.dir, "/out_rawVsnorm")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 3L
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
norm.v <- c("RAW_primary_cohort", "HiCNorm_primary_cohort", "HiCNorm_QQ_primary_cohort")
print(paste0(norm.v[1], " as reference for checkpoints."), quote=FALSE)

chr.v.len <- length(chr.v)
norm.v.len <- length(norm.v)

toExport <- c("melt.dir", "out.dir", "chr.v", "norm.v", "norm.v.len")
#### PARALLEL EXECUTION #########
IJ.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                 .inorder=FALSE, .export=toExport, 
                 .noexport=ls()[!ls()%in%toExport]
) %op% {
  for(i in itr){
    chr <- chr.v[i]
    for(n in 1:norm.v.len){
      
      nrm <- norm.v[n]
      id <- paste0(chr, " ", nrm)
      
      load(paste0(melt.dir,"/", nrm, "/human_", chr, "_allcontacts.RData"))
      MELT.MX$upper.tri.nocontact <- MELT.MX$rest <- NULL
      
      row.means   <- rowMeans(MELT.MX$upper.tri[,-c(1,2)])
      row.sds     <- apply(X=MELT.MX$upper.tri[,-c(1,2)], MARGIN=1, FUN=sd)
      row.medians <- apply(X=MELT.MX$upper.tri[,-c(1,2)], MARGIN=1, FUN=median)
      row.num0s   <- apply(X=MELT.MX$upper.tri[,-c(1,2)], MARGIN=1, FUN=function(rw){sum(rw==0)})
      ind.0s <- lapply(X=MELT.MX$upper.tri[,-(1:2)], FUN=function(col){
        which(col==0)
      })
      
      if(n==1){
        ij.ref <- MELT.MX$upper.tri[,c("i", "j")]
        ct.v <- dimnames(MELT.MX$upper.tri[,-(1:2)])[[2]]
        row.num0s.ref <- row.num0s
        ind.0s.ref <- ind.0s; rm(ind.0s)
      } else {
        if( !identical(ij.ref, MELT.MX$upper.tri[,c("i", "j")]) ){
          print(paste0(id, ": Different set of contacts!"), quote=FALSE)
        } 
        if( !identical(row.num0s.ref, row.num0s) ){
          print(paste0(id, ": Different Cp of contacts!"), quote=FALSE)
        } 
        if( !identical(ind.0s.ref[ct.v], ind.0s[ct.v]) ){
          print(paste0(id, ": Different cell/tissue pattern of contact formation!"), 
                quote=FALSE)
        } 
        rm(ind.0s)
      }
      rm(MELT.MX); gc()
      IJSTAT.MX <- cbind.data.frame(chr=chr, row.means=row.means, 
                                    row.medians=row.medians,
                                    row.sds=row.sds, row.num0s=row.num0s, 
                                    stringsAsFactors=FALSE)
      rm(row.means, row.sds, row.medians, row.num0s)
      save(IJSTAT.MX, file=paste0(out.dir, "/human_", chr, "_uppertri_", nrm, 
                                  "_ijstat.RData"))
      rm(IJSTAT.MX, nrm); gc()
      print(paste0(id, " done!"), quote=FALSE)
      
    } # norm.v.len for loop end
    rm(chr, ij.ref, row.num0s.ref, ind.0s.ref, ct.v); gc()
  } # itr for loop end
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()