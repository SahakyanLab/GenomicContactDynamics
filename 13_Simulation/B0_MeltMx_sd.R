################################################################################
# Scale Cs values of cell/tissue by sd
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
Cs.norm.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
out.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_sd_primary_cohort")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################

################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:cp.v.len, chunks=nCPU), 
        .inorder=TRUE, .combine="rbind",
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
}
### END OF PARALLEL EXECUTION ###

################################################################################
library(compiler)
mapToHiCcontactPersistBins <- cmpfun(mapToHiCcontactPersistBins, options=list(suppressUndefined=TRUE))
################################################################################
#-------------------------------------------------------------------------------
#---------------------------------------
#-------------------

end.time <- Sys.time()
end.time-start.time 

#rm(list=ls()); gc()