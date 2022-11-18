################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")

val1.dir = paste0(wk.dir, "/out_constraints_GfreeSingleNorm/merged_final")  
val2.dir = paste0(wk.dir, "/out_constraints_hg19_rm_GfreeSingleNorm/merged_final") 
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr18"
compl.type = "kmer"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(val1.dir, "/", chr, "_", compl.type, "_", gcb, ".RData"))
a <- CII.MX
load(paste0(val2.dir, "/", chr, "_", compl.type, "_", gcb, ".RData"))
b <- CII.MX

# Objects
xxx.v  # vector
xxx.mx # matrix
xxx.df # dataframe
xxx.id # id 

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:cp.v.len, chunks=nCPU), 
        .inorder=T, .combine="rbind",
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

# rm(list=ls()); gc()