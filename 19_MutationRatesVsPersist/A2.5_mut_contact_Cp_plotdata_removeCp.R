################################################################################
# Remove Cp column from IJ.MUT to decrease size
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
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/19_MutationRatesVsPersist")
src.dir = out.dir = paste0(wk.dir, "/out_mut_contact_Cp_plotdata")
### OTHER SETTINGS #############################################################
nCPU = 3
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
src.files <- list.files(path=src.dir)
src.files.len <- length(src.files)

print(paste0("Orig number of files: ", src.files.len), quote=F)

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:src.files.len, chunks=nCPU), .inorder=F
        
) %op% {
  
  
  for(i in itr){
    
    fle <- src.files[[i]]
    load(paste0(src.dir, "/", fle))
    file.remove(paste0(src.dir, "/", fle))
    
    IJ.MUT <- IJ.MUT[,"value", drop=F]
    save(IJ.MUT, file=paste0(out.dir, "/", fle))
    rm(IJ.MUT, fle)
    
  }
  
}
### END OF PARALLEL EXECUTION ###

out.files.len <- length(list.files(path=out.dir))
print(paste0("New number of files: ", out.files.len), quote=F)

if(src.files.len != out.files.len){
  stop("Number of files not the same.")
}

# rm(list=ls()); gc()