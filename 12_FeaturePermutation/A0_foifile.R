################################################################################
# Generate foifile rest (all foi minus priority foi)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
foi.dir = paste0(wk.dir, "/foifile")
foigroup.v = c("hg19", "HMbp", "TF")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(foigroup in foigroup.v){
  foi.v <- readLines(paste0(foi.dir, "/foifile_", foigroup))
  foipr.v <- readLines(paste0(foi.dir, "/foifile_", foigroup, "_pr"))
  foirest.v <- setdiff(foi.v, foipr.v)
  if(length(foirest.v)==0){
    print(paste0(foigroup, " skipped!"), quote=FALSE)
    next
  }
  write(foirest.v, file=paste0(foi.dir, "/foifile_", foigroup, "_rest"))
  print(paste0(foigroup, " done!"), quote=FALSE)
}

# rm(list=ls())
