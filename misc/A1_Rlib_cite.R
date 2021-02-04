################################################################################
# Get citation of R packages used
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/misc"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir <- paste0(wk.dir, "/out_Rlib_cite")
### OTHER SETTINGS #############################################################
pckgsPath = paste0(wk.dir, "/Rpckgs.csv")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
pckgs.v <- as.character(read.csv(file=pckgsPath, header=TRUE)[,1])
  
#write(x=NULL, file=paste0(out.dir, "/gcd_supp_Rpckgs.bib"))
for(p in pckgs.v){
  x <- citation(p)
  write(x=toBibtex(x), file=paste0(out.dir, "/gcd_supp_Rpckgs.bib"), append=TRUE)
}

# rm(list=ls()); gc()