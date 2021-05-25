################################################################################
# Generate a lighter PERSIST.MX containing only i,j columns (no rownames) and
# Cp value. Keep original list format.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/3_AnnotationVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/3_AnnotationVsPersist"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc") 
out.dir = paste0(persist.dir, "/out_persistmx_ijonly") 
### OTHER SETTINGS #############################################################
gcb = "min05Mb"
chr.v = paste0("chr", c(1:22, "X"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste(gcb, "..."), quote=F)

for(chr in chr.v){
  
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  PERSIST.MX <- list(hits=PERSIST.MX$hits[,c("i","j")], ntis=PERSIST.MX$ntis)
  rownames(PERSIST.MX$hits) <- NULL
  save(PERSIST.MX, file=paste0(out.dir, "/", chr, "_Persist_", gcb, ".RData"))
  #load(file=paste0(out.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  print(paste(chr, " done!"), quote=F)
  
}

# rm(list=ls()); gc()