################################################################################
# Subset of PERSIST.MX containing only contacts of chosen Cp and cell/tissue
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_basePersist")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste0("chr", c(1:22, "X"), sep="")
# The code will also include topCP dynamic contacts; the code automatically
# includes Cp=1 contacts for the study of dynamic hubs. 
topCP = 4; cp.v = 1:21
ct.v = c("FC","ESC", "LC")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cp.v <- sort(as.numeric(cp.v), decreasing=FALSE)
cp.v <- c(1, tail(cp.v, n=topCP))

for(gcb in gcb.v){
  for(chr in chr.v){
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    incl.TF <- PERSIST.MX$ntis%in%cp.v
    # Remove unnecessay data from PERSIST.MX
    drop.nme <- names(PERSIST.MX)[!names(PERSIST.MX)%in%c("hits", "ntis")]
    for(nme in drop.nme){
      PERSIST.MX[[nme]] <- NULL
    }
    PERSIST.MX$hits <- PERSIST.MX$hits[incl.TF, c("i","j", ct.v)]
    PERSIST.MX$ntis <- PERSIST.MX$ntis[incl.TF]
    rm(incl.TF); gc()
    save(PERSIST.MX, file=paste0(out.dir, "/", chr, "_Persist_", gcb, "_topCP",
                                 topCP, ".RData"))
    print(paste0(gcb, "_", chr, " done!"), quote=FALSE)
  }
}
         
# rm(list=ls())