################################################################################
# Count long-range contacts per cell line/tissue using PERSIST.MX
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/1_Count_contacts"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_count")
### OTHER SETTINGS #############################################################
gcb.v <- c("min2Mb", "min05Mb")
chr.v <- paste("chr", c(1:22, "X"), sep="") 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(gcb in gcb.v){
  
  outPerChr <- foreach( chr=chr.v, 
                        .inorder=TRUE, .combine="rbind",
                        .export=c("gcb", "persist.dir"), 
                        .noexport=ls()[!ls()%in%c("gcb", "persist.dir")]
  ) %do% {
    
    
    # Load PERSIST.MX
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    
    # Count contacts
    countPerCp <- rep(0, times=21L)
    cp.v <- 1:21
    for(cp in cp.v){
      
      countPerCp[cp] <- sum(PERSIST.MX$ntis==cp)
      
    }
    rm(PERSIST.MX); gc()
    
    print(paste0(chr, " done!"), quote=FALSE)
    
    return( data.frame(matrix(countPerCp, ncol=21)) )
    
  }
  
  outPerChr <- rbind(outPerChr, colSums(outPerChr))
  rownames(outPerChr) <- c(chr.v, "chrALL")
  colnames(outPerChr) <- cp.v
  
  write.csv(x=outPerChr, file=paste0(out.dir, "/", gcb, "_HiCLRcontactsCount.csv"),
            row.names=TRUE, quote=FALSE) 
  
} # gcb.v for end

# rm(list=ls())