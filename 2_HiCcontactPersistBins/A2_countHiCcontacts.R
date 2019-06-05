################################################################################
# Count the number of HiC contacts per contact persistence (Cp)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start_time <- Sys.time()

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    objective.dir = "/Users/ltamon/DPhil/GenomicContactDynamics"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    feat.dir = "/Users/ltamon/Database/topo_ChIPseq/final"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    objective.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics"
    persist.dir = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
    feat.dir = "/t1-data/user/ltamon/Database/topo_ChIPseq/final"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
### OTHER SETTINGS #############################################################
gcb.v <- c("min2Mb", "min05Mb")
chr.v <- paste("chr", c(1:22, "X"), sep="") #c("chr21", "chr20", "chr19") 
nCPU <- 23L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)

out <- lapply(X=gcb.v, FUN=function(gcb){
  
  #### PARALLEL EXECUTION #########
  
  outPerChr <- foreach( itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                             .inorder=TRUE, .combine="rbind",
                             .export=c("chr.v", "gcb", "persist.dir"), 
                             .noexport=ls()[!ls()%in%c("chr.v", "gcb", "persist.dir")]
  ) %op% {
    
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      chr <- chr.v[i]
      # Load PERSIST.MX
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      
      # Count contacts
      countPerCp <- rep(0, times=21L)
      cp.v <- 1:21
      for(cp in cp.v){
        
        countPerCp[cp] <- sum(PERSIST.MX$ntis==cp)
        
      }
      rm(PERSIST.MX); gc()
      return( data.frame(matrix(countPerCp, ncol=21)) )
    }) # sapply end
    
    return(do.call("rbind", chunk))
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  outPerChr <- cbind( c(chr.v, "chrALL"),
                      rbind(outPerChr, colSums(outPerChr)) )
  return( cbind(rep(gcb), outPerChr) )
    
}) # gcb.v lapply end

out <- do.call("rbind", out)
colnames(out) <- c("gcb", "chr", as.character(1:21))
write.table(x=out, file=paste0(persist.dir, "/HiCcontactsCount"),
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE) 

end_time <- Sys.time()
end_time-start_time 

# rm(list=ls())