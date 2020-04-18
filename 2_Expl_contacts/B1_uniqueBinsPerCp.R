################################################################################
# Make a dataframe (per minimum distance between contacts) of unique bins with 
# their Cp and chromosome info. Note that a bin can have multiple Cps. 
# deva, R/3.5.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
version

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    # PERSIST.MX and BINKSTRINV.MX directory
    data.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    # PERSIST.MX and BINKSTRINV.MX directory
    data.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = data.dir
### OTHER SETTINGS #############################################################
gcb.v = c("min05Mb")
chr.v = paste0("chr", c(22:1, "X"), sep="")
nCPU = 6L # Max vmem~16G
kmer.len = 7
# Identifier of BINKMER7.MX
affix = ""
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
toExport <- c("data.dir", "out.dir", "chr.v")
chr.v.len <- length(chr.v)

for(gcb in gcb.v){
  
  print(paste0(gcb, "..."), quote=FALSE)
  
  toExport <- unique(c(toExport, "gcb"))
  
  #### PARALLEL EXECUTION #########
  
  UNIQBIN.DF <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU),
                        .combine="rbind", .inorder=FALSE,
                        .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      chr <- chr.v[i]
      
      print(paste(chr, "..."), quote=FALSE)
      
      # Load PERSIST.MX
      load(paste0(data.dir, "/", chr, "_Persist_", gcb, ".RData"))
      cp.v <- sort(as.numeric(unique(PERSIST.MX$ntis)), decreasing=FALSE)
     
      mx <- list()
      for(cp in cp.v){
        
        log.cp <- PERSIST.MX$ntis==cp
        # Get unique bins per Cp
        uniq.bin <- unique(
          c( unique(PERSIST.MX$hits[log.cp, "i"]), unique(PERSIST.MX$hits[log.cp, "j"]) )
        )
        mx[[as.character(cp)]] <- cbind( cp=rep(cp), bin=as.numeric(uniq.bin) )
        
        rm(log.cp, uniq.bin)
        
      } # cp.v for loop end
      
      df <- data.frame(chr=rep( strsplit(x=chr, split="chr")[[1]][2] ),
                       do.call("rbind", mx),
                       row.names=NULL, stringsAsFactors=FALSE) 
      
      # Load BINKMER.MX
      load(file=paste0(data.dir, "/", chr, "_BinKmer", kmer.len, "_", gcb, 
                       affix, ".RData"))
      # Bins populated by NAs in k-meric count matrix due to Ns in sequences
      binsWNs <- BINKMER7.MX[is.na(BINKMER7.MX[,"AAAAAAA"]), "bins"]
      rm(BINKMER7.MX); gc()
      
      print(paste(chr, "done!"), quote=FALSE)
      
      return( within(data=df,
                     {Ns <- NA
                     Ns[bin%in%as.numeric(binsWNs)] <- 1}
      ) )
     
    })
    
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  if( all(unique(UNIQBIN.DF$Ns)%in%c(1, NA)) ){
    save(UNIQBIN.DF, file=paste0(out.dir, "/chrALL_Uniqbin_", gcb, ".RData"))
    rm(UNIQBIN.DF); gc()
  }

  print(paste0(gcb, "done!"), quote=FALSE)

} # gcb.v for loop end

# rm(list=ls())