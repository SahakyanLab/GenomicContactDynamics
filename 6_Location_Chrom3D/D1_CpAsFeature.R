################################################################################
# Use IJ.FEAT.MX to make a Cp featurefile compatible for mapping to DOMXYZR.DF.
# Featurefile should have rows corresponding to features and columns for chr, start
# and end coordinates of bin and feature name (bin). To reflect the number of contacts
# per Cp, unique bins per Cp are not taken. Output is cell/tissue specific. 
# because IJ.FEAT.MX is cell/tissue specific. Also, only non-promiscous bins included.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    os = "Mac"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
# IJ.FEAT.MX directory
ijfeat.dir = paste0(wk.dir, "/out_ContactRadDist/", model.id)
out.dir = paste0(wk.dir, "/out_CpAsFeature")
### OTHER SETTINGS #############################################################
ploidy.v = "haploid"
HiC.res = 4e4L
gcb  = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 3L # ~5G
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(itertools) 
library(doParallel) 
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)

for(ploidy in ploidy.v){
  
  out.name <- paste0(model.id, "_", ploidy)
  
  toExport=c("ijfeat.dir", "out.dir", "HiC.res", "gcb", "chr.v", "out.name") 
  
  #### PARALLEL EXECUTION #########
  
  FEATURE.DF <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                        .inorder=FALSE, .combine="rbind",
                        .export=toExport, 
                        .noexport=ls()[!ls()%in%toExport]
                        
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      chr <- chr.v[i]
      
      load(file=paste0(ijfeat.dir, "/", chr, "_", gcb, 
                       "_", out.name, "_ContactRadDist.RData"))
      
      # Get only included contacts (those with bins mapping to only one domain);
      # First 3 columns are i, j, Cp
      IJ.FEAT.MX <- IJ.FEAT.MX[ !is.na(IJ.FEAT.MX[,"iRadDist"]) & !is.na(IJ.FEAT.MX[,"jRadDist"]),]
      if( any(is.na(IJ.FEAT.MX[,"iRadDist"])) ){ stop(paste0(chr, ":NA iRadDist")) }
      if( any(is.na(IJ.FEAT.MX[,"jRadDist"])) ){ stop(paste0(chr, ":NA jRadDist")) }
      IJ.FEAT.MX <- IJ.FEAT.MX[,c("i", "j", "ntis")]
      ntis.uniq <- sort( as.numeric(unique(IJ.FEAT.MX[,"ntis"])) )
      
      # Classify bins based on Cp and convert them into coordinates in the output
      lst <- sapply(X=ntis.uniq, simplify=FALSE, FUN=function(ntis){
        bins <- as.integer(
          unlist(IJ.FEAT.MX[IJ.FEAT.MX[,"ntis"]==ntis,c("i", "j")], use.names=FALSE)
                           )
        # Length of bins should always be even
        if( (length(bins)%%2)!=0 ){ stop(paste0(chr, ":Number of bins not even.")) }
        bin.end   <- bins*HiC.res 
        bin.start <- bin.end-HiC.res+1 
        cbind.data.frame(chr=rep(chr), start=bin.start, end=bin.end,
                         Cp=rep(ntis),
                         stringsAsFactors=FALSE)
      })
      
      return(do.call("rbind", lst))
    
    }) # itr sapply end
    
    return(do.call("rbind", chunk))
  }
    
  ### END OF PARALLEL EXECUTION ###
  
  save(FEATURE.DF, file=paste0(out.dir, "/", out.name,"_", gcb, "_",
                               "CpFeat.RData"))
} # gcb.v for loop end

# rm(list=ls())