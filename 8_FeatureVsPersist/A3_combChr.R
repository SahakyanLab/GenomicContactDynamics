################################################################################
# Consolidate bin data from FETA.MX into Cp categories. 
# Only contacting bins present in given cell line / tissue were considered.
# deva, R/3.5.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=TRUE)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/8_FeatureVsPersist")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/ltamon" #"/t1-data/user/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/10_ChromatinFeatures")
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

# PERSIST.MX directory
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
# feta.dir is FETA.MX chrALL directory
feta.dir = paste0(wk.dir, "/out_FETA_raw_repeats")
out.dir = paste0(wk.dir, "/out_FETACP_raw_repeats")
bincount.dir = paste0(wk.dir, "/out_bincount")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw_repeats")
# If foifile = NULL, all files in foi.dir
foifile = NULL #paste0(wk.dir, "/foifile/foifile_sharedisofpcoding")
nCPU = 2L # stick to 5 NCPU, ~300G
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)

foi.v.len <- length(foi.v)
toExport <- c("foi.v", "persist.dir", "gcb", "feta.dir", "out.dir")

#### PARALLEL EXECUTION #########

foreach(itr=isplitVector(1:foi.v.len, chunks=nCPU), .inorder=FALSE, 
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    foi <- tail(x=strsplit(x=foi.v[i], split="\\/")[[1]], n=1)
    celltiss <- strsplit(x=foi, split="ct\\_|\\_foi")[[1]][2]
    foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
    
    print(paste0(foi, "..."), quote=FALSE)
    
    # Load FETA.MX (chrALL); not all features are present in all chr
    load(file=paste0(feta.dir, "/chrALL_", gcb, "_", foi, ".RData"))
    chr.v <- unique(FETA.MX$chr)
    chr.v.len <- length(chr.v)
    
    # FETA.MX values are the number of overlaps of feature in each bin.
    # Convert these values (>0) to 1 as we are interested only in %bin overlapping
    # with feature.
    
    # chr and bin columns
    temp <- FETA.MX[,1:2]
    FETA.MX[FETA.MX>0 & !is.na(FETA.MX)] <- 1L
    FETA.MX[,1:2] <- temp; rm(temp)
    
    for(c in 1:chr.v.len){
      
      chr <- chr.v[c]
      
      # Load PERSIST.MX
      load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      # Filter contacting bins by cell line/tissue
      if( celltiss%in%colnames(PERSIST.MX$hits[,-(1:2)]) ){
        TF <- PERSIST.MX$hits[,celltiss]>0
      } else {
        TF <- rep(TRUE, times=length(PERSIST.MX$ntis)) 
      }
              
      binsPerCp.ct <- by(data=PERSIST.MX$hits[TF,c("i","j")],
                         INDICES=PERSIST.MX$ntis[TF],
                         FUN=function(x){ unique( c(unique(x$i), unique(x$j)) ) }
      )
      rm(TF, PERSIST.MX); gc()
      #---------------------------------------
      if( !exists("FETACP.MX") ){
        # Initialize output matrix
        FETACP.MX <- matrix(data=0, nrow=21, ncol=ncol(FETA.MX)-2)
        dimnames(FETACP.MX) <- list( 1:21, colnames(FETA.MX)[-(1:2)] )
        
        bincount <- vector(mode="list", length=2)
        names(bincount) <- c("exp", "act")
        bincount[["exp"]] <- rep(0, times=21)
        names(bincount[["exp"]]) <- names(binsPerCp.ct)
        bincount[["act"]] <- FETACP.MX
      }
      
      chr.TF <- FETA.MX[,"chr"]==chr
      for( cp in names(binsPerCp.ct) ){
        ind <- which(FETA.MX[chr.TF,"bin"]%in%binsPerCp.ct[[cp]])
        FETACP.MX[as.numeric(cp),] <- FETACP.MX[as.numeric(cp),] + colSums(FETA.MX[chr.TF, -(1:2)][ind,], na.rm=TRUE) 
        #---------------------------------------
        bincount[["exp"]][as.character(cp)] <-  bincount[["exp"]][as.character(cp)] + length(binsPerCp.ct[[cp]])
        bincount[["act"]][as.numeric(cp),] <- bincount[["act"]][as.numeric(cp),] + apply( X=FETA.MX[chr.TF, -(1:2)][ind,], MARGIN=2, 
                                                                                          FUN=function(x){ sum(!is.na(x)) } )
        rm(ind)
      }
      rm(binsPerCp.ct); gc()
      #print(chr, quote=FALSE)
      
    } # chr.v for loop end
    
    FETACP.MX <- FETACP.MX/bincount[["act"]]*100
    save(FETACP.MX, file=paste0(out.dir, "/chrALL_", gcb, "_", foi, "_fetacp.RData"))
    
    if( !file.exists(paste0(bincount.dir, "/chrALL_", gcb, "_", celltiss, 
                            "_bincountplot.RData")) ){
      save(bincount, file=paste0(bincount.dir, "/chrALL_", gcb, "_", celltiss,
                                 "_bincountplot.RData"))
    }
    
    rm(FETACP.MX); gc()
    
    print(paste0(foi, " done!"), quote=FALSE)
    
  } # itr for loop end
    
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()










