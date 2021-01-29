################################################################################
# Merge FETA.MX from all chr.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_mutPerBinSiteperChr")
feta.dir = paste0(out.dir, "/temp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", 1:22, sep="")
foi.v = c("All", "AToC", "AToG", "AToT", "CToA", "CToT", "CToG")
nCPU = 1L # stick to 5 NCPU, ~300G
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(3.doParallel)
library(itertools)
library(reshape)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
foi.v.len <- length(foi.v)
toExport <- c("foi.v", "chr.v.len", "chr.v", "gcb", "feta.dir", "out.dir")

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:foi.v.len, chunks=nCPU), .inorder=FALSE, 
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    foi <- tail(x=strsplit(x=foi.v[i], split="\\/")[[1]], n=1)
    celltiss <- strsplit(x=foi, split="ct\\_|\\_foi")[[1]][2]
    foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
    
    print(paste0(foi, "..."), quote=FALSE)
    
    if(plotOnly==FALSE){
      chr.v.len <- length(chr.v)
      FETA.MX <- foreach(c=1:chr.v.len, .combine="rbind") %do% {
        chr <- chr.v[c]
        # Load FETA.MX; not all features are present in all chr
        fle <- paste0(feta.dir, "/", chr, "_", gcb, "_", foi, ".RData")
        if(file.exists(fle)){
          load(file=fle)
          return(cbind.data.frame(chr=chr, FETA.MX, stringsAsFactors=FALSE))
        } else {
          return(NULL)
        }
      } 
      save(FETA.MX, file=paste0(out.dir, "/chrALL_", gcb, "_", foi, ".RData"))
    } else {
      load(file=paste0(out.dir, "/chrALL_", gcb, "_", foi, ".RData"))
    }
    #---------------------------------------
    ## Boxplot of overlap values (absolute)
    #BP.DF <- melt(FETA.MX[,!colnames(FETA.MX)%in%c("chr", "bin")])
    #rm(FETA.MX); gc()
    #colnames(BP.DF) <- c("bin", "numOlaps")
    #max.val <- max(BP.DF$numOlaps, na.rm=TRUE)
    
    #bitmap(file=paste0(out.dir, "/max", max.val, "_chrALL_", gcb, "_", foi, ".png"), 
    #       height=10, width=10, units="in", res=500, type="png256")
    #par(oma=c(0,0,0,0), mar=c(6.5, 6, 4, 3.5) + 0.1, mgp=c(3, 1.2, 0))
    
    #boxplot(numOlaps~bin, data=BP.DF, xlab="Bin position", ylab="Number of overlaps", 
    #        boxwex=0.5, outline=TRUE, cex.axis=1.5, cex=1, col="#FDC776", 
    #        main=paste0("chrALL_", gcb, "_", foi, "_max=", max.val)
    #        )
    #rm(BP.DF, max.val); gc()
    
    #dev.off()
   
    #print(paste0(foi, " done!"), quote=FALSE)
  
  } # itr for loop end
    
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()










