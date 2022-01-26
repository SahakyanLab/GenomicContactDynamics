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
out.dir = paste0(wk.dir, "/out_FETA_sharedisofpcoding")
# feta.dir is FETA.MX directory
feta.dir = paste0(out.dir, "/temp")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced")
# If foifile = NULL, all files in foi.dir
foifile = paste0(wk.dir, "/foifile/foifile_sharedisofpcoding")
#foigroup = "topo" #"HMnp"
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
# foi
nCPU = 1L # stick to 5 NCPU, ~300G
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(reshape)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)

foi.v.len <- length(foi.v)
toExport <- c("foi.v", "chr.v.len", "chr.v", "persist.dir", "gcb", 
              "feta.dir", "out.dir")

#png(file=paste0(out.dir, "/chrALL_", gcb, "_", foigroup, "_BP.png"),
#    height=6000, width=10000)
#par(oma=c(0,0,0,0), mar=c(6.5, 6, 4, 3.5) + 0.1, mgp=c(3, 1.2, 0), mfrow=c(4,8))

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
    # Boxplot of overlap values (absolute)
    BP.DF <- melt(FETA.MX[,!colnames(FETA.MX)%in%c("chr", "bin")])
    rm(FETA.MX); gc()
    colnames(BP.DF) <- c("bin", "numOlaps")
    max.val <- max(BP.DF$numOlaps, na.rm=TRUE)
    
    bitmap(file=paste0(out.dir, "/max", max.val, "_chrALL_", gcb, "_", foi, ".png"), 
           height=10, width=10, units="in", res=500, type="png256")
    par(oma=c(0,0,0,0), mar=c(6.5, 6, 4, 3.5) + 0.1, mgp=c(3, 1.2, 0))
    
    boxplot(numOlaps~bin, data=BP.DF, xlab="Bin position", ylab="Number of overlaps", 
            boxwex=0.5, outline=TRUE, cex.axis=1.5, cex=1, col="#FDC776", 
            main=paste0("chrALL_", gcb, "_", foi, "_max=", max.val)
            )
    rm(BP.DF, max.val); gc()
    
    dev.off()
   
    print(paste0(foi, " boxplot done!"), quote=FALSE)
  
  } # itr for loop end
    
}
### END OF PARALLEL EXECUTION ###

#dev.off()

# rm(list=ls())










