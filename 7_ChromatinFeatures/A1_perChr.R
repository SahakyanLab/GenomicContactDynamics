################################################################################
# Define Xkb-windows spanning each of the unique contacting bins per chr and 
# check the overlap of chromatin features (supplied as bed file) at Y-positions 
# in the window. 
# For making metaplot of chromatin features 
# Note: Counts more than one are not converted here. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
#options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/10_ChromatinFeatures"
    data.dir = "/Users/ltamon/Database"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/10_ChromatinFeatures"
    data.dir = "/t1-data/user/ltamon/Database"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/10_ChromatinFeatures"
    data.dir = "/home/ltamon/Database"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# PERSIST.MX directory
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
# Chromatin features directory
# Remove column names in bed file
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced_b2b3")
# List of filenames of features of interest (refer to foi.dir)
# If foifile = NULL, all files in foi.dir
foifile = paste0(wk.dir, "/foifile/foifile_ABscomp")
# Number of bins
nCPU = 4L 
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
out.dir = paste0(wk.dir, "/out_FETA_b2b3")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
HiC.res = 4e4L
pos.v = -12:12
chr.v = paste("chr", c(1:22, "X"), sep="")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(IRanges)
library(GenomicRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(wk.dir, "/lib/olapPerGenomicPos.R"))
source(paste0(lib, "/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  # Load PERSIST.MX
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  ubins <- unique( c(unique(PERSIST.MX$hits$i), unique(PERSIST.MX$hits$j)) )
  data.frame(chr=chr, bin=ubins, stringsAsFactors=FALSE, row.names=NULL)
})
BIN.MX <- do.call("rbind", BIN.MX)
rownames(BIN.MX) <- NULL
#-------------------------------------------------------------------------------
# Chromosome length file
chrLen.df <- fread(file=chrLenfile, colClasses=list(character=1, integer=2), 
                   stringsAsFactors=FALSE, data.table=FALSE)
# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)

for(foi in foi.v){
  
  foi.bed <- fread(file=paste0(foi.dir, "/", foi),
                   colClasses=list(character=1, integer=2:3), 
                   stringsAsFactors=FALSE, data.table=FALSE, header=FALSE)
  foi <- tail(x=strsplit(x=foi, split="\\/")[[1]], n=1)
  foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
  
  chr.v <- chr.v[chr.v%in%intersect(paste("chr", c(1:22, "X"), sep=""), unique(foi.bed[,1]))]
  if( length(chr.v)==0 ){ "No chr to interrogate." }
  pos.v.len <- length(pos.v)
  
  for(chr in chr.v){
    
    ubins <- BIN.MX[BIN.MX$chr==chr,"bin"]
    bin.end <- ubins*HiC.res
    bin.start <- bin.end-HiC.res+1
    
    # Subset bed file based on chr
    #bed <- foi.bed[foi.bed[,1]==chr,]
    chr.TF <- foi.bed[,1]==chr
    
    if(sum(chr.TF)==0L){
      print(paste0(chr, " skipped"), quote=FALSE)
      next
    } else {
      # Get only bins (extended) that overlap with foi.bed (to make it faster for 
      # smaller bed files)
      # Query = extended HiC contact bins
      # Subject = feature 
      olap <- WhichOverlap(start.query=bin.start+(HiC.res*min(pos.v)), 
                           end.query=bin.end+(HiC.res*max(pos.v)), 
                           space.query=rep("a", length(ubins)),
                           start.subject=foi.bed[chr.TF,2], 
                           end.subject=foi.bed[chr.TF,3], 
                           space.subject=rep("a", sum(chr.TF)),
                           maxgap=-1L, minoverlap=1L,
                           type="any")
      ubins.ind <- unique(olap[,1])
      rm(olap); gc()
      
      # Initialize output matrix
      FETA.MX <- matrix(data=0, nrow=length(ubins), ncol=pos.v.len)
      dimnames(FETA.MX)[[2]] <- as.character(format(x=pos.v, scipen=FALSE))
      
      minPos <- 1L
      maxPos <- as.integer(chrLen.df[chrLen.df[,"chromosome"]==chr, "length.bp"])
      
      if(length(ubins.ind)>0L){
        mx <- olapPerGenomicPos(
          # Feature of interest
          foi.id=paste0(chr, "_", foi),
          foi.df=foi.bed[chr.TF,],
          start=bin.start[ubins.ind],
          end=bin.end[ubins.ind],
          res=HiC.res,
          pos.v=pos.v,
          nCPU=nCPU,
          minPos=minPos,
          maxPos=maxPos,
          # Overlap parameters
          min.olap=1L,
          # This default is fixed regardless of os for R.3.5.0 and up
          max.gap=-1L, 
          type.olap="any"
        )
        # Fill bins with overlap to bed file  
        FETA.MX[ubins.ind,] <- mx; rm(mx)
      }
      # Fill the bins with no overlaps to bed file
      for( ind in setdiff(1:length(ubins),ubins.ind) ){
        start.actual <- bin.start[ind]+(pos.v*HiC.res)
        end.actual <- bin.end[ind]+(pos.v*HiC.res)
        incl.TF <- !(start.actual<minPos | start.actual>maxPos | end.actual<minPos | end.actual>maxPos)
        
        foicount <- rep(NA, times=pos.v.len)
        foicount[incl.TF] <- 0L
        
        FETA.MX[ind,] <- foicount
        rm(start.actual, end.actual, foicount)
      }
      FETA.MX <- cbind(bin=ubins, FETA.MX)
      rm(ubins, ubins.ind, bin.start, bin.end); gc()
      save(FETA.MX, file=paste0(out.dir, "/temp/", chr, "_", gcb, "_", foi, ".RData"))
      rm(FETA.MX)
    }
    rm(chr.TF); gc()
  } # chr.v for loop end
  
  rm(foi.bed); gc()
  print(paste0(foi, " done!"), quote=FALSE)
  
} # foi.v for loop end

# rm(list=ls())
