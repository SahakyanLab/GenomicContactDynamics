################################################################################
# Fraction of foi overlapping with Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/7_FeaturePermutation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw_ALL_associated")
foifile = paste0(wk.dir, "/foifile/foifile_TOP2B")
binmx.dir = paste0(wk.dir, "/binmx/out_bindata_1perc_HiCNorm")
out.dir = paste0(wk.dir,"/out_frfoi")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
# Chr allowed for both sets of ranges
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
# Consider bins in either of the Cp
Cp.v = 1:21 
CT.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", 
         "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC",
         "FC", "LC")
# Cell line to be used for cell-type independent features. If allCT.ref="allCT",
# use union of regions from all cell lines.
allCT.ref = "allCT"

# If CpBedFoiChrfilter=TRUE, take out cp regions for chr not in foi bed, this
# should be done if the feature was deliberately not measured in that chr
# CpBedFoiChrfilter=TRUE to be more strict
CpBedFoiChrfilter = FALSE
# Reduce cp bed to combine consecutive regions
CpBedReduce = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(GenomicRanges)
library(regioneR)
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(wk.dir, "/lib/evalfunx.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(CpBedReduce){ print("Reducing Cp bed...", quote=FALSE) }

#
Cp.v <- sort(unique(Cp.v))
CT.v <- sort(CT.v)

print(paste0("Choosing ", allCT.ref, " regions for cell-type independent feature."), 
      quote=FALSE)
if(allCT.ref=="allCT"){ allCT.ref <- CT.v }

# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
ct.v <- unlist(
  lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2])
)

# Prepare contacting region coordinates filtering by Cp.v
if( is.null(chr.v) ){ chr.v <- paste0("chr", c(1:22, "X")) }
BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  BIN.MX <- cbind.data.frame(chr=chr, BIN.MX, stringsAsFactors=FALSE)
})
BIN.MX <- do.call("rbind", BIN.MX); rownames(BIN.MX) <- NULL
CP.bed <- BIN.MX[, c("chr", "start", "end")]
BIN.MX <- data.matrix(BIN.MX[,!colnames(BIN.MX)%in%c("chr","start","end")])
# Remove regions not forming contact in any tissue
rws.TF <- rowSums(BIN.MX)>0
CP.bed <- CP.bed[rws.TF,]; BIN.MX <- BIN.MX[rws.TF,]; rm(rws.TF); gc()

# Association of contacting region and feature bed file

i.len <- length(foi.v)
for(i in 1:i.len){
  
  ct <- ct.v[i]
  if(ct%in%CT.v){
    col.CT <- ct
  } else if(ct=="hg19") {
    col.CT <- allCT.ref
  } else {
    stop("Invalid ct of feature.")
  }
  
  # Foi bed (Default B)
  foi <- foi.v[i]
  foi.bed <- read.table(file=paste0(foi.dir, "/", foi), stringsAsFactors=FALSE, 
                        header=FALSE)[,1:3]
  colnames(foi.bed) <- c("chr", "start", "end")
  foi.bed <- foi.bed[foi.bed$chr%in%chr.v,]
  if(CpBedFoiChrfilter){
    foichr.v <- unique(foi.bed$chr)
    print("Filtering cp bed for chromosomes in foi bed...", quote=FALSE)
  }
  if(nrow(foi.bed)==0){
    next; print(paste0(foi, " skipped because no foi regions."), quote=FALSE)
  }
  foi <- tail(x=strsplit(x=foi, split="\\/")[[1]], n=1)
  pat <- paste(paste0("ct_", ct, "_"), "foi\\_|desc\\_|\\.bed", sep="|")
  foi <- paste0(ct, "_", gsub(x=foi, pattern=pat, replacement="")); rm(pat)
  if( nrow(foi.bed)!=length(reduce(toGRanges(foi.bed))) ){
    print(paste0(foi, ": Not reduced."), quote=FALSE)
  }
  
  B <- toGRanges(foi.bed); rm(foi.bed); gc()
  foinum <- length(B)
  foibp <- sum(width(B))
  
  FOLAP.DF <- list()
  for(Cp in Cp.v){
    
    # Cp bed (Default A)
    col.nme <- paste("s_Cp_", rep(Cp, each=length(col.CT)), 
                     "_ct_", rep(col.CT, times=length(Cp)), 
                     "_e", sep="")
    temp <- as.matrix(BIN.MX[,col.nme])
    cp.bed <- CP.bed[rowSums(temp)>0,]; rm(temp); gc()
    if(nrow(cp.bed)==0){
      next; print(paste0(foi, " skipped because no RS."), quote=FALSE)
    }
    
    #
    if(CpBedFoiChrfilter){
      cp.bed <- cp.bed[cp.bed$chr%in%foichr.v,]
    } 
    
    A <- toGRanges(cp.bed); rm(cp.bed)
    
    if(CpBedReduce){ 
      A <- reduce(A) 
      print( paste0("Max width after reducing is ", max(width(A)), ".") )
    }
    
    #
    FOLAP.DF[[Cp]] <- data.frame(foi=foi, Cp=Cp, 
                                 numOlapB=eval.f.lst$numOlapB(A, B)/foinum,
                                 comOlap=eval.f.lst$comOlap(A, B)/foibp)
    
  } # Cp.v for loop end
  
  FOLAP.DF <- do.call("rbind", FOLAP.DF)

  rm(A, B, out.name); gc()
  print(paste0(foi, " done!"), quote=FALSE)
  
} # i.len for loop end

# rm(list=ls()); gc()

