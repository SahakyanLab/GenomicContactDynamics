################################################################################
# Get p-value of aAssociation of region sets with features
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw")
foifile = paste0(wk.dir, "/C1_associate/foifile/foifile613")
mask.dir = paste0(wk.dir, "/mask")
binmx.dir = paste0(wk.dir, "/binmx/out_bindata_1perc_HiCNorm")
out.dir = paste0(wk.dir,"/out_associate_Cp21")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
# Chr allowed for both sets of ranges
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
# Consider bins in either of the Cp
Cp.v = 21 
CT.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", 
         "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC",
         "FC", "LC")
# Cell line to be used for cell-type independent features. If allCT.ref="allCT",
# use union of regions from all cell lines.
allCT.ref = "allCT"
Cs.v = c(1, 0.01) # Depends on BIN.MX values; 1-Bin forms contact; 5-Bin part of top 5%
id = "Cp21" # c("Cp21", "CptopCP3", "Cp1", "CpAll", "CpAllCs1perc", "CpAllCsmatch")

# Permutation test parameters
NTIMES = 10000
SEED = 662
nCPU = 3L
# If supplyMask=TRUE, use mask from mask.dir.
supplyMask = FALSE
maxmaskOlapFr = 0
# Background
bgr = "HiCAll" # "HiCAll" | "genome"
randomfunx = "resample" # "circular2" | "random" | "resample"

eval.f.v = c("numOlapA", "numOlapAwithin",
             "numOlapB", "numOlapBwithin",
             "comOlap", "meandist",
             "meanLenOlapB") 

# Local z-score calculation parameters
zs.window=2e6
zs.step=2e4

# If CpBedFoiChrfilter=TRUE, take out cp regions for chr not in foi bed, this
# should be done if the feature was deliberately not measured in that chr
CpBedFoiChrfilter = TRUE
# Reduce cp bed to combine consecutive regions
CpBedReduce = FALSE
# A=feature, B=contact regions. If switchAB=TRUE, switch A and B. 
switchAB = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(regioneR)
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(wk.dir, "/lib/evalfunx.R"))
source(paste0(wk.dir, "/lib/permute_", randomfunx, ".R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Print info about run
out.id <- paste0("nperm", NTIMES, "_seed", SEED, "_mxmskfr", maxmaskOlapFr, "_", id)
print(paste0("Cs value/s: ", paste(Cs.v, collapse=";")), quote=FALSE)
print(paste0("Cp value/s: ", paste(Cp.v, collapse=";")), quote=FALSE)
print(out.id, quote=FALSE)
print(paste0("Running on ", nCPU, " cores..."), quote=FALSE)
if(CpBedReduce){ print("Reducing Cp bed...", quote=FALSE) }
if(switchAB){
  print("Switching A and B. A=feature, B=contact region.", quote=FALSE) 
} else {
  print("A=contact region, B=feature.", quote=FALSE) 
}
print(paste0("Max allowed % overlap with mask: ", maxmaskOlapFr*100), quote=FALSE)
#---------------------------------------
Cp.v <- sort(unique(Cp.v))
CT.v <- sort(CT.v)
Cs.v <- sort(as.numeric(unique(Cs.v)))

print(paste0("Choosing ", allCT.ref, " regions for cell-type independent feature."), 
      quote=FALSE)
if(allCT.ref=="allCT"){ allCT.ref <- CT.v }
#---------------------------------------
# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
ct.v <- unlist(
  lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2])
)
#---------------------------------------
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
#---------------------------------------
# Prepare mask to be applied to genome before association
# During randomisation, no regions will come from these areas
# Masking with blacklist because usually this is also done during ChIP-seq analysis 
# (where most of our feature data come from)
if(supplyMask){
  print("Supply mask...", quote=FALSE)
  path.v <- list.files(path=mask.dir, full.names=TRUE)
  path.v <- path.v[!grepl(x=path.v, pattern="ignore")]
  print(path.v, quote=FALSE)
  mask.bed <- sapply(X=path.v, simplify=FALSE, FUN=function(path){
    bed <- read.table(file=path, stringsAsFactors=FALSE, 
                      header=FALSE, row.names=NULL)[,1:3]
  })
  rm(path.v)
  mask.bed <- do.call("rbind", mask.bed)
  rownames(mask.bed) <- NULL
  mask.bed <- mask.bed[mask.bed[,1]%in%chr.v,]
  # By default, reduce()takes into account strand. 
  mask.bed <- reduce(toGRanges(mask.bed[,1:3]))
} else {
  mask.bed <- GRanges()
  print("Don't supply mask...", quote=FALSE)
}
#---------------------------------------
# Prepare genome
if( randomfunx%in%c("circular2", "random") ){
  genome <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                       colClasses=c("character", "integer", "integer"))[1:2]
  genome <- genome[genome$chromosome%in%chr.v,]
  genome <- cbind.data.frame(chr=genome$chromosome, start=1, end=genome$length.bp)
  if(bgr=="genome"){ print("Background is whole genome...", quote=FALSE) }
} else if(randomfunx=="resample" & bgr=="HiCAll") {
  print("Background is regioneR::resampleRegions of HiCAll...", quote=FALSE)
} else {
  stop("Invalid randomfunx/bgr argument.")
}
#---------------------------------------
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
  #-------------------HiCAll bgr
  if(bgr=="HiCAll"){
    col.nme <- paste("s_Cp_", rep(1:21, each=length(col.CT)), 
                     "_ct_", rep(col.CT, times=length(Cp.v)), 
                     "_e", sep="")
    HiCAll <- as.matrix(BIN.MX[,col.nme])
    HiCAll <- CP.bed[rowSums(HiCAll)>0,]
    if(nrow(HiCAll)==0){
      next; print(paste0(foi, " skipped because no HiCAll."), quote=FALSE)
    }
    if(randomfunx=="resample"){
      genome <- HiCAll
      mask.bed2 <- NULL; rm(mask.bed)
    } else {
      mask.bed1 <- setdiff(toGRanges(genome), toGRanges(HiCAll))
      mask.bed2 <- reduce(union(mask.bed, mask.bed1))
      rm(mask.bed1); gc()
      if( !supplyMask & length(commonRegions(toGRanges(HiCAll), mask.bed2))>0 ){
        stop("Intersection of HiCAll and mask not empty.")
      }
      rm(col.nme, HiCAll); gc()
    }
    print("Background is HiCAll...", quote=FALSE)
  }
  #-------------------Cp bed (Default A)
  col.nme <- paste("s_Cp_", rep(Cp.v, each=length(col.CT)), 
                   "_ct_", rep(col.CT, times=length(Cp.v)), 
                   "_e", sep="")
  temp <- as.matrix(BIN.MX[,col.nme])
  temp[ !temp%in%as.numeric(Cs.v) ] <- 0
  cp.bed <- CP.bed[rowSums(temp)>0,]; rm(temp); gc()
  if(nrow(cp.bed)==0){
    next; print(paste0(foi, " skipped because no RS."), quote=FALSE)
  }
  #-------------------Foi bed (Default B)
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
  # Make sure the bed file is reduced because number of random regions drawn is
  # the same as the number of intervals after reduce()
  foi <- tail(x=strsplit(x=foi, split="\\/")[[1]], n=1)
  pat <- paste(paste0("ct_", ct, "_"), "foi\\_|desc\\_|\\.bed", sep="|")
  foi <- paste0(ct, "_", gsub(x=foi, pattern=pat, replacement="")); rm(pat)
  if( nrow(foi.bed)!=length(reduce(toGRanges(foi.bed))) ){
    print(paste0(foi, ": Not reduced."), quote=FALSE)
  }
  #-------------------Finalise A, B and background regions
  if(CpBedFoiChrfilter){
    cp.bed <- cp.bed[cp.bed$chr%in%foichr.v,]
  } 
  genome <- toGRanges(genome)
  A <- toGRanges(cp.bed); rm(cp.bed)
  B <- toGRanges(foi.bed); rm(foi.bed); gc()
  
  if(switchAB){
    A.temp <- A; B.temp <- B
    A <- B.temp; B <- A.temp
    rm(A.temp, B.temp); gc()
  }
  
  if(CpBedReduce){ 
    A <- reduce(A) 
    print( paste0("Max width after reducing is ", max(width(A)), ".") )
  }
  #-------------------Permutation test
  out.name <- paste0(gcb, "_", foi, "_", out.id)
  permute(out.dir=out.dir, out.name=out.name, SEED=SEED, A=A, B=B, 
          NTIMES=NTIMES, genome=genome, eval.f.lst=eval.f.lst[eval.f.v], 
          mask.bed=mask.bed2, maxmaskOlapFr=maxmaskOlapFr, nCPU=nCPU, 
          # Largest size of consecutive 40-kb cp=21 bins is 1.2 Mb
          zs.window=zs.window, zs.step=zs.step)
  
  if( exists("mask.bed2") ){ rm(mask.bed2) }
  rm(A, B, out.name); gc()
  print(paste0(foi, " done!"), quote=FALSE)
  
} # i.len for loop end

# rm(list=ls()); gc()

