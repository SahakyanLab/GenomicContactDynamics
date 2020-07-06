################################################################################
# Calculate how much of the genome is covered by bins per Cp
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
    binmx.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Circos/out_bindata"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Circos/out_bindata"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir,"/out_associate_runB_topo_genesLTr")
mask.dir = paste0(wk.dir, "/mask")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced_temp")
foifile = paste0(wk.dir, "/foifile/foifile_topo_genesLTr")
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
# If chr.v = NULL, chr.v = chr1:22, chrX
chr.v = NULL
# Consider bins in either of the Cp
Cp.v = REPLACE1
cp.id = "REPLACE2"

# Permutation
NTIMES = 10000
SEED = 342
nCPU = 5
masking = FALSE
# if length(A)<min.parallel - single thread
MIN.PARALLEL = 1000
eval.f.v <- c("numOlapA", "numOlapAwithin", 
              "numOlapB", "numOlapBwithin",
              #"comOlap", #"comOlapA", "comOlapB", 
              "meandist"
              #, "meanLenOlapB"
              ) 

## Local z-score calculation. If NULL, default values are used. 
#zs.window = 1e6
#zs.step = 10000

# If CpBedChrfilter=TRUE, take out cp regions for chr not in foi bed, this
# should be done if the feature was deliberately not measured in that chr
CpBedChrfilter = FALSE
# Reduce cp bed to combine consecutive regions
CpBedReduce = FALSE
# A=feature, B=contact regions. If switchAB=TRUE, switch A and B. 
switchAB = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(regioneR)
source(paste0(lib, "/finaliseFOI.R"))

# Custom evaluation functions
numOlapA <- function(A, B, ...) {
  num <- regioneR::numOverlaps(A=A, B=B, count.once=TRUE)
  return( num/length(A)*100 )
}
numOlapAwithin <- function(A, B, ...) {
  num <- length(
    subsetByOverlaps(x=A, ranges=B, maxgap=-1L, minoverlap=1L, type="within",
                     invert=FALSE)
  )
  return( num/length(A)*100 )
}
numOlapB <- function(A, B, ...) {
  num <- regioneR::numOverlaps(A=B, B=A, count.once=TRUE)
  return( num/length(B)*100 )
}
numOlapBwithin <- function(A, B, ...) {
  num <- length(
    subsetByOverlaps(x=B, ranges=A, maxgap=-1L, minoverlap=1L, type="within",
                     invert=FALSE)
  )
  return( num/length(B)*100 )
}
comOlap <- function(A, B, ...) {
  return(sum( width(commonRegions(A,B)) )
         )
}
comOlapA <- function(A, B, ...) {
  total <- sum( width(A) )
  common <- sum( width(commonRegions(A,B)) )
  return( common/total*100 )
}
comOlapB <- function(A, B, ...) {
  total <- sum( width(B) )
  common <- sum( width(commonRegions(A,B)) )
  return( common/total*100 )
}
meandist <- function(A, B, ...) {
  return( -(regioneR::meanDistance(A,B)) )
}
meanLenOlapA <- function(A, B, ...) {
  x <- mean(width(
    subsetByOverlaps(x=A, ranges=B, maxgap=-1L, minoverlap=1L, type="any", 
                     invert=FALSE)
  ))
  return(x)
}
meanLenOlapB <- function(A, B, ...) {
  x <- mean(width(
    subsetByOverlaps(x=B, ranges=A, maxgap=-1L, minoverlap=1L, type="any", 
                     invert=FALSE)
  ))
  return(x)
}
#---------------------------------------
# Selected evaluation functions
ind <- match(eval.f.v, c("numOlapA", "numOlapAwithin", 
                         "numOlapB", "numOlapBwithin",
                         "comOlap", "comOlapA", "comOlapB", 
                         "meandist", "meanLenOlapA", "meanLenOlapB"))
eval.f.lst <- list(`% number of A regions overlapping`=numOlapA,
                   `% number of A regions overlapping within`=numOlapAwithin,
                   `% number of B regions overlapping`=numOlapB,
                   `% number of B regions overlapping within`=numOlapBwithin,
                   `total length of intersection`=comOlap, 
                   `% intersection relative to A`=comOlapA,
                   `% intersection relative to B`=comOlapB,
                   `mean distance A from B`=meandist,
                   `mean length of A regions overlapping`=meanLenOlapA,
                   `mean length of B regions overlapping`=meanLenOlapB
                   )
eval.f.lst <- eval.f.lst[ind[!is.na(ind)]]
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0("Running on ", nCPU, " cores..."), quote=FALSE)
out.id <- paste0("nperm", NTIMES, "_cp", cp.id, "_seed", SEED)
# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
ct.v <- unique(unlist(
  lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2])
))
#---------------------------------------
# Prepare contacting region coordinates
if( is.null(chr.v) ){ chr.v <- paste("chr", c(1:22, "X"), sep="") }
BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  BIN.MX <- cbind.data.frame(chr=chr, 
                             BIN.MX[,colnames(BIN.MX)%in%c("start", "end", ct.v, Cp.v)], 
                             stringsAsFactors=FALSE)
  # Select regions with Cp in Cp.v (union)
  BIN.MX <- BIN.MX[rowSums(as.matrix(BIN.MX[,as.character(Cp.v)]))>0,]
})
BIN.MX <- do.call("rbind", BIN.MX)
rownames(BIN.MX) <- NULL
#---------------------------------------
# Prepare mask to be applied to genome before association
# During randomisation, no regions will come from these areas
# Masking with blacklist because usually this is also done during ChIP-seq analysis 
# (where most of our feature data come from)
if(masking){
  print("Masking...", quote=FALSE)
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
  mask.bed <- as.data.frame(reduce(GRanges(seqnames=mask.bed$V1, 
                                        IRanges(start=mask.bed$V2, end=mask.bed$V3))
  ))[,1:3]
  mask.bed[,1] <- as.character(mask.bed[,1])
} else {
  mask.bed <- NULL
  print("No masking...", quote=FALSE)
}
#---------------------------------------
# Prepare genome
genome <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                     colClasses=c("character", "integer", "integer"))[1:2]
genome <- genome[genome$chromosome%in%chr.v,]
genome <- GRanges(cbind.data.frame(chr=genome$chromosome, start=1, end=genome$length.bp))
#---------------------------------------
if(CpBedReduce){
  print("Reducing Cp bed...", quote=FALSE)
}
if(switchAB){
  print("Switching A and B...", quote=FALSE)
}
# Association of contacting region and feature bed file
print("Association...", quote=FALSE)
for(ct in ct.v){
  cp.bed <- BIN.MX[, c("chr", "start", "end")]
  if(ct%in%c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", 
             "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC",
             "FC", "LC")){
    cp.bed <- cp.bed[BIN.MX[[ct]]==1,]
  }
  foi.ct.ind <- grep(x=foi.v, pattern=paste0("ct_", ct, "_foi"), fixed=TRUE)
  for(ind in foi.ct.ind){
    foi <- foi.v[ind]
    foi.bed <- read.table(file=paste0(foi.dir, "/", foi), stringsAsFactors=FALSE, 
                          header=FALSE)[,1:3]
    colnames(foi.bed) <- c("chr", "start", "end")
    foi.bed <- foi.bed[foi.bed[,1]%in%chr.v,]
    if(CpBedChrfilter){
      foichr.v <- unique(foi.bed$chr)
      print("Filtering cp bed for chromosomes in foi bed...", quote=FALSE)
    }
    if(nrow(foi.bed)==0){ next; print(paste0(foi, " skipped."), quoet=FALSE) }
    
    # Make sure the bed file is reduced because number of random regions drawn is
    # the same as the number of intervals after reduce()
    foi <- tail(x=strsplit(x=foi, split="\\/")[[1]], n=1)
    foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
    if( nrow(foi.bed)!=length(reduce(GRanges(foi.bed))) ){
      stop(paste0(foi, ": Not reduced."), quote=FALSE)
    }
   
    # Do A regions lie on B regions in the genome?
    # Since size range of most features are below the Hi-C resolution.
    # A -> feature B -> contacting regions
    A <- toGRanges(foi.bed); rm(foi.bed); gc()
    if(CpBedChrfilter){
      B <- toGRanges(cp.bed[cp.bed$chr%in%foichr.v,])
    } else {
      B <- toGRanges(cp.bed)
    }
    
    if(CpBedReduce){
      B <- reduce(B)
    }
    
    if(switchAB){
      A.temp <- A
      B.temp <- B
      A <- B.temp
      B <- A.temp
      rm(A.temp, B.temp); gc()
    }
    
    set.seed(SEED)
    PERMT <- permTest(A=A, B=B, ntimes=NTIMES, alternative="auto", 
                      verbose=FALSE, genome=genome, 
                      evaluate.function=eval.f.lst,
                      # Randomisation function
                      randomize.function=circularRandomizeRegions, 
                      per.chromosome=TRUE, max.mask.overlap=1, 
                      # None of the ranges overlap after using reduce() from IRanges
                      non.overlapping=TRUE, 
                      mask=mask.bed,
                      # Parallel processing; if length(A)<min.parallel - single thread
                      min.parallel=MIN.PARALLEL, 
                      # For reproducibility
                      mc.set.seed=FALSE, mc.cores=nCPU)
    
    LZSCOR <- localZScore(A=A, B=B, PERMT)
    #LZSCOR <- localZScore(A=A, B=B, PERMT, window=1.5e6, step=10000)
    rm(A, B); gc()
    if( !( identical(names(PERMT), names(LZSCOR)) & identical(names(PERMT), names(eval.f.lst)) ) ){
      stop("PERMT, LZSCOR and eval.f.lst have different orders.")
    }
    
    # Evalusation function that gave zscore=NaN
    ind <- sapply(X=names(LZSCOR), simplify=TRUE, USE.NAMES=FALSE, FUN=function(eval.f){
      PERMT[[eval.f]][["zscore"]]
    })
    ind <- which(!is.finite(ind))
    if(length(ind)!=0){
      print(paste0(foi, ":", length(ind), " NaN/NA zscore/s."), quote=FALSE)
      # Remove NaN/NAs in LZSCOR to avoid error and termination of script
      for(i in ind){
        len <- length( LZSCOR[[i]][["shifted.z.scores"]] )
        LZSCOR[[i]][["shifted.z.scores"]] <- rep(100, times=21)
      }; rm(ind)
    }
    
    pdf(file=paste0(out.dir, "/", gcb, "_", foi, "_", out.id, "_plots.pdf"), width=10, height=10)
    plot(PERMT)
    plot(LZSCOR)
    dev.off()
    
    save(PERMT, file=paste0(out.dir, "/", gcb, "_", foi, "_", out.id, "_permtest.RData"))
    save(LZSCOR, file=paste0(out.dir, "/", gcb, "_", foi, "_", out.id, "_zscore.RData"))
    rm(PERMT, LZSCOR); gc()
    
    print(paste0(foi, " done!"), quote=FALSE)
  }
  rm(cp.bed); gc()
}

# rm(list=ls())


