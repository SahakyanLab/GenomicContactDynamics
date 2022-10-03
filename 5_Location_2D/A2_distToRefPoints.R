################################################################################
# Distance of bin midpoints to reference points in the genome marking certain 
# areas i.e. chr ends, centromere midpoints. Note that R rounds up >= XXXXXXX.5 
# values to nearest integer when printing which could also happen when saving.
# If last bin shorter than bin resolute, populate output matrix row with NAs.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/5_Location_2D")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")

out.dir  = paste0(wk.dir, "/out_distToRefPoints")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt") 
# 1-based, constant centromere length of 3Mb
centro.file = paste0(data.dir, "/ucsc_tables/hsa_centromere/ct_hg19_foi_centromereonly_desc_DNA")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr21" #paste0("chr", c(1:22, "X"))
nCPU = 1
bin.len = 40000
# Should match ref.points order
ref.points.name = c("chr.start.bp", "chr.end.bp", "centromere.midP.bp") 
options(digits=15) # Digits beyond 15 could be unreliable
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
#################################################################################
# Build list of reference points

centro.df <- read.delim(file=centro.file, header=F, sep="\t", stringsAsFactors=F)
centro.len <- unique( centro.df$V3 - centro.df$V2 + 1 ) # = 3e+06 bp
centro.points <- centro.df$V2 + ((centro.len - 1) / 2) # R rounding off
centro.points <- as.list(centro.points)
names(centro.points) <- centro.df$V1

chrlen.df <- read.delim(file=chrlen.file, header=T, sep="\t", stringsAsFactors=F)
chrlen.df <- chrlen.df[chrlen.df$chromosome != "chrMT",]
chrlen.points <- sapply(names(centro.points), simplify=F, USE.NAMES=T, FUN=function(chr){
  c(1, chrlen.df$length.bp[chrlen.df$chromosome == chr])
})

ref.points <- Map(chrlen.points[names(centro.points)], centro.points[names(centro.points)], f=append)

#

ref.points.id <- paste(ref.points.name, collapse="_")
#chr.v <- unique(centro.df$V1)
chr.v.len <- length(chr.v)
getMidPbin <- (bin.len - 1) / 2 
toExport <- c("chr.v", "chrlen.df", "bin.len", "getMidPbin", "ref.points", "ref.points.name",
              "persist.dir", "gcb", "ref.points.id")

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), .inorder=F,
        .export=toExport, .noexport=ls()[!ls()%in%toExport]

) %op% {
  
  #chunk <- sapply(X=itr, simplify=F, FUN=function(i){

  for(i in itr){
    
    chr <- chr.v[i]
    chr.end.bp <- chrlen.df$length.bp[chrlen.df$chromosome == chr]
    chr.bin.end <- ceiling(chr.end.bp / bin.len)
    ubins <- 1:chr.bin.end
    bin.start <- (ubins * bin.len) - bin.len + 1
    bin.midP <- bin.start + getMidPbin
    
    disttoref.mx <- matrix(data=NA, nrow=length(ubins), ncol=length(ref.points[[1]]), 
                           dimnames=list(NULL, ref.points.name))
    for(bin.num in ubins){
      disttoref.mx[bin.num, ] <- bin.midP[bin.num] - ref.points[[chr]] 
    }
    
    # If last bin shorter than bin resolute, populate with NAs.
    if( abs(bin.start[chr.bin.end] - chr.end.bp + 1) < bin.len ){
      disttoref.mx[chr.bin.end,] <- NA
      print(paste0(chr, ": Populating last bin with NAs in disttoref.mx."))
    }
    
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
    rm(PERSIST.MX)
    
    IJDISTTOREF.MX <- rbind( disttoref.mx[ ij.mx[,"i"], ], 
                              disttoref.mx[ ij.mx[,"j"], ] )
    IJDISTTOREF.MX <- cbind( IJDISTTOREF.MX, Cp=c(ij.mx[,"Cp"], ij.mx[,"Cp"]) )
    
    save(IJDISTTOREF.MX, file=paste0(out.dir, "/", gcb, "_", chr, "_ijdisttoref_", ref.points.id, ".RData"))
    
    rm(disttoref.mx, IJDISTTOREF.MX, ij.mx)
    gc()
    
  }

}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()