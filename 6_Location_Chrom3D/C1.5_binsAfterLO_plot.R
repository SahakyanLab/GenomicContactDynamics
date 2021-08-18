################################################################################
# Make plot showing the length of contact bins after hg19 -> hg38 liftover
# and the percentage of bins per chromosome that were split into multiple ranges
# after liftover.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/6_Location_Chrom3D"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
model.id = "H1-hESC_LMNB1_hg38"
LOdata.dir  = paste0(wk.dir, "/out_mapToHiCcontactBins_Chrom3D/H1-hESC_LMNB1_hg38") 
out.dir  = paste0(wk.dir, "/out_binsAfterLO_plot") 
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt") 
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
bin.len = 40000
chr.v = paste0("chr", c(1:22, "X"))
LOchain = "hg19ToHg38"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/plotLengthDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste(gcb, LOchain, model.id, sep="_")

chrlen.df <- read.delim(file=chrlen.file, header=T, sep="\t", stringsAsFactors=F)
totbin.v <- setNames(object=chrlen.df$bins.40kb, nm=chrlen.df$chromosome)

if( !identical(as.numeric(totbin.v), as.numeric(chrlen.df$bins.40kb)) ){
  stop("Checkpoint 1.")
}

LO <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  
  load(paste0(LOdata.dir, "/", chr, "_Persist_min2Mb_", LOchain, "_RAW.RData"))
  LO.mx$seqnames <- as.character(LO.mx$seqnames)
    
  LO.mx <- LO.mx[!(duplicated(x=LO.mx$group, fromLast=F) | duplicated(x=LO.mx$group, fromLast=T)),
                 c("group", "seqnames", "width")]
  
  if( any(duplicated(LO.mx$group)) | any(LO.mx$group>as.numeric(totbin.v[chr])) ){
    stop(paste0(chr, ": Checkpoint 2."))
  }
  
  return(LO.mx)
  
})

LO <- do.call("rbind", LO)
rownames(LO) <- NULL

LO$group <- as.numeric(as.character(LO$group))
LO$width <- as.numeric(as.character(LO$width))

if( any(!chr.v%in%LO$seqnames) | any(!LO$seqnames%in%chr.v) ){
  stop("Checkpoint 3.")
}
  
# Percentage of bins split into multiple ranges after liftover
binPerChr.v <- table(LO$seqnames)
dropout <- 100*(1-binPerChr.v/totbin.v[names(binPerChr.v)])
write.csv(x=cbind.data.frame(binValidPerChr=as.numeric(binPerChr.v), stack(dropout)), 
          file=paste0(out.dir, "/", out.name, "_percentBinWithMultipleRangesAfterLO.csv"),
          row.names=F)

print(range(LO$width, na.rm=T))
#[1]   324 40000

# Plot total length of converted ranges per bin   
plotLengthDist(df=cbind.data.frame(variable="LO", value=LO$width/10000),
               vline.v=NULL,
               col.v="deepskyblue3",
               out.name=paste0(out.name, "_widthDensAfterLO.pdf"),
               out.dir=out.dir,
               label.x="Length/10000")

# rm(list=ls()); gc()