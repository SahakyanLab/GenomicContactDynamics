################################################################################
# Count insertion sites of repeats in ALL bins of a chromosome 
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/18_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam" # "fam" | "subfam"
repmaskfile = paste0(data.dir, "/ucsc_tables/hsa_RepeatMasker/RepeatMasker_hg19")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
out.dir = paste0(wk.dir, "/out_RepeatOverlapPerBinALL/", rep.group)
### OTHER SETTINGS #############################################################
hic.resol = 40000L
chr.v = paste("chr", c(1:22, "X"), sep="")
overlapBy = "repName" # repFamily for families | "repName" for sub-families
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel) 
library(itertools) 
library(GenomicRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/TrantoRextr/GEN_getRepMaskSubset.R"))
source(paste0(wk.dir, "/lib/GEN_getBinsOverlapRepMasker2.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrLen.df <- read.delim(file=chrLenfile, header=TRUE)

for(chr in chr.v){
  
  bin.num <- as.integer(chrLen.df[chrLen.df[,"chromosome"]==chr, "length.bp"])
  bin.num <- ceiling(bin.num/hic.resol)
  if( bin.num!=chrLen.df$bins.40kb[chrLen.df$chromosome==chr]  ){
    stop(paste0(chr, ": Checkpoint 1."))
  }
  bins.uniq <- 1:bin.num
 
  BINREP.MX <- getBinsOverlapRepMasker(
    chr=chr,
    hic.resol=hic.resol,
    bins=bins.uniq,
    repmask.Filepath=repmaskfile,
    overlapBy=overlapBy,
    min.olap = 1L,
    max.gap = -1L,
    type.olap = "any"
  )
  
  save(BINREP.MX,  file=paste0(out.dir,"/", chr, "_BinRep.RData"))
  
  rm(BINREP.MX, bins.uniq, bin.num); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}

# rm(list=ls()); gc()
 