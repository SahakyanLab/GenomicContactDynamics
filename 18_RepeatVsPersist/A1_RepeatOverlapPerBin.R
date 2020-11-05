################################################################################
# Count insertion sites of repeats in each contact region (any type of overlap).
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "fam" # "fam" | "subfam"
repmaskfile = paste0(data.dir, "/ucsc_tables/hsa_RepeatMasker/RepeatMasker_hg19")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_RepeatOverlapPerBin/", rep.group)
### OTHER SETTINGS #############################################################
hic.resol = 40000L

#2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gcb = "min05Mb" 
chr.v = paste("chr", c(1:22, "X"), sep="")

overlapBy = "repFamily" # repFamily for families | "repName" for sub-families
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
print(paste0(gcb, "..."), quote=FALSE)

for(chr in chr.v){
  
  load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                         unique(PERSIST.MX$hits[,"j"])) )
  rm(PERSIST.MX); gc()
  
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
  
  save(BINREP.MX,  file=paste0(out.dir,"/", gcb, "_", chr, 
                               "_BinRep.RData"))
  
  rm(BINREP.MX, bins.uniq); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}

# rm(list=ls()); gc()
 