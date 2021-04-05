################################################################################
# Calculate base content per bin with or without masking genome of a feature or 
# combination thereof. If a given sequence (after masking, if applicable) contains 
# at least one "N" witin, all the values in that row will be populated with NA. 
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_Expl_contacts"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
lib.getMaskedKmerCount = paste0(lib, "/getMaskedKmerCount")
binkmer.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_allbins")
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
maskbed.dir = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/regions_anno/out_hierarchyRegions")
# List masks/combinations of masks (;-separated)
maskfile = NULL
out.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binBaseContent/maskfilearr2.repl")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr14"
bin.len = 4e4L
kmer.len = 1L
# Unique HiC bins per chr, parallel execution by getKmerCountsPerInterval()
nCPU = 1L 

genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome."
fastafile.ending = ".fa"

# If seq has fraction of m >= maskLimit, populate with NAs
#maskLimit = 0.5
maskingChar = NULL # "m" | NULL; If maskingChar is NULL, no masking.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(foreach)
library(doParallel)
library(itertools)
library(compiler)
library(RColorBrewer)

source(paste0(lib, "/UTL_doPar.R"))

source(paste0(lib.getMaskedKmerCount, "/getMaskAffix.R"))

source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 

source(paste0(lib.getMaskedKmerCount, "/maskGenome.R"))

source(paste0(lib.TrantoRextr, "/GEN_getKmers.R"))                 
source(paste0(lib.TrantoRextr, "/GEN_getKmerHybridisationGs.R")) 
source(paste0(lib.getMaskedKmerCount, "/getKmerCountsPerIntervalMASK.R"))  

source(paste0(lib.getMaskedKmerCount, "/getMaskedKmerCount.R"))  
source(paste0(lib.getMaskedKmerCount, "/getMergedReducedMask.R"))  
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( !dir.exists(out.dir) ){
  dir.create(path=out.dir)
}

if( !is.null(maskingChar) ){
  # List of masks or combinations
  maskcombi.v <- readLines(con=maskfile)
} else {
  maskcombi.v <- ""
}

# Mask or combination of masks
for(combi in maskcombi.v){
  
  if( !is.null(maskingChar) ){
    
    # Get affix and cell/tissue type from maskfile
    affix <- getMaskAffix(mask=combi, splitChar=";")
    # Reduce bed to remove overlapping ranges
    mergedmask <- getMergedReducedMask(maskbed.dir=maskbed.dir,
                                       mask=combi,
                                       splitChar=";")[,c("chr", "start", "end")]
    
  } else {
    
    affix <- "" 
    mergedmask <- NULL
  
  }
  
  for(chr in chr.v){
    
    if( !is.null(maskingChar) ){
      
      chrmask.TF <- mergedmask[,1]==chr
      
      # Skip if bed file empty or if chr not in mask
      if( sum(chrmask.TF)==0 ){
        print(paste0(combi, " skipped."))
        rm(chrmask.TF)
        next
      }
      
    } else {
      chrmask.TF <- NULL
    }
    
    getMaskedKmerCount(
      
      # Get bins and their coordinates
      binkmer.dir=binkmer.dir,
      out.dir=out.dir,
      
      gcb=gcb,
      chr=chr, 
      kmer.len=kmer.len,
      HiC.res=HiC.res,
      # For identifying outputs;_<affix>
      affix=affix,
      nCPU=nCPU,
      
      PATH.genome=genome.dir,
      genome.prefix=genome.prefix,
      fastafile.ending = ".fa",
      split = FALSE, 
      
      maskbed=mergedmask[chrmask.TF,],
      maskingChar=maskingChar
      
    ) 
    
    rm(chrmask.TF)
    
  } # chr.v for loop end
  
  print(paste0(chr, " done!"), quote=FALSE)
  
} # maskcombi.v for loop end

# rm(list=ls()); gc()
