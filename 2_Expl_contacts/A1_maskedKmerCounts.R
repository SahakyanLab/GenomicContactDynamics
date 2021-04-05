################################################################################
# Mask a or any combination of features in genome and then re-calculate k-mer 
# counts per unique contacting bin (BINKMER.MX)
# Bins exceeding the limit of acceptable fraction of masked characters and those
# not forming contacts in given cell/tissue are populated by NAs
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/12_MaskingFeatures"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/12_MaskingFeatures"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# TrantoRextr functions directory
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
binkmerAll.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_divLen_all")
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
# Mask bed files directory
maskbed.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced")
# List masks/combinations of masks (;-separated)
maskfile = paste0(wk.dir, "/A1_masking/maskfile/maskfileREPLACE")
#maskfile = paste0(wk.dir, "/A1_masking/test/maskfile_test")
out.dir = paste0(wk.dir, "/out_binkmer_divLen")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr1"
kmer.len = 7L
HiC.res = 4e4L
# Unique HiC bins per chr ( getKmerCountsPerInterval() )
nCPU = 2L # used 5L, ~100G

# If seq has fraction of m >= maskLimit, populate with NAs
#maskLimit = 0.5
maskingChar = "m"

genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome."
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(IRanges)
library(Biostrings)
library(foreach)
library(doParallel)
library(itertools)
library(compiler)
library(RColorBrewer)

source(paste0(wk.dir, "/lib/getMaskAffix.R"))

source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 

source(paste0(wk.dir, "/lib/maskGenome.R"))

source(paste0(lib.TrantoRextr, "/GEN_getKmers.R"))                 
source(paste0(lib.TrantoRextr, "/GEN_getKmerHybridisationGs.R")) 
source(paste0(wk.dir, "/lib/getKmerCountsPerIntervalMASK.R"))  

source(paste0(wk.dir, "/lib/getMaskedKmerCount.R"))  
source(paste0(wk.dir, "/lib/getMergedReducedMask.R"))  
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# List of masks or combinations
maskcombi.v <- readLines(con=maskfile)

# Mask or combination of masks
for(chr in chr.v){

  for(combi in maskcombi.v){
    
    # Get affix and cell/tissue type from maskfile
    affix <- getMaskAffix(mask=combi, splitChar=";")
    
    # Reduce bed to remove overlapping ranges
    mergedmask <- getMergedReducedMask(maskbed.dir=maskbed.dir,
                                       mask=combi,
                                       splitChar=";")
    
    # Skip if bed file empty or if chr not in mask
    if( !chr%in%names(mergedmask) ){
      print(paste0(combi, " skipped.")); next
    }
    
    getMaskedKmerCount(
      
      # PERSIST.MX directory
      binkmerAll.dir=binkmerAll.dir,
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
      
      maskbed=mergedmask[[chr]],
      maskingChar=maskingChar
    )
    
  } # maskcombi.v for loop end
  
  print(paste0(chr, " done!"), quote=FALSE)
  
} # chr.v for loop end

# rm(list=ls())
