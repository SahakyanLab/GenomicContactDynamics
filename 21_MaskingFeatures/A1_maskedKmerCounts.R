################################################################################
# Mask a or any combination of features in genome and then re-calculate k-mer 
# counts per unique contacting bin (BINKMER.MX). Bins with missing sequence (N)
# populated with NAs unless masking character is N.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    lib = paste0(home.dir, "/DPhil/lib")
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/21_MaskingFeatures")
    data.dir = paste0(home.dir, "/Database")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    lib = paste0(home.dir, "/DPhil/lib")
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/12_MaskingFeatures")
    data.dir = paste0(home.dir, "/Database")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# TrantoRextr functions directory
lib.TrantoRextr = paste0(lib, "/TrantoRextr")
binkmerAll.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_allbins") #binkmer_divLen_all")
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
# Mask bed files directory
maskbed.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced")
# List masks/combinations of masks (;-separated)
maskfile = paste0(wk.dir, "/A1_masking/maskfile/maskfile101")
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

source(paste0(lib, "/UTL_doPar.R"))

source(paste0(lib, "/getMaskedKmerCount/getMaskAffix.R"))

source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 

source(paste0(lib, "/getMaskedKmerCount/maskGenome.R"))

source(paste0(lib.TrantoRextr, "/GEN_getKmers.R"))                 
source(paste0(lib.TrantoRextr, "/GEN_getKmerHybridisationGs.R")) 
source(paste0(lib, "/getMaskedKmerCount/getKmerCountsPerIntervalMASK.R"))  

source(paste0(lib, "/getMaskedKmerCount/getMaskedKmerCount.R"))  
source(paste0(lib, "/getMaskedKmerCount/getMergedReducedMask.R")) 
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# List of masks or combinations
maskcombi.v <- readLines(con=maskfile)

# Mask or combination of masks
for(combi in maskcombi.v){

  # Get affix and cell/tissue type from maskfile
  affix <- getMaskAffix(mask=combi, splitChar=";")
  
  # Reduce bed to remove overlapping ranges
  mergedmask <- getMergedReducedMask(maskbed.dir=maskbed.dir,
                                     mask=combi,
                                     splitChar=";")[,c("chr", "start", "end")]
  
  # Skip if bed file empty or if chr not in mask
  if( length(mergedmask$chr) == 0 ){
    
    rm(mergedmask)
    print(paste0(combi, " mask empty.")); next
    
  }
  
  for(chr in chr.v){
    
    chrmask.TF <- mergedmask[,1]==chr
    
    # Skip if bed file empty or if chr not in mask
    if( sum(chrmask.TF)==0 ){
      print(paste0(combi, "_", chr, ": skipped.")); next
    } else {
      
      getMaskedKmerCount(
        
        # PERSIST.MX directory
        binkmer.dir=binkmerAll.dir,
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
      
    }
    
    print(paste0(chr, " done!"), quote=FALSE)
    
  } # chr.v for loop end
  
  print(paste0(combi, " done!"), quote=FALSE)
  
} # maskcombi.v for loop end

# rm(list=ls()); gc()
