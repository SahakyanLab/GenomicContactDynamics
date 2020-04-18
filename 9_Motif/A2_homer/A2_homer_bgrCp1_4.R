################################################################################
# Motif analyses with homer via R API package, marge
# https://github.com/robertamezquita/marge/blob/master/R/find_motifs_genome.R
# Info about homer motif enrichment analysis and the findMotifsGenome.pl command
## http://homer.ucsd.edu/homer/ngs/peakMotifs.html
## http://homer.ucsd.edu/homer/motif/index.html 
## http://homer.ucsd.edu/homer/motif/practicalTips.html
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
version

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/15_Motif"
    data.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/15_Motif"
    data.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# UNIQBIN.DF directory
data.dir = data.dir
out.dir = paste0(wk.dir, "/out_motif")
### OTHER SETTINGS #############################################################
HiCres = 40000L
out.name = "chrALL_uniquebins"
gcb = "min2Mb"
nCPU = 4L
targCp = 21
bgrCp = 1 # "HiCAll" (all HiC contacting bins)
samp.size = 0.8 #use 50% of target/bgr sequences
SEED = 845
options('homer_path' = "/t1-data/user/ltamon/homer")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(marge)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
suffix <- paste0(out.name, "_targCp", targCp, "_bgrCp", bgrCp)
out.dir <- paste0(out.dir, "/", suffix, "_seed", SEED, "_sampsize", samp.size)
if(!dir.exists(out.dir)){
  dir.create(path=out.dir)
}

v <-  as.character(c(bgrCp, targCp))
names(v) <- c("bgr", "targ")
# Load UNIQBIN.DF, dataframe of all Hi-C contacting bins
load(file=paste0(data.dir, "/chrALL_Uniqbin_", gcb, ".RData"))
# Remove bins with Ns because they were not considered in discordance measurement
# min2Mb, 9114/1489257 (0.6%) bins with Ns

UNIQBIN.DF <- UNIQBIN.DF[is.na(UNIQBIN.DF$Ns),]

set.seed(SEED)
for( type in names(v) ){
  
  subset <- v[type]
  
  if(subset=="HiCAll"){
    log <- rep(TRUE, length(UNIQBIN.DF$cp))
  } else if( subset%in%as.character(1:21) ){
    log <- UNIQBIN.DF$cp==as.numeric(subset)
  }
  
  bin.end <- UNIQBIN.DF$bin[log]*HiCres
  ind <- sample(x=1:sum(log), size=sum(log)*samp.size, replace=FALSE)
  bed <- data.frame(chr=paste("chr", UNIQBIN.DF$chr[log][ind], sep=""), 
                    start=bin.end[ind]-HiCres+1L, end=bin.end[ind], stringsAsFactors=FALSE)
  write.table(x=bed, file=paste0(out.dir, "/", suffix, "_", type), col.names=FALSE,
              row.names=FALSE, sep="\t", quote=FALSE)
  
  rm(subset, log, bin.end, ind); gc()
  
  print(paste0(type, " done!"), quote=FALSE)
  
}

# Motif analysis (de novo and known)
find_motifs_genome(
  x=bed,
  path = out.dir,
  # or add "r" for repeat-masked version e.g. "hg19r"
  genome = "hg19",
  # Vector of motif lengths to consider [default is c(8, 10, 12)]
  # Not advisable to look for motifs longer than 12 bp. Start with shorter
  # lengths and once you get promising motifs, try longer ones
  motif_length = 8:12,
  scan_size = "given",
  # Default optimize_count = 8
  # Specifies the number of motifs of each length to find
  optimize_count = 5,
  background = paste0(out.dir, "/", suffix, "_bgr"),
  local_background = FALSE,
  only_known = FALSE, 
  # Default fdr_num = 0; for de novo discovery,
  only_denovo = FALSE,
  fdr_num = 3,
  cores = nCPU, 
  cache = 100,
  overwrite = TRUE, 
  keep_minimal = FALSE
)

# Homer command build by marge
#cmd <- paste(
#  paste0(homer_base, "findMotifsGenome.pl"),
#  target_bed, genome, path,
#  "-len", paste0(motif_length, collapse = ","),
#  "-size", scan_size,
#  "-S", optimize_count,
#  "-p", cores,
#  "-cache", cache,
#  "-fdr", fdr_num
#)

## Default Parameters:
## Sequences not masked for repeats
## -mis Allowed mismatches set to 2
## -gc (use GC% for sequence content normalization, now the default)
## CG normalization
## Homer looks for motifs on both strands


## HOMER will still try to normalize the background to remove GC-bias and 
## will also perform autonormalization, normalization can be turned off
## Extract sequences from the genome corresponding to the regions in the 
## input file, filtering sequences that are >70% "N"
## Custom background regions provided are still GC-normalized and HOMER
## makes sure regions don't overlap with target regions

# Interpretation of results
## Enrichment p-values reported by HOMER should be very very significant
## (i.e. << 1e-50).  If this is not the case, there is a strong possibility
## that the experiment may have failed in one way or another. 
## In principle, in a motif is present in less than 5% of the targets 
## sequences, there may be a problem

# rm(list=ls())
 
