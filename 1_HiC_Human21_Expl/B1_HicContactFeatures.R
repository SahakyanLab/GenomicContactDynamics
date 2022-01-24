################################################################################
# Count k-mers for all bins of a chr
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/11_Complementarity")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/11_Complementarity")
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/11_Complementarity")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
data.dir = paste0(home.dir, "/Database")

LIB.TRANTOR = paste0(home.dir, "/DPhil/lib/TrantoRextr")
combinedContactDir = paste0(data.dir, "/rice_HiC/Liu2017_50kb_combined_contacts/ICE_Liu")
genomePath = paste0(data.dir, "/rice_genome_unmasked_IRGSP1.0") 
outputPath = paste0(data.dir, "/HiC_features_Liu2017oryza_HiC_NORMpc/features_50kb")  
repmaskPath = paste0(data.dir, "/ucsc_tables/hsa_RepeatMasker/RepeatMasker_hg19")
### OTHER SETTINGS #############################################################
CHR      = paste0("chr", 1:12)  #rep( paste0("chr",c(1:22,"X")), 2 )
MIN.DIST = rep(0, times=length(CHR))  #c( rep(500000, 23), rep(2000000, 23) )
SUFFIX   = rep("min0Mb", times=length(CHR))  #c( rep("min05Mb"outputPath         = "./OUT"
genome.prefix = "Oryza_sativa_Nipponbare.IRGSP1.0.dna.chromosome." 
fastafile.ending = ".fa"
hic.resol = 50000 #2000 #40000
top.perc  = 100
min.tiss  = 1
species.id = "osa"
## FUNCTION ####################################################################
HicContactFeatures <- function(
  # Chromosome for which the data are to be pulled:
  chr  = "chr10",

  # Output path and filename suffix:
  outputPath = "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  suffix     = "min2Mb",

  # Combined contact directory path for Hi-C data and filtering configuration:
  combinedContactDir =
         "/Volumes/Data/Database/GSE87112/combined_contacts/RAW_primary_cohort",
  hic.resol = 40000,
  min.dist  = 2000000,
  top.perc  = 100,
  min.tiss  = 1,

  # RepeatMasker database path:
  repmaskPath = "/Volumes/Data/Database/RepeatMasker_hg19/RepeatMasker_hg19",

  # Human genome path and filename configuration:
  genomePath       = "/Volumes/Data/Database/human_genome_unmasked_37.73",
  genome.prefix    = "Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending = ".fa",

  # Path to the TrantoR library:
  LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR",
  
  species.id = "human"
){
################################################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings)
library(gtools)
# Needs an installation of data.table for fread.
library(data.table)
library(IRanges)
library(GenomicRanges)
source(paste0(LIB.TRANTOR, "/GEN_readfasta.R"))
source(paste0(LIB.TRANTOR, "/UTIL_readLinesFast.R"))
source(paste0(LIB.TRANTOR, "/GEN_loadGenome.R"))                #***
source(paste0(LIB.TRANTOR, "/GEN_getGenomicSeq.R"))             #***
source(paste0(LIB.TRANTOR, "/GEN_getKmers.R"))                  #***
source(paste0(LIB.TRANTOR, "/GEN_getKmerHybridisationGs.R"))    #***
source(paste0(LIB.TRANTOR, "/GEN_getMELTMXpersist.R"))          #***
source(paste0(LIB.TRANTOR, "/GEN_getRepMaskSubset.R"))          #***
source(paste0(LIB.TRANTOR, "/GEN_WhichOverlap.R"))
source(paste0(LIB.TRANTOR, "/GEN_getBinsOverlapRepMasker.R"))   #***
source(paste0(LIB.TRANTOR, "/GEN_getKmerCountsPerInterval.R"))  #***
################################################################################
# Loading the object as MELT.MX
load(paste0(combinedContactDir,"/",species.id,"_",chr,"_allcontacts.RData"))

#***#

PERSIST.MX <- getMELTMXpersist(MELT.MX=MELT.MX, min.dist=min.dist,
                               hic.resol=hic.resol, top.perc=top.perc,
                               min.tiss=min.tiss)
# > names(PERSIST.MX) > [1] "hits"    "ntis"    "control"
rm(MELT.MX); gc()

unique.bins <- unique(as.vector(as.matrix(PERSIST.MX$hits[,c("i","j")])))
save(PERSIST.MX,  file=paste0(outputPath,"/",chr,"_Persist_",suffix,".RData"))
rm(PERSIST.MX); gc()

#***#

BINREP.MX <- getBinsOverlapRepMasker(
  chr=chr,
  hic.resol=hic.resol,
  bins=unique.bins,
  repmask.Filepath=repmaskPath
)
bins     <- BINREP.MX[,"bins"]
startpos <- BINREP.MX[,"startpos"]
endpos   <- BINREP.MX[,"endpos"]

save(BINREP.MX,  file=paste0(outputPath,"/",chr,"_BinRep_",suffix,".RData"))
rm(BINREP.MX); gc()

#***#

BINKMER4.MX <- getKmerCountsPerInterval(
  chr=chr,
  startpos=startpos,
  endpos=endpos,
  K=4,
  PATH.genome=genomePath,
  genome.prefix=genome.prefix,
  fastafile.ending=fastafile.ending,
  silent=TRUE
)
BINKMER4.MX <- cbind(bins=bins, startpos=startpos, endpos=endpos, BINKMER4.MX)
save(BINKMER4.MX,  file=paste0(outputPath,"/",chr,"_BinKmer4_",suffix,".RData"))
rm(BINKMER4.MX); gc()

#***#

BINKMER7.MX <- getKmerCountsPerInterval(
  chr=chr,
  startpos=startpos,
  endpos=endpos,
  K=7,
  PATH.genome=genomePath,
  genome.prefix=genome.prefix,
  fastafile.ending=fastafile.ending,
  silent=TRUE
)
BINKMER7.MX <- cbind(bins=bins, startpos=startpos, endpos=endpos, BINKMER7.MX)
save(BINKMER7.MX,  file=paste0(outputPath,"/",chr,"_BinKmer7_",suffix,".RData"))
rm(BINKMER7.MX); gc()

#***#

kmerGfree.par <- getKmerHybridisationGs(k = 4, plot = TRUE)
  write.table(kmerGfree.par, file=paste0("Gfree_",4,"mer.par"), row.names=FALSE)

kmerGfree.par <- getKmerHybridisationGs(k = 7, plot = TRUE)
  write.table(kmerGfree.par, file=paste0("Gfree_",7,"mer.par"), row.names=FALSE)

print("HicContactFeatures is DONE!", quote=FALSE)

}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(i in 1:length(CHR)){

  print(paste0("**** ",CHR[i]," -- ",SUFFIX[i]," ****"))

  HicContactFeatures(
    # Chromosome for which the data are to be pulled:
    chr  = CHR[i],
    # Output path and filename suffix:
    outputPath = outputPath,
    suffix     = SUFFIX[i],
    # Combined contact directory path for Hi-C data and filtering configuration:
    combinedContactDir = combinedContactDir,
    hic.resol = hic.resol,
    min.dist  = MIN.DIST[i],
    top.perc  = top.perc,
    min.tiss  = min.tiss,
    # RepeatMasker database path:
    repmaskPath = repmaskPath,
    # Human genome path and filename configuration:
    genomePath       = genomePath,
    genome.prefix    = genome.prefix,
    fastafile.ending = fastafile.ending,
    # Path to the TrantoR library:
    LIB.TRANTOR = LIB.TRANTOR,
    
    species.id = species.id
    
  )

}
################################################################################

# rm(list=ls()); gc()