## FUNCTION ####################################################################
HicContactFeatures <- function(
  # Chromosome for which the data are to be pulled:
  chr  = "chr10",

  # Output path and filename suffix:
  outputPath = "/Users/alex/Desktop/Chrom_seq",
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
  LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR"
){
################################################################################
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
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings)
library(gtools)
# Needs an installation of data.table for fread.
library(IRanges)
################################################################################


# Loading the object as MELT.MX
load(paste0(combinedContactDir,"/human_",chr,"_allcontacts.RData"))

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

# Execution of the function writen above...
################################################################################
CHR = c("chr10","chr10","chr1","chr1")
MIN.DIST = c(2000000,500000,2000000,500000)
SUFFIX = c("min2Mb","min05Mb","min2Mb","min05Mb")

for(i in 3:4){ #length(CHR)
  
  print(paste0("**** ",CHR[i]," -- ",SUFFIX[i]," ****"))

  HicContactFeatures(
    # Chromosome for which the data are to be pulled:
    chr  = CHR[i],
    # Output path and filename suffix:
    outputPath = "/Users/alex/Desktop/Chrom_seq",
    suffix     = SUFFIX[i],
    # Combined contact directory path for Hi-C data and filtering configuration:
    combinedContactDir = 
      "/Volumes/Data/Database/GSE87112/combined_contacts/RAW_primary_cohort",
    hic.resol = 40000,
    min.dist  = MIN.DIST[i],
    top.perc  = 100,
    min.tiss  = 1,
    # RepeatMasker database path:
    repmaskPath = "/Volumes/Data/Database/RepeatMasker_hg19/RepeatMasker_hg19",
    # Human genome path and filename configuration:
    genomePath       = "/Volumes/Data/Database/human_genome_unmasked_37.73",
    genome.prefix    = "Homo_sapiens.GRCh37.73.dna.chromosome.",
    fastafile.ending = ".fa",
    # Path to the TrantoR library:
    LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR"
  )
  
}
################################################################################
