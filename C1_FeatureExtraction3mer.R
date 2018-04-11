## FUNCTION ####################################################################
FeatureExtraction3mer<- function(
  # Chromosome for which the data are to be pulled:
  chr  = "chr1",

  # Database path and filename suffix:
  featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT",
  # "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc"
  suffix        = "min2Mb",

  # Human genome path and filename configuration:
  genomePath       = "/home/alex/Desktop/CHROMSEQ/human_genome_unmasked_37.73",
  # "/Volumes/Data/Database/human_genome_unmasked_37.73"
  genome.prefix    = "Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending = ".fa",

  # Path to the TrantoR library, used only if BINKMER3 is not yet calculated:
  LIB.TRANTOR = "/home/alex/Desktop/CHROMSEQ/TrantoR",
  # "/Users/alex/GIT/GITrepo/TrantoR"

  # Number of CPU cores to be used
  nCPU = 5, # optimal for the task on bucephalus (takes ~5GB per core for chr1)
  scale = 40000 # the length of DNA/Hi-C resolution for Gfree scaling to per-nt
){
################################################################################
library(doMC)
library(foreach)
library(itertools)
registerDoMC(cores=nCPU)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings)
################################################################################

################################
################################
generate.BINKMER3 <- function(){
  BINREP.filepath   <- paste0(featureDBPath,"/",chr,"_BinRep_",suffix,".RData")
  BINKMER3.filepath <- paste0(featureDBPath,"/",chr,"_BinKmer3_",suffix,".RData")

  if(!file.exists(BINKMER3.filepath)){
    #----------------------------------------------
    source(paste0(LIB.TRANTOR, "/GEN_readfasta.R"))
    source(paste0(LIB.TRANTOR, "/UTIL_readLinesFast.R"))
    source(paste0(LIB.TRANTOR, "/GEN_loadGenome.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getGenomicSeq.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getKmers.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getKmerCountsPerInterval.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getKmerHybridisationGs.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getMELTMXpersist.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getRepMaskSubset.R"))
    source(paste0(LIB.TRANTOR, "/GEN_WhichOverlap.R"))
    source(paste0(LIB.TRANTOR, "/GEN_getBinsOverlapRepMasker.R"))
    library(gtools)
    # Needs an installation of data.table for fread.
    library(IRanges)
    #----------------------------------------------
    kmerGfree.par <- getKmerHybridisationGs(k=3, plot=TRUE)
     write.table(kmerGfree.par, file=paste0("Gfree_",3,"mer.par"), row.names=FALSE)

    #----------------------------------------------
    if(file.exists(BINREP.filepath)){
      load(BINREP.filepath)
    } else {
      BINREP.MX <- getBinsOverlapRepMasker(chr=chr, hic.resol=hic.resol,
                                           bins=unique.bins,
                                           repmask.Filepath=repmaskPath)
      save(BINREP.MX, file=BINREP.filepath)
    }
    bins     <- BINREP.MX[,"bins"]
    startpos <- BINREP.MX[,"startpos"]
    endpos   <- BINREP.MX[,"endpos"]
    rm(BINREP.MX); gc()
    #----------------------------------------------

    BINKMER3.MX <- getKmerCountsPerInterval(
      chr=chr,
      startpos=startpos,
      endpos=endpos,
      K=3,
      PATH.genome=genomePath,
      genome.prefix=genome.prefix,
      fastafile.ending=fastafile.ending,
      silent=TRUE
    )
    BINKMER3.MX <- cbind(bins=bins, startpos=startpos, endpos=endpos, BINKMER3.MX)
    save(BINKMER3.MX,  file=BINKMER3.filepath)
  }
  return(NULL)
}
################################
################################

################################
generate.BINKMER3()
rm(generate.BINKMER3); gc()
load(paste0(featureDBPath,"/",chr,"_BinKmer3_",suffix,".RData"))
################################

################################
k     <- 3
id    <- paste0(chr,"_",suffix,"_kmer",k)
BINKMER.name <- load(paste0(featureDBPath,"/",chr,"_BinKmer",k,"_",suffix,".RData"))
eval(parse(
 text=paste0("BINKMER.MX <- ",BINKMER.name,"; rm(",BINKMER.name,"); rm(BINKMER.name); gc()")
))
load(paste0(featureDBPath,"/",chr,"_Persist_",suffix,".RData"))
gfree <- read.table(file=paste0("./Gfree_",k,"mer.par"), header=TRUE, stringsAsFactors=FALSE)

############################
kmer.vec <- dimnames(BINKMER.MX)[[2]][-c(1:3)]
kmer.vec.revcomp <- as.vector(reverseComplement(DNAStringSet(kmer.vec)))
kmer.revcomp.ind <- match(kmer.vec, kmer.vec.revcomp)
unique.kmer.ind <- unique(pmin(1:64, kmer.revcomp.ind))
rm(kmer.vec, kmer.vec.revcomp)
############################
ij.mx    <- cbind(i=PERSIST.MX$hits[,"i"], j=PERSIST.MX$hits[,"j"])
ntis.vec <- PERSIST.MX$ntis
foreach.export <- c("BINKMER.MX","ij.mx","kmer.revcomp.ind","gfree","scale")
rm(PERSIST.MX)
############################
gc()
################################

############################
# Eliminating ij pairs where at least one sequence segment has N or Ns.
binsNA <- BINKMER.MX[is.na(BINKMER.MX[,4]), "bins"]
for(b in binsNA){
  banned.rows <- which( ij.mx[,"i"] == b | ij.mx[,"j"] == b)
  if(length(banned.rows)!=0){
    ij.mx <- ij.mx[-banned.rows,]
    ntis.vec <- ntis.vec[-banned.rows]
  }
}; rm(b, binsNA)
############################
hits.len <- length(ij.mx[,1])
dimnames(BINKMER.MX)[[2]] <- NULL

gc()

#### FOREACH EXECUTION #########
PAIR.MX <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                   .combine="rbind", .inorder=TRUE,
                   .export=foreach.export,
                   .noexport=ls()[!ls()%in%foreach.export]
                  ) %dopar% {

  pair.data <- sapply(itr,
   FUN=function(i){
     kmer.mx <- BINKMER.MX[which(BINKMER.MX[,1] %in% ij.mx[i,]), -c(1:3) ]

     kmerAall <- (kmer.mx[1,] + kmer.mx[1,kmer.revcomp.ind]) #(kmerA + kmerArc)
     kmerBall <- (kmer.mx[2,] + kmer.mx[2,kmer.revcomp.ind]) #(kmerB + kmerBrc)

     #countsA <- kmerAall[unique.kmer.ind]
     #names(countsA) <- paste("a",names(countsA), sep="")
     #countsB <- kmerBall[unique.kmer.ind]
     #names(countsB) <- paste("b",names(countsB), sep="")

     return(c( kmerAall[unique.kmer.ind],
               kmerBall[unique.kmer.ind],
               sum(pmin(kmerAall, kmerBall)*gfree[,2])/(2*scale), # e=
               sum(abs(kmerAall - kmerBall)) # sumabsdif=
           ))
   }, simplify=TRUE) # USE.NAMES=FALSE

   return(t(pair.data))
}
### END OF FOREACH EXECUTION ###

dimnames(PAIR.MX)[[2]] <-
   c("aAAA","aAAC","aAAG","aAAT","aACA","aACC","aACG","aACT","aAGA","aAGC","aAGG",
     "aATA","aATC","aATG","aCAA","aCAC","aCAG","aCCA","aCCC","aCCG","aCGA","aCGC",
     "aCTA","aCTC","aGAA","aGAC","aGCA","aGCC","aGGA","aGTA","aTAA","aTCA",
     "bAAA","bAAC","bAAG","bAAT","bACA","bACC","bACG","bACT","bAGA","bAGC","bAGG",
     "bATA","bATC","bATG","bCAA","bCAC","bCAG","bCCA","bCCC","bCCG","bCGA","bCGC",
     "bCTA","bCTC","bGAA","bGAC","bGCA","bGCC","bGGA","bGTA","bTAA","bTCA",
     "e","sumabsdif")

#save(PAIR.MX, file="PAIR_MX_TEST.RData")
#load("PAIR_MX_TEST.RData")

PAIR.MX <- cbind(ij.mx, ntis=ntis.vec, PAIR.MX)
save(PAIR.MX, file=paste0(featureDBPath,"/",chr,"_PairMx_",suffix,".RData"))

}
################################################################################

# Execution of the above function:
FeatureExtraction3mer(
  chr = "chr1",
  featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT", # "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc"
  suffix = "min2Mb",
  genomePath = "/home/alex/Desktop/CHROMSEQ/human_genome_unmasked_37.73", # "/Volumes/Data/Database/human_genome_unmasked_37.73"
  genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending = ".fa",
  LIB.TRANTOR = "/home/alex/Desktop/CHROMSEQ/TrantoR", # "/Users/alex/GIT/GITrepo/TrantoR"
  nCPU = 5, # optimal for the task on bucephalus (takes ~5GB per core for chr1)
  scale = 40000 # the length of DNA/Hi-C resolution for Gfree scaling to per-nt
)
