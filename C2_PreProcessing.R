## FUNCTION ####################################################################
InitialProcessing3mer <- function(
  featureDBPath = "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # "/home/alex/Desktop/CHROMSEQ/OUT"
  suffix        = "min2Mb",
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Path to the TrantoR library:
  LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR", #"/home/alex/Desktop/CHROMSEQ/TrantoR"
  difFeatures = TRUE # whether difference features should be created, rather
                     # than absolute features for each interacting DNA spans.
){
################################################################################
  source(paste0(LIB.TRANTOR, "/ML_SanityCheckAndPartition.R"))

  loadPairMx <- function(chr=chr, featureDBPath=featureDBPath, suffix=suffix){
    load(paste0(featureDBPath,"/",chr,"_PairMx_",suffix,".RData"))
    return(PAIR.MX)
  }

  data   <- loadPairMx(chr=chr, suffix=suffix, featureDBPath=featureDBPath)
  data   <- as.data.frame(data, stringsAsFactors=FALSE)
  data$j <- paste(as.integer(data$i), as.integer(data$j), sep="_")
  data   <- data[,-1]
  names(data)[1:2] <- c("ID","Value")

  if(difFeatures==FALSE){
    SanityCheckAndPartition(data = data, partition = 40.0, seed = 185,
                            naValueFilter = TRUE, sdZeroFilter = TRUE,
                            corHighFilter = 0.98, outDir = "./C")
  } else {
    nms <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG",
             "ATA","ATC","ATG","CAA","CAC","CAG","CCA","CCC","CCG","CGA","CGC",
             "CTA","CTC","GAA","GAC","GCA","GCC","GGA","GTA","TAA","TCA")

    data <- cbind(data[,1:2],
                  abs(data[,paste0("a", nms)] - data[,paste0("b", nms)]),
                  data[,67:68])
    names(data)[3:34] <- paste0("dif", nms)

    SanityCheckAndPartition(data = data, partition = 40.0, seed = 185,
                            naValueFilter = TRUE, sdZeroFilter = TRUE,
                            corHighFilter = 0.98, outDir = "./C")
  }

}
################################################################################

InitialProcessing3mer(
  featureDBPath = "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # "/home/alex/Desktop/CHROMSEQ/OUT"
  suffix        = "min2Mb",
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Path to the TrantoR library:
  LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR", #"/home/alex/Desktop/CHROMSEQ/TrantoR"
  difFeatures = TRUE # whether difference features should be created, rather
                     # than absolute features for each interacting DNA spans.
)

# dimnames(PAIR.MX)[[2]] <-
#   c("aAAA","aAAC","aAAG","aAAT","aACA","aACC","aACG","aACT","aAGA","aAGC","aAGG",
#     "aATA","aATC","aATG","aCAA","aCAC","aCAG","aCCA","aCCC","aCCG","aCGA","aCGC",
#     "aCTA","aCTC","aGAA","aGAC","aGCA","aGCC","aGGA","aGTA","aTAA","aTCA",
#     "bAAA","bAAC","bAAG","bAAT","bACA","bACC","bACG","bACT","bAGA","bAGC","bAGG",
#     "bATA","bATC","bATG","bCAA","bCAC","bCAG","bCCA","bCCC","bCCG","bCGA","bCGC",
#     "bCTA","bCTC","bGAA","bGAC","bGCA","bGCC","bGGA","bGTA","bTAA","bTCA",
#     "e","sumabsdif")
