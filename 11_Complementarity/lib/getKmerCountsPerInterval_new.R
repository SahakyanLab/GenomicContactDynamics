## FUNCTION ####################################################################
# This function takes the chroosome ID, the vectors of starting and ending po- #
# sitions, k for k-mers and path/characteristics of the genome files. It re-   #
# turns a matrix in the same row order as the supplied startpos and endpos,    #
# with the column names being the requested k-mer permutation names, and the   #
# individual values representing the counts of the exact k-mers of a given (by #
# row) sequence. If a given sequence contains at least one "N" witin, all the  #
# values in that row will be populated with NA.                                #
#                                                                              #
# REQUIRES:                                                                    #
# For its getGenomicSeq() function:                                            #
# LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR"                              #
# library(foreach)                                                             #
# library(itertools)                                                           #
# library(doParallel)                                                          #
# source(paste0(LIB.TRANTOR, "/GEN_loadGenome.R"))                             #
# source(paste0(LIB.TRANTOR, "/GEN_readfasta.R"))                              #
# source(paste0(LIB.TRANTOR, "/UTIL_readLinesFast.R"))                         #
# The underlying functions may require fread() from <data.table> library.      #
# Next, the functions directly needed for this one:                            #
# source(paste0(LIB.TRANTOR, "/GEN_getGenomicSeq.R"))                          #
# source(paste0(LIB.TRANTOR, "/GEN_getKmers.R"))                               #
################################################################################
getKmerCountsPerInterval <- function(
chr           = "chr10",
startpos      = BINDATA.MX[,"startpos"], # vector of starting positions
endpos        = BINDATA.MX[,"endpos"],   # vector of ending positions
K             = 4,                       # k in the k-mer
PATH.genome   = "/Volumes/Data/Database/human_genome_unmasked_37.73",
genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
fastafile.ending = ".fa",
silent        = FALSE,
nCPU          = 2L
){
  
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  kmer.names <- names(getKmers(seq.string=DNAStringSet("AAAAAAAAAAAAAAA"),
                               k=K, method="Biostrings"))
  kmer.len <- length(kmer.names)
  
  toExport <- c("silent", "PATH.genome", "genome.prefix", "fastafile.ending",
                "chr", "silent", "startpos", "endpos", "K", "kmer.len")
  
  bin.len <- length(startpos)
  
  #### PARALLEL EXECUTION #########
  
  KMER.MX <- foreach(itr=isplitVector(1:bin.len, chunks=nCPU), 
                     .inorder=TRUE, .combine="rbind",
                     .export=toExport, .noexport=ls()[!ls()%in%toExport]
          
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      seq <- getGenomicSeq(PATH.genome=PATH.genome,
                           genome.prefix=genome.prefix,
                           fastafile.ending=fastafile.ending,
                           chr.id=strsplit(chr,"chr")[[1]][2],
                           remove.other.loads = TRUE,
                           silent=FALSE, split = FALSE,
                           borders=c(startpos[i], endpos[i]) )
      
      if( length(grep("N", seq))==1 ){
        
        # If there is N in the sequence populate the counts with NA
        kmers <- rep(NA, times=kmer.len) 
        
      } else {
        
        # Kmers in order already
        kmers <- getKmers(seq.string=DNAStringSet(seq), k=K, method="Biostrings")
        
      }

      # Number of unmasked character in bin sequence and kmer counts 
      return( c(numUMChar=nchar(seq), kmers) )
        
    })
    
    return( do.call("rbind", chunk) )
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  dimnames(KMER.MX)[[2]] <- c("numUMChar", kmer.names)
  
  return(KMER.MX)
    
}
################################################################################
