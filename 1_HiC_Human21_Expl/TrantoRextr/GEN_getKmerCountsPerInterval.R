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
silent        = FALSE
){
    
    kmer.names <- names(getKmers(seq.string=DNAStringSet("AAAAAAAAAAAAAAA"),
    k=K, method="Biostrings"))
    ncol       <- length(kmer.names)
    nrow       <- length(startpos)
    
    KMER.MX    <- matrix(0, ncol=ncol, nrow=nrow)
    dimnames(KMER.MX)[[2]] <- kmer.names
    
    for(rw in 1:nrow){
        if(!silent){ print(rw) }
        seq <- getGenomicSeq(PATH.genome = PATH.genome,
        genome.prefix = genome.prefix,
        fastafile.ending = fastafile.ending,
        chr.id = strsplit(chr,"chr")[[1]][2],
        remove.other.loads = TRUE,
        silent = silent, split = FALSE,
        borders = c(startpos[rw], endpos[rw]) )
        if(length(grep("N", seq))==1){
            KMER.MX[rw,] <- NA # if there is N in the sequence, populate the counts with NA
        } else {
            kmers <- getKmers(seq.string=DNAStringSet(seq), k=K, method="Biostrings")
            KMER.MX[rw, ] <- kmers # already in order, hence no need to
            # write: KMER.MX[rw,names(kmers)] <- kmers
        }
    }
    
    return(KMER.MX)
    
}
################################################################################
