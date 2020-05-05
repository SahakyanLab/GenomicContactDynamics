################################################################################
# "table" is a single processor method that partitions the genome at once and uses
# the built-in table function. The method is simple but efficient for only the
# short sequences.
#
# "index" method uses the complex vector updating method that utilises the
# Seq2Index function of TrantoR, that is written based on the "Bioinformatics
# Algorithm" by Pavzner and Companeue.
#
# "Biostrings" method uses the ready functions made available through the
# Biostrings library. Is the fastest method. The outputs are compatible, but the
# input sequences are of different format. The "Biostrings" method requires
# XStringSet object as a sequence input, which can be created while, for instance
# reading the input fasta file via the following Biostrings function:
# > readDNAStringSet(file=test.fasta, format="fasta", seek.first.rec=TRUE)
#
# DEPENDENCY:
# The method="index" requires GEN_Seq2Index.R and "gtools" libraries.
# The method="Biostrings" requires the Biostrings library from bioconductor.
################################################################################
getKmers <- function(seq.string = seq,    # Sequence as a vector of letters.
                                          #  if method="Biostrings", should be a
                                          #  XStringSet object.
                     k          = 7,      # Length of the k-mers.
                     method     = "index" # The method to be used.
                                          # can also be "table" and "Biostrings".
                    ){

  #------------------------
  if(method == "Biostrings"){
    library("Biostrings")
    #seq.string <- readDNAStringSet(file=test.fasta, format="fasta", seek.first.rec=TRUE)
    count.vec <- oligonucleotideFrequency(x=seq.string, width=k, step=1,
                                          as.array=FALSE, as.prob=FALSE,
                                          fast.moving.side="right", with.labels=TRUE,
                                          simplify.as="collapsed")
    return(count.vec)
  } else {

    end   <- length(seq.string)-k+1

    if(method == "table"){
      kmers <- substring( paste(seq.string, collapse=""), first=1:end, last=(1:end)+k-1 )
      return(table(kmers))
    }

    #*********************
    if(method == "index"){
      library(gtools)
      #-------------------------------------------------------------------------------
      print("Generating the permutations of the query k-mer...", quote=F)
      query.permuts <- permutations(n=4, r=k, v=c("A","G","T","C"), repeats.allowed=TRUE)
      #-------------------------------------------------------------------------------
      print("Optimising the object that holds the query segment permutations...", quote=F)
      query.permuts <- sapply(1:length(query.permuts[,1]), FUN=function(i){
                                return(paste(query.permuts[i,],collapse=""))
                              }, simplify=TRUE, USE.NAMES=FALSE)
      count.vec <- rep(0, length(query.permuts))
      names(count.vec) <- query.permuts

      #first <- 1:end
      last  <- (1:end)+k-1

      print("Scanning the sequence...", quote=F)
      pb <- txtProgressBar(min=1, max=end, style=3)
      for(i in 1:end){
        ##print(paste(i," of ", end, sep=""))
        query.seq <- seq.string[i:last[i]]
        # Checking the the letters in the pattern are only either A, or C, or G, or T,
        # hence skipping the ones that contain N and other possible letters.
        if(sum(is.na(match(unique(query.seq),c("A","C","G","T"))))==0){
          ##print(paste(query.seq, collapse=""))
          indx      <- Seq2Index(query.seq=query.seq) + 1
          count.vec[indx] <- count.vec[indx] + 1
          setTxtProgressBar(pb, i)
        }
      }
      close(pb)

      return(count.vec)
    }
    #*********************
  }
  #------------------------
}
################################################################################

