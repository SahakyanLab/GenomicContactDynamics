################################################################################
# This functiontakes two sequences as string objects (such as "AGTTGGGCCTGA"), #
# one being the reference sequence (pattern, for instance human.cdna), the other
# being the alignment subject sequence (for instance guineapig.cdna). It also  #
# needs the alignment type specification:                                      #
#                                                                              #
# type - type of alignment. One of "global", "local", "overlap", "global-local",
#        and "local-global" where "global" = align whole strings with end gap  #
#        penalties, "local" = align string fragments, "overlap" = align whole  #
#        strings without end gap penalties, "global-local" = align whole strings
#        in pattern with consecutive subsequence of subject, "local-global" =  #
#        align consecutive subsequence of pattern with whole strings in subject.
# the rest of the arguments are passed to pairwiseAlignmemnt function as ...   #
#                                                                              #
# It then performs am alignment of the subject on the reference sequence.      #
# align subject (for instance org.cdna) on the pattern (for instance human.cdna)
# pattern is the reference sequence, correctly indexing the alignment pattern  #
# relative to the reference sequence. The function returns an output with the  #
# following components:                                                        #
# $ref.seq  - supplied reference sequence (pattern)                            #
# $subj.seq - supplied subject sequence to align onto the reference            #
# $ref.afteralign  - reference sequence after the alignment (with possible indels)
#                    returns the same ref.seq, in case there is no alignment   #
# $subj.afteralign - subject sequence after the alignment (with possible indels)
#                    returns a string of <->s, in case there is no alignment   #
# $ref.seqind      - the vector of original ref.seq reference position indices #
#                    with NAs for indels.                                      #
# For this function to work, the R environment should have "Biostrings"        #
# installed:                                                                   #
# > source("http://bioconductor.org/biocLite.R")                               #
# > biocLite("Biostrings")                                                     #
################################################################################
SequenceAlignerTwo <- function(ref.seq ="AGGTGCGGGTACATGCTCCTTGAACTTGATCATTCGG",
                          subj.seq ="CCAGGTGCGGGTACATGCTGAACTTGAAAAAATCATTCGG",
                          type="global", # "local","overlap","global-local","local-global",
                          ...){

  library(Biostrings)

  align <- pairwiseAlignment(subject=subj.seq, pattern=ref.seq, type="global", ...)
  seq.i <- attr(attr(align, "pattern"),"range")[[1]]

  #####################
  if(length(seq.i)!=0){ # there is an outcome from this alignment

    ref.afteralign  <- as.character(attr(align, "pattern"))
    subj.afteralign <- as.character(attr(align, "subject"))

    indel <- attr(attr(align, "pattern"),"indel")[[1]]
    l.indel <- length(indel)
    # --------------- finding the positions that are <-> indels
    if(l.indel != 0){ # there are indels for reference sequence
      indel.ind <- NULL
      shift     <- 0
      for(h in 1:l.indel){
        indel.ind <- c(indel.ind,
                       (attr(indel[h],"start")+shift) :
                         (attr(indel[h],"start")+attr(indel[h],"width")-1+shift))
        shift <- shift + attr(indel[h],"width")
      }
      # correcting seq.i to ignore NAs and reset index counting
      seq.ind <- 1:nchar(ref.afteralign)
      seq.ind[indel.ind] <- NA
      seq.ind[-indel.ind] <- seq.i
      seq.i <- seq.ind
    }
    # ---------------------------------------------------------

  } else {
    # there is no outcome from this alignment
    ref.afteralign  <- ref.seq
    subj.afteralign <- paste(rep("-", nchar(ref.seq)), collapse="")
  }
  #####################

  RESULT                 <- NULL
  RESULT$ref.seq         <- ref.seq
  RESULT$subj.seq        <- subj.seq
  RESULT$ref.afteralign  <- ref.afteralign
  RESULT$subj.afteralign <- subj.afteralign
  RESULT$ref.seqind      <- seq.i

  return(RESULT)

}
################################################################################
