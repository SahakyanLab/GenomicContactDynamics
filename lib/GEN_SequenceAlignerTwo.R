################################################################################
# This function takes two sequences as string objects (such as "AGTTGGGCCTGA"), #
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
# Returns only the pairwise alignment score
# For this function to work, the R environment should have "Biostrings"        #
# installed:                                                                   #
# > source("http://bioconductor.org/biocLite.R")                               #
# > biocLite("Biostrings") 
# > library(Biostrings)
################################################################################
SequenceAlignerTwo <- function(ref.seq = "AGGTGCGGGTACATGCTCCTTGAACTTGATCATTCGG",
                               # RC:     "CCGAATGATCAAGTTCAAGGAGCATGTACCCGCACCT"
                               subj.seq = "CCAGGTGCGGGTACATGCTGAACTTGAAAAAATCATTCGG",
                               # RC:       CCGAATGATTTTTTCAAGTTCAGCATGTACCCGCACCTGG
                               type = "local", # "local","overlap","global-local","local-global",
                               submat = nucleotideSubstitutionMatrix(match=5, mismatch=-4, 
                                                                     baseOnly=TRUE, type="DNA"), 
                               # c(gapOpening, gapExtension)
                               gapCosts = c(10,6), 
                               ...){
                               
  align <- pairwiseAlignment(subject=subj.seq, pattern=ref.seq, type=type,
                             substitutionMatrix=submat,
                             gapOpening=gapCosts[1], gapExtension=gapCosts[2])
  
  return(attr(align, "score"))
  
}

################################################################################
