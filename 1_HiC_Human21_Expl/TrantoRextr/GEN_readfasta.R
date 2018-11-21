###############################################################################
## This function takes the name/path of the fasta file and parses it, retur- ##
## ning the parsed.fasta composite object, that has two components:          ##
## $header - for the fasta header stored as a single string (NULL if none),  ##
## $seq    - for the vector of all the sequence letters (NULL if none),      ##
## $length - a numeric of the length of the sequence (or 0 if $seq is NULL). ##
## By default, it will always return upper case characters in the read seq., ##
## however, this can be controlled by the case="UPPER"/"LOWER" argument.     ##
## If the file is multifasta, it will returns all the headers in the $header ##
## section, but will merge the sequences without the separator identification##
## in the $seq section, with the joint length in the $length section.        ##
## If the argument <fasta> is not NULL, the function assumes that an already ##
## read in fasta file is passed (prior read with readLines) and treats the   ##
## passed object as such.                                                    ##
## If split = FALSE, the sequence will be read as one long string.           ##
###############################################################################
readfasta <- function(filename="filename.fasta", case="UPPER", fasta=NULL,
                       split=TRUE, fastread=TRUE){

  if(is.null(fasta)){
    if(fastread==TRUE){

      if(exists("readLinesFast", mode="function")){
        fasta <- readLinesFast(filename, freadchar="@")
      } else {
          stop("readfasta: Fast reading is requested. Make readLinesFast() available first!")
      }

    } else {
      fasta <- readLines(filename)
    }
  }

  parsed.fasta <- NULL

  discard.ind <- grep(">", fasta, fixed=TRUE)
  if(length(discard.ind)==0){
    parsed.fasta$header <- NULL
  } else {
    parsed.fasta$header <- fasta[discard.ind]
  }

  discard.ind <- c(discard.ind, which(fasta==""|fasta==" "|fasta=="   "))
  if(length(discard.ind)!=0){
    fasta <- fasta[-discard.ind]
  }

  parsed.fasta$seq <- NULL
  if(length(fasta)!=0){
    if(split==TRUE){
      parsed.fasta$seq <- unlist(strsplit(paste(fasta,collapse=""),""))
      parsed.fasta$length <- length(parsed.fasta$seq)
    } else {
      parsed.fasta$seq    <- paste(fasta,collapse="")
      parsed.fasta$length <- nchar(parsed.fasta$seq)
    }

    if(case=="UPPER"){
      parsed.fasta$seq <- toupper(parsed.fasta$seq)
    } else if(case=="LOWER"){
      parsed.fasta$seq <- tolower(parsed.fasta$seq)
    }

  }

  return(parsed.fasta)

}
###############################################################################
