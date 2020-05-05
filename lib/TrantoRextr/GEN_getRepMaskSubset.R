## FUNCTION ####################################################################
# This function, given the file path to the repat masker database and the name #
# of the chromosome, reads the repeat data and returns the entries for only the#
# requested chromosome and requested column names. It directly accesses fread  #
# from data.table for fast reading.                                            #
################################################################################
getRepMaskSubset <- function(
 repmask.Filepath = "/Volumes/Data/Database/RepeatMasker_hg19/RepeatMasker_hg19",
 freadchar = "\t",
 chr = "chr1", # genoName column in repmask file. Returns all, if NULL.
 repmask.col = NULL # NULL returns all, otherwise, returns the columns of
                    # matching names, such as:
            # c("genoStart","genoEnd","strand","repName","repClass","repFamily")
){

  repmask <- data.table::fread(repmask.Filepath, sep = freadchar, header = TRUE,
                               stringsAsFactors = FALSE, data.table = FALSE)
  if(is.null(repmask.col[1])){
    col.ind <- 1:dim(repmask)[2]
  } else {
    col.ind <- match(repmask.col, names(repmask))
  }
  
  if(is.null(chr[1])){
    return(repmask[, col.ind])
  } else {
    return(repmask[repmask$genoName==chr, col.ind])
  }

}
################################################################################
