################################################################################
# Make sure file has chr-start-end and take only proper chromosomes
# REQUIRES:
# remove unnecessary lines in the file
################################################################################
turnToBed <- function(file="/Users/ltamon/Database/maskingLib/hg19_CpG",
                          chr=paste("chr", c(1:22, "X", "Y", "MT"), sep=""),
                          chrCol=1, startCol=2, endCol=3, sep="\t", header=FALSE){
  library(data.table)
  tbl <- fread(file=file, header=header, data.table=FALSE)
  setnames(x=tbl, old=c(chrCol, startCol, endCol), new=c("chr", "start", "end")) 
  return( tbl[tbl$chr%in%chr, ] )
}

