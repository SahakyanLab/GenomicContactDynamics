################################################################################
# This function takes the starting/ending coordinates and a stratifying factor #
# (space) for the query segments and the subject segments. It then returns a   #
# matrix with all the indices of the overlapping segments.                     #
################################################################################
WhichOverlap <- function(start.query=start, end.query=end, space.query=chr,
                         start.subject=anno$txStart, end.subject=anno$txEnd,
                         space.subject=anno$chrom, maxgap=0L, minoverlap=1L){

  #> source("http://bioconductor.org/biocLite.R")
  #> biocLite("IRanges")
  library(IRanges)

  # RangedData function generates the corresponding objects, where it also reor-
  # ders the data, hence we provide subject.ind and query.ind index object in
  # order to account for the original indices for the data at the supplied order.
  subject <- RangedData( IRanges(start = start.subject, end = end.subject),
                         space = space.subject, subject.ind = 1:length(start.subject) )

  query <- RangedData( IRanges(start = start.query, end = end.query),
                       space = space.query,  query.ind = 1:length(start.query) )

  ol <- as.matrix(findOverlaps(query=query, subject=subject, type="any",
                            maxgap=maxgap, minoverlap=minoverlap, select="all"))

  return( cbind(query=query[ol[,"queryHits"], ]$query.ind,
                subject=subject[ol[,"subjectHits"], ]$subject.ind) )

}
################################################################################
