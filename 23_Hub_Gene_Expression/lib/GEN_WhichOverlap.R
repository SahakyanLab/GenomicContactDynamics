################################################################################
# This function takes the starting/ending coordinates and a stratifying factor 
# (space) for the query segments and the subject segments. It then returns a   
# matrix with all the indices of the overlapping segments.         

# For R < 3.5.0, the default maxgap differs across operating systems:
# Mac = -1L, Linux = 0L
# For R > 3.5.0, the default value is the same for Mac and Linux.
# Mac = Linux = -1L (one range has its start or end strictly inside the other 
# (i.e. non-disjoint ranges))
# Default minoverlap = 0L
# If type="any", one of minoverlap and maxgap should be set to default

# RangedData() is deprecated R>=3.6.0
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(IRanges)
# library(GenomicRanges)

WhichOverlap <- function(start.query=start, end.query=end, space.query=chr,
                         start.subject=anno$txStart, end.subject=anno$txEnd,
                         space.subject=anno$chrom, maxgap=-1L, minoverlap=1L,
                         type=c("any", "start", "end", "within", "equal")){  
                     
  # GRanges (used to be RangedData but was deprecated in 2019) function generates 
  # the corresponding objects. RangedData also reorders the data, hence we provide 
  # subject.ind and query.ind index object in order to account for the original 
  # indices for the data at the supplied order. This is kept even for the GRanges 
  # to make sure proper correspondence.
  subject <- GRanges( IRanges(start = start.subject, end = end.subject),
                      seqnames = space.subject, subject.ind = 1:length(start.subject) )
 
  query <- GRanges( IRanges(start = start.query, end = end.query),
                    seqnames = space.query, query.ind = 1:length(start.query))

  # maxgap, minoverlap and type identical to their counterpart in IRanges package
  ol <- as.matrix(findOverlaps(query=query, subject=subject, type=type,
                               maxgap=maxgap, minoverlap=minoverlap, select="all"))

  return( cbind(query=query[ol[,"queryHits"], ]$query.ind,
                subject=subject[ol[,"subjectHits"], ]$subject.ind) )

}
################################################################################
