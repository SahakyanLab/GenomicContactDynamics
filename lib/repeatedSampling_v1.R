################################################################################
# Function specific for randomiseContact.R
################################################################################
################################################################################
repeatedSampling <- function(pool = bins,
                             numContOrig = 173,
                             gapmin = gapmin,
                             ijorig = ij.sub){
  
  # Get random pairs
  pair.lst <- getRandomPairs(pool=pool, replacement=FALSE)
  ij.new <- do.call("rbind", pair.lst)
  rm(pair.lst); gc()
  
  # Duplicated contacts
  test.dup <- duplicated(ij.new) 
  # Contacts not satisfying gap requirement
  gap.new <- ij.new[,2] - ij.new[,1] 
  test.gap <- gap.new < gapmin
  ij.new <- ij.new[!test.dup | !test.gap, ]
  colnames(ij.new) <- c("i", "j")
  # Number of randomised contacts (filtered for duplicates and gap req)
  numcont.new <- nrow(ij.new)
  
  ij.new.str <- paste(ij.new[,"i"], ij.new[,"j"], sep=",")
  # Contact overlaps between original and randomised sets
  olap <- intersect( paste(ijorig[,"i"], ijorig[,"j"], sep=","),
                     ij.new.str ) 
  info <- c( 100*(numcont.new/numContOrig),
             100*(length(olap)/numcont.new) )
  rm(ij.new.str, olap, numcont.new);gc()
  return( list(ij.new, info) )
  
}