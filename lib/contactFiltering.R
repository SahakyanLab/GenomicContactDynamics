################################################################################
# Function specific for randomiseContact.R
################################################################################
################################################################################
contactFiltering <- function(numContOrig = numcont.orig,
                             gapmin = gapmin,
                             ijorig = ij.sub,
                             ijnew = ij.new){
  
  test.dup <- duplicated(cbind(ij.new$i, ij.new$i))
  test.belowGap <- ij.new$j-ij.new$i < gapmin 
  test.rev <- ij.new$j < ij.new$i
  
  
  ijnew <- ijnew[(!test.dup & !test.gap), ]
  colnames(ijnew) <- c("i", "j")
  # Number of randomised contacts (filtered for duplicates and gap req)
  numcont.new <- nrow(ijnew)
  
  ijnew.str <- paste(ijnew[,"i"], ijnew[,"j"], sep=",")
  # Contact overlaps between original and randomised sets
  olap <- intersect( paste(ijorig[,"i"], ijorig[,"j"], sep=","),
                     ijnew.str ) 
  info <- c( 100*(numcont.new/numContOrig),
             100*(length(olap)/numcont.new) )
  rm(ijnew.str, olap, numcont.new);gc()
  return( list(ijnew, info) )
  
}