################################################################################
# Function specific for randomiseContact.R
################################################################################
################################################################################
contactFiltering <- function(numContOrig = numcont.orig,
                             gapmin = gapmin,
                             ijorig = ij.sub,
                             ijnew = ij.new){
  
  # Duplicated contacts
  test.dup <- duplicated(ijnew) 
  # Contacts not satisfying gap requirement
  gap.new <- (ijnew[,2]*40000) - (ijnew[,1]*40000)
  test.gap <- gap.new < gapmin
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