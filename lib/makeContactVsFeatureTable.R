################################################################################
# Assign features (from FEATURE.BIN.MX) to contacts (from PERSIST.MX) with filtering
################################################################################
################################################################################
makeContactVsFeatureTable <- function(setCountPerBin = c(1,1),
                                      featureBinObj = FEATURE.BIN.MX,
                                      persistObj = PERSIST.MX,
                                      featureToAdd =  featureToAdd){
  
  # Filter for bins that satisfy setCountPerBin
  if( !is.null(setCountPerBin) ){
    
    if( length(setCountPerBin)==1 ){
      incl.bin <- as.numeric( featureBinObj[featureBinObj[,"countPerBin"]>=setCountPerBin,"bin"] )
    } else if( length(setCountPerBin)==2 ){
      test1 <- as.numeric(featureBinObj[,"countPerBin"])>=setCountPerBin[1] 
      test2 <- as.numeric(featureBinObj[,"countPerBin"])<=setCountPerBin[2]
      incl.bin <- as.numeric( featureBinObj[test1 & test2,"bin"] )
    } else {
      stop("If all contacts to be included, set setCountPerBin to 1.")
    }
    
  }
  
  # Initialize final table
  mx <- matrix( ncol=(length(featureToAdd)+4), 
                nrow=(length(persistObj$ntis)),
                dimnames=list(NULL, 
                              c("include", "i", "j", "ntis", featureToAdd)) )
  mx[,c("i", "j", "ntis")] <- cbind(i=as.numeric(persistObj$hits[,"i"]),
                                    j=as.numeric(persistObj$hits[,"j"]),
                                    ntis=as.numeric(persistObj$ntis))
  test3 <- as.numeric(persistObj$hits[,"i"]%in%incl.bin)
  test4 <- as.numeric(persistObj$hits[,"j"]%in%incl.bin)
  # 1 = yes included; NA = not included
  mx[test3 & test4, "include" ] <- 1
  
  return(mx)
}
