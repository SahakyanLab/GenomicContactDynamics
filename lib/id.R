################################################################################
#Assign features (from FEATURE.BIN.MX) to contacts (from PERSIST.MX) 
################################################################################
################################################################################
function(setAgeCountPerBin = 1,
         featBinObj = FEATURE.BIN.MX,
         persistObj = PERSIST.MX){
  
  #filter for bins that satisfy setAgeCountPerBin
  if( !is.null(setAgeCountPerBin) ){
    if( length(setAgeCountPerBin==1) ){
      incl.bin <- FEATURE.BIN.MX[FEATURE.BIN.MX[,"countPerBin"]==setAgeCountPerBin,"bin"]
    } else if( length(setAgeCountPerBin==2) ){
      test1 <- as.numeric(FEATURE.BIN.MX[,"countPerBin"])>=setAgeCountPerBin[1]
      test2 <- as.numeric(FEATURE.BIN.MX[,"countPerBin"])>=setAgeCountPerBin[2]
      incl.bin <- FEATURE.BIN.MX[test1 & test2,"bin"]
    } else {
      stop("Invalid input for setAgeCountPerBin.")
    }
  }
 
  #initialize final table
  CONTACT.VAR.MX <- matrix( ncol=(2*length(feature))+3, nrow=(length(PERSIST.MX$ntis)),
                            dimnames=list(rown=NULL, 
                                          coln=c("i", "j", "ntis", final.feature)) )
  CONTACT.VAR.MX[,c("i", "j", "ntis")]<- cbind(i=as.numeric(PERSIST.MX$hits[,"i"]),
                                               j=as.numeric(PERSIST.MX$hits[,"j"]),
                                               ntis=PERSIST.MX$ntis)
  numOfCont <- length(PERSIST.MX$ntis)
  rm("PERSIST.MX"); gc()
  
  return( list(CONTACT.VAR.MX=CONTACT.VAR.MX, 
               #included contacts based on setAgeCountPerBin
               incl.ij.ind=which(CONTACT.VAR.MX[,"i"]%in%incl.bin & CONTACT.VAR.MX[,"i"]%in%incl.bin),
               drop.ij.perc=paste0( round(100-(length(incl.ij.ind)/numOfCont)*100, digits=2), "%")) )
}