################################################################################
# Per chr, return anv including anv distance value (default is euclidean) per contact
### FUNCTION ###################################################################
feat_anv <- function(ANV.MX, contact.mx = 'i-j 2-column matrix', 
                     distanceOnly = F, 
                     distFUN=function(a, b){ sqrt(sum((a - b)^2)) } # https://www.statology.org/euclidean-distance-in-r/
                     ){
  
  is_wrong <- c(any(!is.numeric(contact.mx[,1])),
                any(!is.numeric(contact.mx[,2])),
                any(!is.numeric(ANV.MX[,"bins"]))
  )
  if( any(is_wrong) ){
    rm(ANV.MX)
    stop("feat_anv(): Checkpoint.")
  }
  
  #
  i.anvInd <- match( contact.mx[,1], table=ANV.MX[,"bins"] )
  j.anvInd <- match( contact.mx[,2], table=ANV.MX[,"bins"] )
  
  ANV.MX <- ANV.MX[,setdiff(dimnames(ANV.MX)[[2]], c("bins", "startpos", "endpos", "numUMChar"))]
  ij.len <- length(contact.mx[,1])
  
  #
  distances <- sapply(X=1:ij.len, simplify=T, FUN=function(ij.ind){
    do.call(distFUN, list(a=ANV.MX[i.anvInd[ij.ind],], 
                          b=ANV.MX[j.anvInd[ij.ind],]))
  })
  
  #
  if(distanceOnly){
    return(distances)
  } else {
    
    feat.mx <- cbind(ANV.MX[i.anvInd,], ANV.MX[j.anvInd,])
    anv.nme <- dimnames(ANV.MX)[[2]]
    dimnames(feat.mx) <- list(dimnames(contact.mx)[[1]], 
                              c(paste0("i_", anv.nme),  paste0("j_", anv.nme)))
    feat.mx <- cbind(ij_anvdist=distances, feat.mx)
    
    return(feat.mx)
    
  }
  
}

# rm(list=ls()); gc()
