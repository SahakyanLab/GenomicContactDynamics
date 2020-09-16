################################################################################
# Calculate true positive, negative and false positive, negative rates of a 
# contact matrix based on a reference contact matrix. 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
# library(foreach)
# library(doParallel)
# library(itertools)
# lib = "/Users/ltamon/DPhil/lib"
# source(paste0(lib, "/UTL_doPar.R"))
# source(paste0(wk.dir, "/lib/filterContactsInMx.R"))
### FUNCTION ###################################################################
compareContactMx <- function(MXsubj, MXref, c.offsubj.v, c.offref.v,
                             # Upper triangle as reference for list of indices
                             incl.bin.x = NULL, # list, 
                             incl.bin.y = NULL, # list, 
                             mask.bin.x = NULL, # list
                             mask.bin.y = NULL, # list
                             gap.range = NULL,  # vector
                             nCPU = 1L){
  
  #-------------------Check and process matrices
  if( ( !is.matrix(MXsubj) | !is.matrix(MXref) ) ){
    stop("Not matrices.")
  }
  if( !is.numeric(MXsubj) | !is.numeric(MXref) | !is.numeric(c.offsubj.v) | !is.numeric(c.offref.v) ){
    stop("Matrices or cut-off vectors not numeric.")
  }
  if( !identical(dim(MXsubj), dim(MXref)) ){
    stop("Dimension of matrices different.")
  }
  
  len  <- unique(dim(MXsubj))
  if(length(len)>1){ stop("Matrix not symmetric.") }
  dimnames(MXsubj) <- dimnames(MXref) <- list(1:len, 1:len)
  
  # Retain only upper triangle
  MXsubj[lower.tri(MXsubj, diag=TRUE)] <- NA_integer_
  MXref[lower.tri(MXref, diag=TRUE)] <- NA_integer_
  
  # Choose and mask contacts
  ij.df <- expand.grid(1:len, 1:len)
  colnames(ij.df) <- c("i", "j")
  ij.df <- ij.df[ij.df$j>ij.df$i,]
  incl.TF <- filterContacts(ij.df=ij.df, gap.range=gap.range,
                            incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y,  
                            mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y)   
  
  ij.df <- ij.df[!incl.TF,]; rm(incl.TF)
  MXsubj[as.numeric(ij.df$i), as.numeric(ij.df$j)] <- NA
  rm(ij.df); gc()
  
  # Check for contacts NA in ref but not in subj. Nothing to compare to.
  nonNA.subj.NA.ref.ij <- !is.na(MXsubj) & is.na(MXref) 
  MXsubj[nonNA.subj.NA.ref.ij] <- NA_integer_
  nonNA.subj.NA.ref.ij <- sum(nonNA.subj.NA.ref.ij)
  
  # Other important numbers
  nonNA.subj.TF <- !is.na( MXsubj[upper.tri(MXsubj, diag=FALSE)] )
  final.subj.tot.ij <- sum(nonNA.subj.TF)
  NA.subj.ij <- sum(!nonNA.subj.TF)
  rm(nonNA.subj.TF)
  
  #-------------------Make cut-off matrix for grid search
  c.off.mx <- data.matrix(expand.grid(x=c.offsubj.v, y=c.offref.v))
  dimnames(c.off.mx)[[2]] <- c("c.offsubj", "c.offref")
  c.off.mx.len <- nrow(c.off.mx)
  rm(c.offsubj.v, c.offref.v)
  
  toExport <- c("c.off.mx", "MXsubj", "MXref", "final.subj.tot.ij")
  #### PARALLEL EXECUTION #########
  out <- foreach(itr=isplitVector(1:c.off.mx.len, chunks=nCPU), 
                 .inorder=TRUE, .combine="rbind",
                 .export=toExport, .noexport=ls()[!ls()%in%toExport]
          
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      c.offsubj <- unname(c.off.mx[i,"c.offsubj"])
      c.offref <- unname(c.off.mx[i,"c.offref"])
      
      #  Modify matrices based on cut-off for contact
      mxsubj <- MXsubj; mxref <- MXref
      mxsubj[MXsubj<=c.offsubj] <- "n"; mxsubj[mxsubj!="n"] <- "y" 
      mxref[MXref<=c.offref] <- "n"; mxref[mxref!="n"] <- "y"
      
      #-------------------Calculate
      P <- sum(mxsubj=="y", na.rm=TRUE)
      N <- sum(mxsubj=="n", na.rm=TRUE)
      
      # True positive
      TP <- sum(mxsubj=="y" & mxref=="y", na.rm=TRUE)
      # False positive
      FP <- sum(mxsubj=="y" & mxref=="n", na.rm=TRUE)
      # True negative
      TN <- sum(mxsubj=="n" & mxref=="n", na.rm=TRUE)
      # False negative
      FN <- sum(mxsubj=="n" & mxref=="y", na.rm=TRUE)
      rm(mxsubj, mxref); gc()
      
      if( (TP+FP)!=P | (TN+FN)!=N | (P+N)!=final.subj.tot.ij ){
        stop("Important numbers don't add up.")
      }
      
      v <- c(c.offsubj=c.offsubj, c.offref=c.offref, 
             P=P, TP=TP, FP=FP, N=N, TN=TN, FN=FN)
      return(v)
      
    }) # chunk sapply end
    
    return(do.call("rbind", chunk))
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  out <- cbind(NA.subj.ij=NA.subj.ij, nonNA.subj.NA.ref.ij=nonNA.subj.NA.ref.ij, 
               final.subj.tot.ij=final.subj.tot.ij, out)
  return(out)
  
}
################################################################################
compareContactMx <- cmpfun(compareContactMx, options=list(suppressUndefined=TRUE))
################################################################################
#  rm(list=ls()); gc()



