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
  MXsubj[lower.tri(MXsubj, diag=TRUE)] <- NA
  MXref[lower.tri(MXref, diag=TRUE)] <- NA
  
  # Choose and mask contacts
  ij.df <- expand.grid(1:len, 1:len)
  colnames(ij.df) <- c("i", "j")
  ij.df <- ij.df[ij.df$j>ij.df$i,]
  incl.TF <- filterContacts(ij.df=ij.df, gap.range=gap.range,
                            incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y,  
                            mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y)   
  # Contacts not included, to be marked NA
  ij.df <- ij.df[!incl.TF,]; rm(incl.TF)
  if(nrow(ij.df)>0){
    MXsubj[as.numeric(ij.df$i), as.numeric(ij.df$j)] <- NA
    MXref[as.numeric(ij.df$i), as.numeric(ij.df$j)] <- NA
  }
  rm(ij.df); gc()

  tot.ij.v <- NA.ij.v <- nonNA.inc.v <- c(subj=NA, ref=NA) 
  for( mx.id in c("subj", "ref") ){
    
    # Check for contacts NA in ref but not in subj and vice versa.
    # Ideally, there should be none.
    if(mx.id=="subj"){
      inc.TF <- !is.na(MXsubj) & is.na(MXref) 
    } else if(mx.id=="ref"){
      inc.TF <- is.na(MXsubj) & !is.na(MXref) 
    }
    eval(parse(text=paste0(
      "MX", mx.id, "[inc.TF] <- NA"
    )))
    nonNA.inc.v[mx.id] <- sum(inc.TF)
    rm(inc.TF)
    
    # Calculate final ij total
    eval(parse(text=paste0(
      "nonNA.TF <- !is.na( MX", mx.id, "[upper.tri(MX", mx.id, ", diag=FALSE)] )"
    )))
    tot.ij.v[mx.id] <- sum(nonNA.TF)
    # Should be 0, except for complementarity matrices
    NA.ij.v[mx.id] <- sum(!nonNA.TF)
    rm(nonNA.TF)
    
  }
  if( (tot.ij.v[1]!=tot.ij.v[2]) | (NA.ij.v[1]!=NA.ij.v[2]) ){
    stop("Total nonNA ij different for subject and matrices.")
  }
  tot.ij.v <- unique(tot.ij.v)
  NA.ij.v <- unique(NA.ij.v)
  if( (tot.ij.v+NA.ij.v)!=((len*len-len)/2) ){
    stop("Final ij and NA count dont add up to length of upper matrix.")
  }
  #-------------------Make cut-off matrix for grid search
  c.off.mx <- data.matrix(expand.grid(x=c.offsubj.v, y=c.offref.v))
  dimnames(c.off.mx)[[2]] <- c("c.offsubj", "c.offref")
  c.off.mx.len <- nrow(c.off.mx)
  rm(c.offsubj.v, c.offref.v)
  
  toExport <- c("c.off.mx", "MXsubj", "MXref", "tot.ij.v")
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
      # Ref positive and negative
      RP <- sum(mxref=="y", na.rm=TRUE)
      RN <- sum(mxref=="n", na.rm=TRUE)
      # Subj positive and negative
      SP <- sum(mxsubj=="y", na.rm=TRUE)
      SN <- sum(mxsubj=="n", na.rm=TRUE)
      
      # True positive
      TP <- sum(mxsubj=="y" & mxref=="y", na.rm=TRUE)
      # False positive
      FP <- sum(mxsubj=="y" & mxref=="n", na.rm=TRUE)
      # True negative
      TN <- sum(mxsubj=="n" & mxref=="n", na.rm=TRUE)
      # False negative
      FN <- sum(mxsubj=="n" & mxref=="y", na.rm=TRUE)
      rm(mxsubj, mxref); gc()
      
      if( (TP+FP)!=SP | (TN+FN)!=SN | 
          (TP+FN)!=RP | (TN+FP)!=RN |
          (SP+SN)!=tot.ij.v | (RP+RN)!=tot.ij.v ){
        stop("Important confusion matrix numbers don't add up.")
      }
      
      v <- c(c.offsubj=c.offsubj, c.offref=c.offref, 
             SP=SP, RP=RP, TP=TP, FP=FP, 
             SN=SN, RN=RN, TN=TN, FN=FN)
      return(v)
      
    }) # chunk sapply end
    
    return(do.call("rbind", chunk))
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  out <- cbind(nonNA.subj.NA.ref.ij=nonNA.inc.v["subj"],
               nonNA.ref.NA.subj.ij=nonNA.inc.v["ref"],
               final.NA.ij=NA.ij.v, final.nonNA.ij=tot.ij.v, out)
  return(out)
  
}
################################################################################
compareContactMx <- cmpfun(compareContactMx, options=list(suppressUndefined=TRUE))
################################################################################
#  rm(list=ls()); gc()



