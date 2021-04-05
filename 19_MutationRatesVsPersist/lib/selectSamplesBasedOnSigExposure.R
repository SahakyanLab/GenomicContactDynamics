################################################################################
# Select donors/samples based on a signature exposure threshold to a signature
# or set of signatures defined by Degasperi et al. 2020.
### FUNCTION ###################################################################
selectSamplesBasedOnSigExposure <- function(sigExposure.df='sigExposure.df',
                                            sig.v='c("RefSig.MMR1","RefSig.MMR2")',
                                            # Only consider samples/donor within 
                                            # this signature exposure percentage 
                                            # range (right closed interval).
                                            sigEperc.limits = 'c(min=50, max=100)'
){
  
  # 14 samples with no exposure to any signature (use as control?)
  
  print("Selecting donors/samples based on signature exposure limits...", 
        quote=F)
  
  sigEperc.limits <- sort(as.numeric(sigEperc.limits), decreasing=F)
  names(sigEperc.limits) <- c("min", "max")
  
  if( any(!sig.v%in%colnames(sigExposure.df)) ){
    stop("Invalid input signatures...")
  } else {
    print(sig.v, quote=F)
  }
  
  # If selecting based on >1 signatures, sample should satisfy sigEperc.limits for
  # both signatures.
  donor.mx <- data.matrix(sigExposure.df[,sig.v, drop=FALSE])
  donor.TF <- rowSums(x=donor.mx>as.numeric(sigEperc.limits["min"]) &
                        donor.mx<=as.numeric(sigEperc.limits["max"]), na.rm=F)
  if( any(is.na(donor.TF)) ){
    stop("Missing signature exposure values.")
  }
  if( any(donor.TF>length(sig.v)) | any(donor.TF<0) ){
    stop("Problem with selecting samples.")
  }
  donor.TF <- donor.TF==length(sig.v) & !is.na(sigExposure.df$alt.ID)
  
  donor <- sigExposure.df$alt.ID[donor.TF]
  
  if( any(duplicated(donor)) ){
    stop("Duplicated donor/s.")
  }
  
  return(donor)
  
}

# rm(list=ls()); gc()

