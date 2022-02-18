################################################################################
# Function to perform correlation tests
### FUNCTION ###################################################################
doCorTest <- function(xval, yval, alt="two.sided", exactpval=T, out.dir, out.name){
  
  tmp <- c(xval, yval)
  if( any(is.na(tmp)) | !is.numeric(tmp) ){
    
    # Used stop to prompt user to manually check missing values that can mess up value-group
    # pairing
    stop("doCorTest: Non-finite values in xval and/or yval or values not numeric.")
    
  } else {
    
    pear <- cor.test(x=xval, y=yval, method="pearson", alternative=alt, exact=exactpval)
    spea <- cor.test(x=xval, y=yval, method="spearman", alternative=alt, exact=exactpval)
    
    TEST <- list(pear=pear, spea=spea, alt=alt)
    save(TEST, file=paste0(out.dir, "/", out.name, "_cortest.RData"))
    
  }
  
}

# rm(list=ls()); gc()