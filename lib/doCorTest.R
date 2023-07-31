################################################################################
# Function to perform correlation tests
### FUNCTION ###################################################################
doCorTest <- function(xval, yval, alt="two.sided", exactpval=NULL, out.dir, out.name,
                      method.names=c("pearson", "spearman")){ # kendall
  
  # Check method.names argument
  
  valid.methods <- c("pearson", "spearman", "kendall")
  if( any(!method.names %in% valid.methods) ){
    warning("doCorTest(): Valid method.names are pearson, spearman, kendall. Ignoring invalid names.")
  }
  
  method.names <- intersect(valid.methods, method.names)
  if( length(method.names) == 0 ){
    stop("doCorTest: No valid method.names.")
  }

  # Check missing values. If none, proceed. 
  
  tmp <- c(xval, yval)
  if( any(is.na(tmp)) | !is.numeric(tmp) ){
    
    # Used stop to prompt user to manually check missing values that can mess up value-group
    # pairing
    stop("doCorTest: Non-finite values in xval and/or yval or values not numeric.")
    
  } else {
    
    pear <- spea <- kend <- NULL
    
    if("pearson" %in% method.names){
      pear <- cor.test(x=xval, y=yval, method="pearson", alternative=alt, exact=exactpval)
    }
    if("spearman" %in% method.names){
      spea <- cor.test(x=xval, y=yval, method="spearman", alternative=alt, exact=exactpval)
    }
    if("kendall" %in% method.names){
      kend <- cor.test(x=xval, y=yval, method="kendall", alternative=alt, exact=exactpval)
    }
    
    TEST <- list(pear=pear, spea=spea, kend=kend, alt=alt)
    
    if( !is.null(out.dir) & !is.null(out.name) ){
      save(TEST, file=paste0(out.dir, "/", out.name, "_cortest.RData"))
    } 
    
    return(TEST)
    
  }
  
}

# rm(list=ls()); gc()