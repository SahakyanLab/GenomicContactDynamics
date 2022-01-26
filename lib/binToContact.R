################################################################################
# Use contact bin values to derive contact values.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
# library(foreach)
# library(doParallel)
# library(itertools)
# source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
binToContact <- function(v = 'vector of values named with bin number',
                         binval.str.TF = 'logical, elements of v are strings? 
                                          If yes, it should be string of numeric 
                                          values separated by ";"',
                         missing.val = 'values for bins without values or bins not 
                                        in v',
                         funx.v = 'string, how to combine bin values of a contact',
                         funx.id.v = 'column name containing output of each function
                                      in funx.v; should be same order as funx.v',
                         mxcalc.TF = 'logical, T if functions in funx.v do row 
                                      matrix calculation',
                         persist.dir = 'directory of PERSIST.MX',
                         gcb = 'min2Mb or min05Mb',
                         chr = 'chromosome',
                         ct = 'cell/tissue, use All to use all contacts',
                         tot.bin = 'total bins of a chr',
                         nCPU # contacts
                         ){
  
  run.id <- paste0(gcb, "_", chr, "_", ct, ":")
  funx.v.len <- length(funx.v)
  
  # CHECKS
  if( funx.v.len!=length(funx.id.v) ){
    stop("funx.v and funx.id.v not equal in length.")
  }
  if( is.null(names(v)) ){
    stop(paste0(run.id, " v not named."))
  }
  if( any(names(v)=="") ){
    stop(paste0(run.id, " Some values in v not named."))
  }
  if( any(duplicated(names(v))) ){
    stop(paste0(run.id, " Names of v duplicated."))
  }
  
  # Prepare bin values for All bins of chr
  if(binval.str.TF){
    
    v <- strsplit(x=v, split=";", fixed=T)
    v <- lapply(X=v, FUN=as.numeric)
    
  } else {
    
    if( !is.numeric(v) ){
      stop(paste0(run.id, " v not a numeric vector."))
    }
     v <- as.list(v) 
  }
  
  V <- rep(x=missing.val, times=tot.bin)
  names(V) <- 1:tot.bin
  V <- as.list(V)
  V[names(v)] <- v
  rm(v); gc()
 
  # Get contacts
  #print("PERSIST.MX loading...")
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  PERSIST.MX$control <- PERSIST.MX$valsum <- NULL
  #print("PERSIST.MX loaded.")
  ct.v <- colnames(PERSIST.MX$hits)[!colnames(PERSIST.MX$hits)%in%c("i","j")]
  
  if(ct%in%ct.v){
    
    ij.TF <- PERSIST.MX$hits[[ct]]>0 
    ij.df <- data.frame(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, 
                        Cp=PERSIST.MX$ntis, stringsAsFactors=F)[ij.TF,]
    print(paste0(run.id, " Using only ", ct, " contacts."))
    rm(ij.TF)
    
  } else if(ct=="All"){
    
    ij.df <- data.frame(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, 
                        Cp=PERSIST.MX$ntis, stringsAsFactors=F)
    print(paste0(run.id, " Using ", ct, " contacts."))
    
  } else {
    stop(paste0(run.id, " Invalid ct."))
  }
  rm(PERSIST.MX, ct.v); gc()
  
  # Calculate contact values
  
  if(!mxcalc.TF){ # Loop calculation but with parallel computing
    
    ij.len <- length(ij.df[,1])
    
    toExport <- c("V", "ij.df", "funx.v")
    #### PARALLEL EXECUTION #########
    values <- foreach(itr=isplitVector(1:ij.len, chunks=nCPU), 
                      .inorder=T, .combine="rbind",
                      .export=toExport, .noexport=ls()[!ls()%in%toExport]
                      
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=F, FUN=function(n){
        
        val <- c( V[[as.character(ij.df$i[n])]], V[[as.character(ij.df$j[n])]] )
        a <- sapply(X=funx.v, simplify=T, FUN=function(funx){
          eval(parse(text=paste0(
            funx, '(c(', paste(val, collapse=","), '))'
          )))
        })
        return(a)
        
      })
        
      return(do.call("rbind", chunk))
      
    }
    ### END OF PARALLEL EXECUTION ###
    
    if(length(values[,1])!=ij.len){
      stop(paste0(run.id, " Error in getting contact values"))
    }
    rm(toExport)
    dimnames(values)[[2]] <- funx.id.v
    
    ij.df <- cbind.data.frame(ij.df, values)
    rm(values)
    
  } else if(mxcalc.TF & !binval.str.TF){ # Vectorised calculation
    
    values <- cbind(as.numeric(V[as.character(ij.df$i)]), 
                    as.numeric(V[as.character(ij.df$j)])
                    )
    
    for(f in 1:funx.v.len){
      
      eval(parse(text=paste0(
        #'ij.df$', funx.id.v[f], ' <- as.numeric(', funx.v[f], 'values))'
        'ij.df$', funx.id.v[f], ' <- as.numeric(', funx.v[f], '(values))'
      )))
      
    } # funx.v.len for loop end
    
  } else {
    stop("Invalid mxcalc.TF and binval.str.TF arguments!")
  }

  return(ij.df)
  
}
################################################################################
binToContact <- cmpfun(binToContact, options=list(suppressUndefined=TRUE))
################################################################################

# rm(list=ls()); gc()