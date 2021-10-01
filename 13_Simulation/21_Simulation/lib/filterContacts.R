################################################################################
# Filter contacts on a symmetric matrix of length L
### FUNCTION ###################################################################
filterContacts <- function(ij.df=ij.df, # i-j dataframe, upper triangle
                           gap.range=gap.range,
                           incl.bin.x=NULL, incl.bin.y=NULL, # list, 
                           mask.bin.x=NULL, mask.bin.y=NULL  # list 
                         ){
  
  ij.df <- ij.df[,1:2]
  colnames(ij.df) <- c("i", "j")
  if( any(ij.df$j<=ij.df$i) ){
    stop("filterContacts(): Not upper triangle contacts only.")
  }
  gap.v <- ij.df$j-ij.df$i-1
  
  # Choose and mask contacts
  len.incl.x <- length(incl.bin.x)
  len.incl.y <- length(incl.bin.y)
  len.mask.x <- length(mask.bin.x)
  len.mask.y <- length(mask.bin.y)
  
  incl.TF <- rep(NA, times=nrow(ij.df))
  
  # In upper triangle perspective, i -> y, j -> x
  
  # Mark selected contacts
  if( is.null(incl.bin.x) & is.null(incl.bin.y) ){ 
    
    incl.TF[is.na(incl.TF)] <- TRUE
    print("filterContacts(): No selection of contacts...", quote=FALSE)
    
  } else if( !is.null(incl.bin.x) & !is.null(incl.bin.y) & 
             len.incl.x==len.incl.y ){
    
    for(i in 1:len.incl.x){
      incl.TF[ ij.df$i%in%incl.bin.y[[i]] & 
               ij.df$j%in%incl.bin.x[[i]] ] <- TRUE
    }
    print("filterContacts(): Selecting contacts with i in incl.bin.y and j
          in incl.bin.x...", quote=FALSE)
    
  } else if( !is.null(incl.bin.x) & is.null(incl.bin.y) ){
    
    # Any contacts with bins in incl.bin.x included
    for(i in 1:len.incl.x){
      incl.TF[ ij.df$i%in%incl.bin.x[[i]] | 
               ij.df$j%in%incl.bin.x[[i]] ] <- TRUE
    }
    print("filterContacts(): Selecting contacts with i or j in 
          incl.bin.x...", quote=FALSE)
  
  } else {
    stop("filterContacts(): incl.bin.x and incl.bin.y should both be NULL 
         or both not NULL or incl.bin.x not null and incl.bin.y NULL. 
         Length of list should be the same if both not NULL.")
  }
  
  # Mark unwanted contacts
  if( !is.null(mask.bin.x) & !is.null(mask.bin.y) & 
      len.mask.x & len.mask.y ){
    
    for(i in 1:len.mask.x){
      incl.TF[ ij.df$i%in%mask.bin.y[[i]] & 
               ij.df$j%in%mask.bin.x[[i]] ] <- FALSE
    }
    print("filterContacts(): Masking contacts with i in mask.bin.y and j 
          in mask.bin.x...", quote=FALSE)
    
  } else if( is.null(mask.bin.x) & is.null(mask.bin.y) ) {
    
    print("filterContacts(): No masking of contacts...", quote=FALSE)
    
  } else if( !is.null(mask.bin.x) & is.null(mask.bin.y) ){
    
    # Any contacts with bins in mask.bin.x excluded
    for(i in 1:len.mask.x){
      incl.TF[ ij.df$i%in%mask.bin.x[[i]] | 
               ij.df$j%in%mask.bin.x[[i]] ] <- FALSE
    }
    print("filterContacts(): Masking contacts with i or j in mask.bin.x...", 
          quote=FALSE)
    
  } else {
    stop("filterContacts(): mask.bin.x and mask.bin.y should both be NULL 
         or both not NULL or mask.bin.x not null and mask.bin.y NULL. 
         Length of list should be the same if both not NULL.")
  }
  
  incl.TF[is.na(incl.TF)] <- FALSE
  if( !is.null(gap.range) & !identical(gap.range, c(0, Inf)) ){
    rm(ij.df)
    incl.TF[ gap.v<gap.range[1] | gap.v>gap.range[2] ] <- FALSE
    print("filterContacts(): Gap filtering...", quote=FALSE)
  }

  return(incl.TF)
  
}
################################################################################
# rm(list=ls()); gc()
