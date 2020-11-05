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
  if( any(ij.df$j<=ij.df$i) ){ stop("Not upper triangle contacts only.") }
  gap.v <- ij.df$j-ij.df$i-1
  
  # Choose and mask contacts
  len.incl.x <- length(incl.bin.x)
  len.incl.y <- length(incl.bin.y)
  len.mask.x <- length(mask.bin.x)
  len.mask.y <- length(mask.bin.y)
  
  incl.TF <- rep(NA, times=nrow(ij.df))
  
  # In upper triangle, i -> y, j -> x
  
  # Mark selected contacts
  if( is.null(incl.bin.x) & is.null(incl.bin.y) ){ 
    
    incl.TF[is.na(incl.TF)] <- TRUE
    print("No selection of contacts...", quote=FALSE)
    
  } else if( !is.null(incl.bin.x) & !is.null(incl.bin.y) & 
             len.incl.x==len.incl.y ){
    
    for(i in 1:len.incl.x){
      incl.TF[ ij.df$i%in%incl.bin.y[[i]] & 
               ij.df$j%in%incl.bin.x[[i]] ] <- TRUE
    }
    print("Selecting contacts...", quote=FALSE)
    
  } else {
    stop("incl.bin.x and incl.bin.y should both be NULL or not NULL. 
         Length of list should be the same too.")
  }
  
  # Mark unwanted contacts
  if( !is.null(mask.bin.x) & !is.null(mask.bin.y) & 
      len.mask.x & len.mask.y ){
    
    for(i in 1:len.mask.x){
      incl.TF[ ij.df$i%in%mask.bin.y[[i]] & 
               ij.df$j%in%mask.bin.x[[i]] ] <- FALSE
    }
    print("Masking contacts...", quote=FALSE)
    
  } else if( is.null(mask.bin.x) & is.null(mask.bin.y) ) {
    
    print("No masking of contacts...", quote=FALSE)
    
  } else {
    stop("mask.bin.x and mask.bin.y should both be NULL or not NULL. 
         Length of list should be the same too.")
  }
  
  incl.TF[is.na(incl.TF)] <- FALSE
  if( !is.null(gap.range) & !identical(gap.range, c(0, Inf)) ){
    rm(ij.df)
    incl.TF[ gap.v<gap.range[1] | gap.v>gap.range[2] ] <- FALSE
    print("Gap filtering...", quote=FALSE)
  }

  return(incl.TF)
  
}
################################################################################
# rm(list=ls()); gc()
