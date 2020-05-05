subsetANDbinTableRowCol <- function(table = MINELM.MX[MINELM.MX[,"ntis"]==21,-1],
                 # vector or NULL if no need for subsetting
                 toSubsetANDbin = unique( fread(file=Giordano364file, header=FALSE, data.table=FALSE, 
                                          stringsAsFactors=FALSE)[[1]] ),
                 defineOrder = TRUE, row = FALSE, binBy = 50 ){

  # Transpose if row = FALSE because succeding code is made to bin rows
  if(row==FALSE){ 
    oldcolNames <- rownames(table)
    table <- t(table)
    print("Subsetting/Binning on columns.")
  } else {
    oldcolNames <- colnames(table)
    print("Subsetting/Binning on rows.")
  }
 
  # Subset table if required
  if( !is.null(toSubsetANDbin[1]) ){
    table <- table[ rownames(table)%in%toSubsetANDbin, ]
    print("Table subsetted.")
    if(defineOrder==TRUE){ 
      table <- table[ order( match(rownames(table), toSubsetANDbin) ), ] 
      print("Table ordered based on toSubsetANDbin.")
    } else {
      print("Original table order kept.")
    }
  }
  
  #----------BINNNING-START----------
  
  ## Two ways to define boundaries
  
  if(length(binBy)==1){
    ##### [1] BIN BY A FIXED LENGTH
    bounds.start <- seq(from=1, to=nrow(table), by=binBy)
    bounds.end <- c(bounds.start[-1]-1, nrow(table))
  
  } else if( length(binBy)>1 ){
    ##### [2] BIN BY A VECTOR OF GROUPS
    
    # this is needed when some toSubsetANDbin elements are not in table
    # toSubsetANDbin and binBy should correspond
    # order kept
    binBy <- binBy[ !is.na( match(toSubsetANDbin, rownames(table)) ) ]
    
    bounds.start <- c()
    counter <- 0L
    binBy.len <- length(binBy)
    for(k in 2:binBy.len){
      if( binBy[k]!=binBy[k-1] ){
        counter <- counter + 1L
        bounds.start[counter] <- k
      }
    }
    
    bounds.end <- c(bounds.start-1, binBy.len)
    bounds.start <- c(1, bounds.start)
    
  } else {
    stop("Invalid input for binBy.")
  }
  
  if( length(bounds.start)==length(bounds.end) ){
    Numbins <- length(bounds.start)
  } else {
    stop("Checkpoint 1")
  } 
  
  # initialize binned table
  binned.table <- matrix(ncol=ncol(table), nrow=Numbins)
  for( bin in 1:length(bounds.start) ){
    #make sure dim > 2
    #BUG! different result if dim > 1 or 2
    #subtable <- rbind(rep(NA), 
    #                  table[bounds.start[bin]:bounds.end[bin],])
    subtable <- matrix(table[bounds.start[bin]:bounds.end[bin],],
                       nrow=binBy, ncol=length(oldcolNames))
    #did not use colSums because CN were converted to characters
    #binned.table[bin,] <- apply(X=subtable, MARGIN=2, 
    #                            FUN=function(x) {sum(as.numeric(x), 
    #                                                 na.rm=TRUE)} )
    binned.table[bin,] <- apply(X=subtable, MARGIN=2, 
                                FUN=function(x){mean(as.numeric(x))} 
                                )
    #binned.table[bin,] <- colSums(x=subtable, na.rm=TRUE)
  }
  
  #----------BINNNING-END------------
  
  # add column names and row names
  if( !is.null(oldcolNames) ){ colnames(binned.table) <- oldcolNames }
  if( length(binBy)==1 ){
    rownames(binned.table) <- paste("Bin", 1:Numbins, sep="")
  } else {
    rownames(binned.table) <- binBy[bounds.start]
  }
  
  # transpose back if required
  if(row==FALSE){ binned.table <- t(binned.table) }
  
  #---------------------
  return(binned.table)
}
  

