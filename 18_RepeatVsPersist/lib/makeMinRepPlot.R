############################################################################### 
# Make boxplot of shared number of repeat elements vs. Cp
# source(paste0(lib, "/compareTwoDist.R"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/compareTwoDist.R"))
### FUNCTION ###################################################################
makeMinRepPlot <- function(MINREPCOUNTS, ntis.v, 
                           affix=paste0(chr, "_", gcb, "_", out.name),
                           out.dir,
                           nCPU
                           ){
  
  if( any(lengths(MINREPCOUNTS)!=length(ntis.v)) ){
    stop("Missing Cp in an element.")
  } else {
    
    print("Plotting...", quote=FALSE)
    
    elements <- names(MINREPCOUNTS)
    elements.len <- length(elements)
    
    #for(elm.ind in 1:elements.len){
    
    toExport <- c("elements", "MINREPCOUNTS", "ntis.v",  "affix", "out.dir")
    
    #### PARALLEL EXECUTION #########
    
    # List of lists
    foreach( itr=isplitVector(1:elements.len, chunks=nCPU),
             .export=toExport, .noexport=ls()[!ls()%in%toExport]
                         
    ) %op% {
      
      for(elm.ind in itr){
        
        elm <- elements[elm.ind]
        lst <- MINREPCOUNTS[[elm]]
        df <- sapply(X=names(lst), simplify=FALSE, 
                     FUN=function(ntis){
                       v <- unname(lst[[ntis]])
                       mincount <- as.numeric( names(lst[[ntis]]) )
                       if( length(mincount)==0 ){
                         return(NULL)
                       } else {
                         
                         df <- data.frame(ntis=rep(as.numeric(ntis)),
                                          mincount=unlist( mapply(rep, x=mincount, times=v) ),
                                          stringsAsFactors=F)
                         
                       }
                       
                     })
        
        rm(lst)
        df <- do.call("rbind", df)
        
        if( ncol(df)==1 ){
          print(paste0(elm, ": Not present in any contact."))
          next
        }
        
        colnames(df) <- c("cp", "mincount")
        affix1 <- paste0( affix, "_", gsub(pattern="[^[:alnum:][:space:]]", 
                                           replacement="", x=elm) )
        affix1 <- paste0(affix1, "_", elm.ind)
        
        #cp.v <- sort(unique(df$cp), decreasing=F)
        TEST <- try(compareTwoDist(x=df$mincount[df$cp%in%1:3],
                                   y=df$mincount[df$cp%in%19:21]))
        
        if( is.atomic(TEST) ){
          TEST <- list(test.id="failed")
        }
        
        df$cp <- factor(as.character(df$cp), levels=as.character(1:21))
        
        png(file=paste0(out.dir, "/", affix1, "_bp.png"), 
            res=300, width=3000, height=3000)
        boxplot(formula=mincount~cp, data=df, col="honeydew3", 
                cex.lab=1.3, cex.axis=1.3,
                outline=FALSE, xlab="", ylab="", main="")
        # x axis
        mtext(side=1, text=expression("c"["p"]), line=3, cex=1.5)
        # y axis
        mtext(side=2, text="Contact min. repeat count", line=2.7, cex=1.5)
        # Diagram title
        mtext(side=3, text=paste0(affix1, TEST$test.id, "_x=Cp1To3,y=Cp19To21"), line=1.5, cex=0.5)
        dev.off()
        
        rm(df)
        
        print(elm)
        
      }
      
    }
    
    ### END OF PARALLEL EXECUTION ###
    
  }
  
}