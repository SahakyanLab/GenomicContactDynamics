############################################################################### 
# Make boxplot of shared number of repeat elements vs. Cp
# source(paste0(lib, "/compareTwoDist.R"))
# source(paste0(lib, "/doVarTest.R"))
# source(paste0(lib, "/doCorTest.R"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/compareTwoDist.R"))
### FUNCTION ###################################################################
makeMinRepPlot <- function(MINREPCOUNTS, ntis.v, 
                           affix=paste0(chr, "_", gcb, "_", out.name),
                           out.dir,
                           nCPU,
                           metric,
                           generatePREELMTISSDYNonly=FALSE
                           ){
  
  if( any(lengths(MINREPCOUNTS)!=length(ntis.v)) ){
    stop("Missing Cp in an element.")
  } else {
    
    elements <- names(MINREPCOUNTS)
    elements.len <- length(elements)
    
    #for(elm.ind in 1:elements.len){
    
    toExport <- c("elements", "MINREPCOUNTS", "ntis.v",  "affix", "out.dir", "metric")
    
    #### PARALLEL EXECUTION #########
    
    # List of lists
    PREELMTISSDYN.MX <- foreach( itr=isplitVector(1:elements.len, chunks=nCPU), .combine="rbind", 
                                 .export=toExport, .noexport=ls()[!ls()%in%toExport]
                         
    ) %op% {
      
      mean.val <- list()
      median.val <- list()
      
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
                                          mincount=unlist( mapply(rep, x=mincount, times=v, SIMPLIFY=FALSE) ),
                                          stringsAsFactors=FALSE)
                         
                       }
                       
                     })
        
        rm(lst)
        df <- do.call("rbind", df)
        
        if( is.null(df) ){
          
          print(paste0(elm, " makeMinRepPlot(): Checkpoint 1."))
          next
          
        } else if( !is.null(df) & ( ncol(df)==1 | any(!is.finite(df$mincount)) ) ){
          
          print(paste0(elm, " makeMinRepPlot(): Checkpoint 2."))
          next
          
        }
        
        colnames(df) <- c("cp", "mincount")
        
        affix1 <- paste0( affix, "_", gsub(pattern="[^[:alnum:][:space:]]", 
                                           replacement="", x=elm) )
        affix1 <- paste0(affix1, "_", elm.ind)
        
        if( !generatePREELMTISSDYNonly ){
          
          print("Plotting...", quote=FALSE)
          
          #cp.v <- sort(unique(df$cp), decreasing=F)
          TEST <- try(compareTwoDist(x=df$mincount[df$cp%in%1:3],
                                     y=df$mincount[df$cp%in%19:21]))
          
          if( is.atomic(TEST) ){
            TEST <- list(test.id="failed")
          }
          
          save(TEST, file=paste0(out.dir, "/", affix1, "_testresultinboxplot.RData"))
          
          df$cp <- factor(as.character(df$cp), levels=as.character(1:21))
          
          # Other statistical tests
          try(doVarTest(xval=df$mincount, grp=df$cp, out.dir=paste0(out.dir, "/correlation"),  
                        out.name=affix1))
          try(doCorTest(xval=as.numeric(as.character(df$cp)), yval=df$mincount, 
                        alt="two.sided", exactpval=F, out.dir=paste0(out.dir, "/correlation"), 
                        out.name=paste0(affix1, "_cortest.RData")))
                    
          # Plot
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
          
          #rm(df)
          
        }
        
        # Calculate PREELMTISSDYN.MX
        
        #if(metric=="_skewrep"){
        #  excl.TF <- df$mincount == 1
        #} else if(metric==""){
        #  excl.TF <- df$mincount == 0
        #}
        
        write.csv(table(df$cp), file=paste0(out.dir, "/", affix1, "_DataPerCp.csv"), row.names=FALSE)
        
        cent <- by(data=df$mincount, INDICES=df$cp, simplify=FALSE, FUN=function(mincountsPcp){ 
          
          # With all NAs, mean(na.rm=TRUE) returns NaN while median(na.rm=TRUE) returns NA
          c(MEAN=mean(mincountsPcp, na.rm=FALSE), MEDIAN=median(mincountsPcp, na.rm=FALSE) )
          
        })
        
        cent[lengths(cent)==0] <- rep(x=list(c(MEAN=NA, MEDIAN=NA)))
        
        if( identical( names(cent), as.character(1:21)) ){
          
          cent <- do.call("rbind", cent)
          mean.val[[elm]] <- cent[,"MEAN"]
          median.val[[elm]] <- cent[,"MEDIAN"]
          
        } else {
          
          rm(df)
          stop(paste0(elm, ": CENTRAL checkpoint."))
          
        }
        
        print(paste0(elm, " done!"), quote=FALSE)
        
      } # itr for loop end
      
      mean.val <- do.call("rbind", mean.val)
      median.val <- do.call("rbind", median.val)
      median.val[is.na(median.val)] <- NaN
      
      return(mean.val)
    
    } 
    
    ### END OF PARALLEL EXECUTION ###
    
    PREELMTISSDYN.MX <- rbind(num.contact=rep(1), PREELMTISSDYN.MX)
    return(PREELMTISSDYN.MX)
    
  }
  
}