################################################################################
# Function to plot density of feature at overlapping radial windows of the model.
# Densities of features are normalized to the amount of DNA at each window.
# Dependencies:
## library(compiler)
## library(foreach)
## library(itertools)
## library(doParallel)
## library(reshape2)
## library(ggplot2)
## library(RColorBrewer)
## source(paste0(lib, "/GG_bgr.R"))
## source(paste0(lib, "/BINorSLIDE.R"))
################################################################################
################################################################################
plotFeatDensPerRadWindow <- function(
  out.dir = "dir",
  out.name = paste0(celltiss, "_", ploidy, "_", m),
  DOMXYZRFEATfile = "filepath",
  # periphery=1, bead constrained to periphery
  # periphery=0, bead not constrained to periphery
  periphery = c(1,0),
  foi.list = as.list( as.character(c(1:21)) ),
  feat.type = expression("c"["p"]),
  foiRData = FALSE,
  # Number of features
  nCPU=21L,
  # SLIDE
  dividingR="SLIDE",
  # Set radius of bin/window
  dr = 0.5,
  multiplier = 100L,
  # Cannot generate plot for raw counts because it exceeds the limit
  count.v = c("raw", "norm"),
  # Pick hist.breaks that will not generate 0 densities
  # Should cover all the rVal and may vary between Chrom3D structure
  # hist.breaks <- seq(from=0, to=max(RPCOUNT.DF$rVal), by=hist.breaks)
  hist.breaks = 0.5,
  regenerateCountData = TRUE,
  regeneratePlotData = TRUE
){
  
  # Initialise the nCPU-driven (number of CPUs) setting for do/dopar
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  # Load foi.list if RData
  if(foiRData==TRUE){
    
    file.nme <- load(file=foi.list)
    eval(parse(text=paste0(
      "foi.list <- ", file.nme, "; rm(", file.nme, "); gc()" 
    )
    ))
    
  }
  
  coul0 <- brewer.pal(11, "Spectral")
  coul <- colorRampPalette(rev(coul0))(length(foi.list))
  
  if(regenerateCountData==TRUE){
    
    # Load DOMXYZRFEAT.DF
    file.nme <- load(file=DOMXYZRFEATfile)
    eval(parse(text=paste0(
      "DTA <- ", file.nme, "; rm(", file.nme, "); gc()" 
    )
    ))
    
    DTA <- DTA[DTA$periphery%in%periphery,]
    print(paste0("Periphery: ", paste(periphery, collapse="&"), "; ", 
                 nrow(DTA), " beads."), quote=FALSE)
    
    rmin <- min(DTA$radDist, na.rm=TRUE)
    rmax <- max(DTA$radDist, na.rm=TRUE)
    
    DTA$domainlen <- abs(as.numeric(DTA$end-DTA$start))
    
    foiNme.v <- names(foi.list)
    
    # Generate reference points (midpoints) of windows
    rVal <- BINorSLIDE(numVec=DTA$radDist, action="SLIDE", 
                       numPoints=1000, dr=dr)
    
    rVal.len <- length(rVal)
    
    #### PARALLEL EXECUTION #########
    
    RPCOUNT.DF <- foreach( itr=isplitVector(1:rVal.len, chunks=nCPU), 
                           .combine="rbind", .inorder=TRUE,
                           .export=c("DTA", "rVal", "foiNme.v"), 
                           .noexport=ls()[!ls()%in%c("DTA", "rVal", "foiNme.v")] 
    ) %op% {
      RPCOUNT.DF.CHUNK <- sapply(X=itr, simplify=FALSE, FUN=function(i){
        
        # Subset beads/rows that lie within the window
        dta <- DTA[( DTA$radDist >= (rVal[i]-dr) & DTA$radDist <= (rVal[i]+dr) ),]
        
        totDNA <- sum(as.numeric(dta$domainlen), na.rm=TRUE)
        
        lst <- sapply(X=foiNme.v, simplify=FALSE, FUN=function(foiNme){
          foi <- foi.list[[foiNme]]
          foiCNT <- as.numeric( unlist( dta[,foi], use.names=FALSE) )
          totfoiCNT <- sum(foiCNT, na.rm=TRUE)
          data.frame( foi=as.character(foiNme), raw=totfoiCNT, 
                      norm=totfoiCNT/totDNA, stringsAsFactors=FALSE )
        })
        cbind.data.frame( rVal=rep(rVal[i]), do.call("rbind.data.frame", lst))
      })
      return(do.call("rbind", RPCOUNT.DF.CHUNK))
    }
    
    ### END OF PARALLEL EXECUTION ###
    
    rm(rVal, DTA); gc()
    
    rownames(RPCOUNT.DF) <- NULL
    # if 1000, vector exhausted
    RPCOUNT.DF$norm <- round( (RPCOUNT.DF$norm/min(RPCOUNT.DF$norm[RPCOUNT.DF$norm!=0], na.rm=TRUE)*multiplier), digits=0 )
    # NaN 0/0
    RPCOUNT.DF$norm[is.nan(RPCOUNT.DF$norm)] <- 0L
    
    save(RPCOUNT.DF, file=paste0(out.dir, "/", out.name, "_times", multiplier, 
                                 "_windowDens.RData"))
    
    print("Count data regenerated!")
    
  } else {
    load(file=paste0(out.dir, "/", out.name, "_times", multiplier, 
                     "_windowDens.RData"))
    
  }
  
  file.nme <- load(file=DOMXYZRFEATfile)
  eval(parse(text=paste0(
    "DTA <- ", file.nme, "; rm(", file.nme, "); gc()" 
  )
  ))
  rmax <- max(DTA$radDist, na.rm=TRUE)
  rm(DTA); gc()
  
  # Should cover all the rVal; values will vary between Chrom3D structures
  # but should be the same regardless if Periphery=0/1 or both
  #hist.breaks <- unique(c(seq(from=0, to=rmax, by=hist.breaks), rmax))
  
  # Generate density plot for raw and normalized counts 
  foiNme.v <- unique(RPCOUNT.DF$foi)
  foiNme.v.len <- length(foiNme.v)
  foilst.start <- head(foiNme.v, n=1)
  
  if(regeneratePlotData==TRUE){
    
    print("Calculating density")
    
    # RPCOUNT.DENS will be saved as a list
    RPCOUNT.DENS <- sapply(X=count.v, simplify=FALSE, FUN=function(cnt){
      print(cnt)
      
      toExport <- c("foiNme.v", "RPCOUNT.DF", "cnt", "hist.breaks")
      
      #### PARALLEL EXECUTION #########
      
      dens.mx <- foreach( itr=isplitVector(1:foiNme.v.len, chunks=nCPU), 
                          .combine="cbind", .inorder=TRUE,
                          .export=toExport, 
                          .noexport=ls()[!ls()%in%toExport] 
      ) %op% {
        dens.col <- sapply(X=itr, simplify=FALSE, FUN=function(i){
          foiNme <- foiNme.v[i]
          print(foiNme)
          df1 <- RPCOUNT.DF[RPCOUNT.DF$foi==foiNme,]
          
          # Check if vector limit reached
          cnt.sum <- sum(df1[[cnt]], na.rm=TRUE)
          if(cnt.sum > .Machine$integer.max){
            stop("Length of vector greater than limit.")
          }
          
          # Generate the long vector
          vec <- unlist(mapply(rep, df1$rVal, times=df1[[cnt]]), use.names=FALSE)
          print("Long vector generated.")
          hst <- hist(vec, breaks=hist.breaks, plot=FALSE) 
          
          rm(df1, vec, cnt.sum); gc()
          
          if( !identical(hist.breaks, hst$breaks) ){
            stop("Breaks used not identical to what was set.")
          }
          
          if(foiNme==foilst.start){
            out <- cbind(hst$mids, hst$density)
          } else {
            out <- matrix(hst$density)
          }
          
          rm(hst); gc()
          
          out
        })
        
        dens.mx.chunk <- do.call("cbind", dens.col)
        return(dens.mx.chunk)
      }
      
      ### END OF PARALLEL EXECUTION ###
      
      dimnames(dens.mx) <- list(NULL, c("mid", foiNme.v))
      dens.mx
    })
    
    rm(RPCOUNT.DF, foiNme.v, hist.breaks); gc()
    
    print("Plot data regenerated!")
    
    save(RPCOUNT.DENS,  file=paste0(out.dir, "/", out.name, "_times", multiplier, 
                                    "_windowDensPlot.RData"))
    
  } else {
    load(file=paste0(out.dir, "/", out.name, "_times", multiplier, 
                     "_windowDensPlot.RData"))
  }
  
  file.nme <- load(file=DOMXYZRFEATfile)
  eval(parse(text=paste0(
    "DTA <- ", file.nme, "; rm(", file.nme, "); gc()" 
  )
  ))
  rmax <- max(DTA$radDist, na.rm=TRUE)
  rm(DTA); gc()
    
  for(cnt in count.v){
    df <- melt(as.data.frame(RPCOUNT.DENS[[cnt]]), id="mid")
    ggplot( data=df, aes(x=mid, y=value, group=variable) ) +
      #geom_line( size=1, aes(colour=factor(variable)) ) +
      geom_point(size=3, aes(colour=factor(variable))) + 
      scale_x_continuous(limits=c(0, ceiling(rmax)), breaks=0:ceiling(rmax)
                         ) + 
      scale_colour_manual(values=coul) +
      guides(colour=guide_legend(ncol=1)) +
      labs(title=paste0(out.name, "_", cnt),
           #x=paste0("r±", dr), y="Density",
           x="r", y="Density",
           colour=feat.type
      ) +
      bgr2 
    
    rm(df); gc()
    
    ggsave(file=paste0(out.dir, "/", out.name, "_times", multiplier, 
                       "_", cnt, "_windowDensPlot.pdf"),
           units="in", width=10, height=10)
    print(paste0(cnt, "plot done!"), quote=FALSE)
  }
  
  gc()
  
}
################################################################################
plotFeatDensPerRadWindow <- cmpfun(plotFeatDensPerRadWindow, options=list(suppressUndefined=TRUE))
################################################################################
