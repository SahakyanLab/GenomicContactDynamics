################################################################################
# Make Hi-Cs and Hi-Cp map of regions of interest
# MX should be ibin-jbin-Cs-Cp-CII
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(ggplot2)
# library(RColorBrewer)
#source(paste0(lib, "/GG_bgr.R"))
# source(paste0(lib, "/multiplot.R"))
# source(paste0(wk.dir, "/lib/hmplot.R"))
################################################################################
visualiseBinRegions <- function(
  out.dir = "/dir",
  out.name = paste0(chr, "_", gcb, "_", ct, "_", region), # KRAS_scaled
  MX = PERSIST.MX$hits,
  metric.v = c("Cs", "Cp", "CII"),
  bins.x = 1203:1205,
  bins.y = 1146:1150,
  format = "square", # "symmetric" | "square"
  scalebr.v = c(xmin=1, xmax=50, ymin=1, ymax=30),
  outtype = outtype
){
  
  metric.v.len <- length(metric.v)
  if( ncol(MX)==(metric.v.len+2) ){
    colnames(MX) <- c("i", "j", metric.v)
    rownames(MX) <- NULL
  } else {
    stop("Incorrect MX.")
  }
  
  # Filter contacts based on specificed bins
  bins <- sort(as.numeric(c(bins.x, bins.y)))
  if( !is.null(bins.x) & !is.null(bins.y) ){
    # Get desired contacts based on bins.x and bins.y
    MX <- MX[ MX[,"i"]%in%bins & MX[,"j"]%in%bins , ]
  } else if( (!is.null(bins.x) & is.null(bins.y)) | (is.null(bins.x) & !is.null(bins.y)) ){
    MX <- MX[ MX[,"i"]%in%bins | MX[,"j"]%in%bins , ]
  } 

  if(nrow(MX)==0){
    stop("No contacts formed by given bins.")
  }
  
  if(format=="symmetric"){
    # Check if MX is the upper triangle; needed for symmetric to work
    if( any(MX[,"j"]-MX[,"i"] <= 0) ){
      stop("MX not upper triangle.")
    }
    diag.mx <- unique(c(MX[,"i"], MX[,"j"]))
    diag.mx <- cbind(i=diag.mx, j=diag.mx,
                     matrix( data=-2, ncol=metric.v.len, nrow=length(diag.mx), 
                             dimnames=list(NULL, metric.v) )
                     )
    lower.tri.mx <- cbind(MX[,c("j", "i")], MX[,metric.v])
    colnames(lower.tri.mx) <- c("i", "j", metric.v)
    MX <- rbind(MX, lower.tri.mx, diag.mx)
    
    rm(diag.mx, lower.tri.mx); gc()
  }
  
  if( format=="square" & !is.null(bins.x) ){
    # Rearrange so that all bins.x will be in i column
    log <- MX[,"j"]%in%bins.x
    if(sum(log) > 0){
      mx <- cbind(i=MX[log,"j"], j=MX[log, "i"])
      MX[log,"j"] <- mx[,"i"]
      MX[log,"j" ] <- mx[,"j"]
      rm(mx); gc()
    }
    rm(log)
  }

  # Plot
  MX <- as.data.frame(MX)
  MX$i <- factor( 
    x=MX$i, levels=unique( sort(as.numeric(MX$i)) )
    )
  MX$j <- factor( 
    x=MX$j, levels=unique( rev(sort(as.numeric(MX$j))) )
    )
  if("Cs"%in%metric.v){
    MX$Cs <- cut(
      x=as.numeric(MX$Cs), 
      breaks=c(-Inf, -2, 0, 1, 2, 3, 4, 5, 10, Inf), 
      include.lowest=FALSE, right=TRUE,
      labels=c("-2", "0", "1", "2", "3", "4", "5", "(5,10]", ">10")
      #labels=c("(0,1]", "(1,2]", "(2,3]", "(3,4]", 
      #          "(4,5]", "(5,10]", ">10")
    )
    MX$Cs <- factor( x=MX$Cs, levels=c("1", "2", "3", "4", "5", "(5,10]", ">10", "0", "-2") )
  }
  
  if("Cp"%in%metric.v){
    MX$Cp <- factor(x=MX$Cp, levels=c(1:21, -2)) 
  }
  
  if("CII"%in%metric.v){
    MX$CII <- factor( x=MX$CII, levels=c(-1, 0, 1, -2) )
  }
  MX$h <- 1
  
  hmplot(df=MX, out.name=out.name, scalebr.v=scalebr.v, outtype=outtype)
  
}