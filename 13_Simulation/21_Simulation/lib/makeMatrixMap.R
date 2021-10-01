################################################################################
# Make a matrix map, either symmetrical or square format from dataframe of i, j
# and value. Note that for metrics other than CII, 0 contact values are converted 
# to NAs but the processForMap() function is built to deal with 0 values by
# making them black on the map to differentiate them from NA values due to 
# invalidity of contact based on contact filtering arguments. Transformation to 
# reduce positive skewness even for zero (was already converted to NA before 
# transformation) and negative values.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(ggplot2)
# source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makeMatrixMap <- function(df = data.frame(x, y, value), # lower/upper triangle only
                          metric = NULL, format = "symmetric", # symmetric | square
                          scalebr.v = c(xmin=100, xmax=350, ymin=1, ymax=30),
                          check.dup = TRUE, plot.title = NULL,
                          mark.x = NULL, mark.y = NULL, 
                          limits.x = c(NA, NA), limits.y = c(NA, NA)
                          ){
  
  if( !is.null(scalebr.v) ){
    
    plot.title <- paste0(plot.title, "_scaleXbyY:", (scalebr.v["xmax"]-scalebr.v["xmin"]),
                         "by", (scalebr.v["ymax"]-scalebr.v["ymin"])
                         )
    
  }
  
  if(check.dup){
    
    notdup.TF <- !duplicated(df)
    if( !all(notdup.TF) ){
      print("makeMatrixMap(): Duplicates found and removed.", quote=FALSE)
      df <- df[notdup.TF,]
    }
    
  }
  
  #-------------------Prepare dataframe for making heatmap
  
  colnames(df) <- c("i", "j", "value")
  
  tri.TF <- all(df$i<df$j) | all(df$j<df$i) | all(df$i==df$j)
  if( all(tri.TF==FALSE) ){
    stop("makeMatrixMap(): df not from one triangle.")
  }
  rm(tri.TF)
  
  print(paste0("makeMatrixMap(): Generating ", format, " map..."), quote=FALSE)
  
  if(format=="symmetric"){
    
    df <- rbind.data.frame(df, cbind.data.frame(i=df$j, j=df$i, value=df$value), 
                           stringsAsFactors=FALSE)
    
  } else if(format=="square"){
    
  } else {
    stop("makeMatrixMap(): Invalid format argument.")
  }
  
  # For metrics other than CII, 0 values converted to NAs
  if( !grepl(x=metric, pattern="CII.cont.|CII.disc.") ){
    df[!is.na(df$value) & df$value==0,"value"] <- NA
  }
  
  # Transformation to reduce positive skewness even for zero (was
  # already converted to NA before transformation) and negative values
  if( grepl(x=metric, pattern="Cs.raw|Cs.norm|CII.cont.") ){
    
    nonNA.TF <- !is.na(df$value)
    df$value[nonNA.TF] <- log(df$value[nonNA.TF] + 1)
    #df$value[nonNA.TF] <- (df$value[nonNA.TF])^(1/3)
    rm(nonNA.TF)
    
  }
  
  #-------------------Parameters of heatmap
  
  # Axes break values, upper triangle perspective
  brk.x <- sort(unique(df$j))
  brk.y <- sort(unique(df$i))
  len.x <- length(brk.x)
  len.y <- length(brk.y)
  
  if( is.null(limits.x) ){
    limits.x <- range(df$j, na.rm=FALSE)
  }
  if( is.null(limits.y) ){
    limits.y <- range(df$i, na.rm=FALSE)
  }
  
  # Define which breaks to show (defines axes labels)
  if( (nrow(df)/2)>50 ){
    
    brk.x <- sort(unique( c(limits.x, mark.x, brk.x[c(1, len.x/2 , len.x)]) ))
    brk.y <- sort(unique( c(limits.y, mark.y, brk.y[c(1, len.y/2 , len.y)]) ))
    
  }
  rm(len.x, len.y)
  
  # Get min and max metric value given all valid contacts (independent of plot limits set)
  if( !grepl(x=metric, pattern="CII.disc.") ){
    minmax <- paste( min(df$value, na.rm=T), ",", max(df$value, na.rm=T), collapse="," )
  }
  
  # Process fill values depending on metric and decide colouring system for map
  lst <- processForMap(x=df$value, metric=metric)
  df$value <- lst$value
  lst$value <- NULL
  
  #-------------------Generate heatmap
  
  p <- ggplot(data=df, aes(x=j, y=i)) + # Upper triangle perspective
    geom_raster(aes(fill=value)) +
    { if(format=="symmetric") geom_abline(slope=-1, intercept=0, colour="gray20") } +
    geom_hline(yintercept=mark.y, colour="gray80", size=0.1) +
    geom_vline(xintercept=mark.x, colour="gray80", size=0.1) +
    scale_x_continuous(breaks=brk.x, limits=limits.x) + 
    scale_y_continuous(breaks=brk.y, limits=rev(limits.y), trans="reverse" ) + 
    lst$col +
    labs(x=NULL, y=NULL, 
         title=paste0(plot.title, "MINMAXna.rmTRUE=", minmax,
                      "_0valuesForMetricsOtherThanCIIsetToNA")) +
    bgr2 + 
    theme(
      plot.title=element_text(size=2),
      axis.text.x=element_text(face="bold", size=1.5, angle=90, colour="black"),
      axis.text.y=element_text(face="bold", size=1.5, angle=360, colour="black"),
      axis.ticks=element_line(size=0.1), 
      axis.ticks.length=unit(.15, "cm"),
      legend.text=element_text(size=5, face="bold"),
      legend.title=element_text(size=10, face="bold"),
      legend.position="bottom",
      legend.direction="horizontal", 
      legend.box="vertical"
    )  
  
  if( !is.numeric(df$value) ){
    p <-  p + guides(fill=guide_legend(nrow=1, label.position="bottom")) 
  }
  
  if( !is.null(scalebr.v) & length(scalebr.v)==4 ){
    
    p <- p +
      geom_rect(colour="black", fill="black", 
                aes(xmin=scalebr.v["xmin"], xmax=scalebr.v["xmax"], 
                ymin=scalebr.v["ymin"], ymax=scalebr.v["ymax"]))
  
  }
  
  return(p)
  
}
################################################################################
# rm(list=ls()); gc()
