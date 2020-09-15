################################################################################
# Make a matrix map, either symmetrical or square format from dataframe of i, j
# and value.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(ggplot2)
# source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makeMatrixMap <- function(df = data.frame(x, y, value), # lower/upper triangle only
                          metric = NULL,
                          format = "symmetric", # symmetric | square
                          check.dup = TRUE,
                          scalebr.v = c(xmin=100, xmax=350, ymin=1, ymax=30),
                          plot.title = NULL
){
  
  if(check.dup){
    notdup.TF <- !duplicated(df)
    if( !all(notdup.TF) ){
      print("Duplicates found and removed.", quote=FALSE)
      df <- df[notdup.TF,]
    }
  }
  #-------------------Prepare dataframe for making heatmap
  colnames(df) <- c("i", "j", "value")
  
  tri.TF <- all(df$i<df$j) | all(df$j<df$i) | all(df$i==df$j)
  if( all(tri.TF==FALSE) ){
    stop("df not from one triangle.")
  }; rm(tri.TF)
  
  print(paste0("Generating ", format, " map..."), quote=FALSE)
  if(format=="symmetric"){
    df <- rbind.data.frame(df, cbind.data.frame(i=df$j, j=df$i, value=df$value), 
                           stringsAsFactors=FALSE)
  }
  
  if( !grepl(x=metric, pattern="CII.cont.|CII.disc.") ){
    df[!is.na(df$value) & df$value==0,"value"] <- NA
  }
  #-------------------Parameters of heatmap
  brk.x <- sort(unique(df$i))
  brk.y <- sort(unique(df$j))
  len <- length(brk.x)
  if( (nrow(df)/2)>100 ){
    seq.v <- ceiling(seq.int(from=1, to=len, length.out=4))
    brk.x <- brk.x[seq.v]
    brk.y <- brk.y[seq.v]
    rm(seq.v)
  }; rm(len)
  
  # Process fill values depending on metric and decide colouring system for map
  lst <- processForMap(x=df$value, metric=metric)
  df$value <- lst$value; lst$value <- NULL
  #-------------------Generate heatmap
  p <- ggplot(data=df, aes(x=i, y=j)) + 
    geom_raster(aes(fill=value)) +
    scale_x_continuous(breaks=brk.x ) + 
    scale_y_continuous(breaks=brk.y, trans="reverse" ) + 
    lst$col +
    labs(x=NULL, y=NULL, title=plot.title) +
    bgr2 + 
    theme(
      plot.title=element_text(size=2),
      axis.text.x=element_text(face="bold", size=10, angle=360, colour="black"),
      axis.text.y=element_text(face="bold", size=10, angle=360, colour="black"),
      axis.ticks.y=element_blank(),
      axis.ticks.x=element_blank(),
      legend.text=element_text(size=10, face="bold"),
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
