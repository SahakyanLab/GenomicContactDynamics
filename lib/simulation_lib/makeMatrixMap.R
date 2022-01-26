################################################################################
# Make a matrix map, either symmetrical or square format from dataframe of i, j
# and value. Note that for metrics other than CII, 0 contact values are converted 
# to NAs. processForMap() function is not built to differentiate 0 values from NA 
# values due to invalidity of contact based on contact filtering arguments. 
# Transformation to reduce positive skewness even for zero (was already converted 
# to NA before transformation) and negative values.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(ggnewscale)
# library(cowplot)
### FUNCTION ###################################################################
makeMatrixMap <- function(df.lst = 'list of 1 (upper) or 2(upper and lower) data.frame(x, y, value)', 
                          symmetric = 'T/F', metric.v = 'metric per dataframe in the list', 
                          scalebr.v = 'c(xmin=100, xmax=350, ymin=1, ymax=30) or NULL ',
                          is.contProb = F, is.scaleContactByDist = F, check.dup = 'T/F', plot.title = '', 
                          mark.x = NULL, mark.y = NULL, limits.x = c(NA, NA), limits.y = c(NA, NA),
                          species.id = "e.g. osa"
                          ){
  
  if( !is.null(scalebr.v) ){
    
    plot.title <- paste0( plot.title, "_scaleXbyY:", (scalebr.v["xmax"]-scalebr.v["xmin"]+1),
                          "by", (scalebr.v["ymax"]-scalebr.v["ymin"]+1) )
    
  }
  
  if(check.dup){
    
    dup.TF <- !duplicated(df)
    if( any(dup.TF) ){
      stop("makeMatrixMap(): Duplicates in df", quote=F)
    }
    
  }
  
  metric.v.len <- length(metric.v)
  if( length(df.lst)!=metric.v.len ){
    stop("makeMatrixMap(): Lengths of df.lst and metric.v not equal.")
  }
  
  col.lst <- minmax.lst <- list()
  for(m in 1:metric.v.len){
    
    # Process fill values depending on metric and decide colouring system for map
    lst <- processForMap(x=df.lst[[m]]$value, metric=metric.v[m], is.contProb=is.contProb,
                         is.scaleContactByDist=is.scaleContactByDist, species.id=species.id)
    df.lst[[m]]$value <- lst$value
    col.lst[[m]] <- lst$col
    
    minmax.lst[[m]] <- NA
    # Get min and max metric value given all valid contacts (independent of plot limits set)
    if( !grepl(x=metric.v[m], pattern="CII.disc.|Cp", fixed=F, ignore.case=F) ){
      
      minmax.lst[[m]] <- paste0( min(df.lst[[m]]$value, na.rm=T), ",", 
                                 max(df.lst[[m]]$value, na.rm=T))
      
    } 
    
  } 
  
  #-------------------Parameters of heatmap
  
  for( col in c("i", "j") ){
    
    tmp <- paste(paste0("df.lst[[", 1:metric.v.len, "]]$", col), collapse=", ")
    eval(parse(text=paste0(
      "brk.", col, " <- sort(unique(c(", tmp, ")), decreasing=F)"
    )))
    
  }
  
  len.x <- length(brk.j)
  len.y <- length(brk.i)
  
  if( is.null(limits.x) ){
    limits.x <- range(brk.j, na.rm=F)
  }
  if( is.null(limits.y) ){
    limits.y <- range(brk.i, na.rm=F)
  }
  
  # Define which breaks to show (defines axes labels)
  if( (nrow(df.lst[[1]]))>50 ){
    
    brk.j <- sort(unique( c(limits.x, mark.x, brk.j[c(1, len.x/2 , len.x)]) ))
    brk.i <- sort(unique( c(limits.y, mark.y, brk.i[c(1, len.y/2 , len.y)]) ))
    
  } 
  rm(len.x, len.y)
  
  #-------------------Generate heatmap

  base.lst <- list(
    
    geom_hline(yintercept=mark.y, colour="gray80", size=0.1),
    geom_vline(xintercept=mark.x, colour="gray80", size=0.1),
      scale_x_continuous(breaks=brk.j, limits=limits.x),
      scale_y_continuous(breaks=brk.i, limits=rev(limits.y), trans="reverse" ),
      labs(x=NULL, y=NULL, 
           title=paste0(plot.title, "_MINMAXna.rmTRUE=", paste(unlist(minmax.lst), collapse="_"),
                        "_0valuesForMetricsOtherThanCIIsetToNA")),
    theme(
      plot.title=element_text(size=2),
      axis.text.x=element_text(face="bold", size=5, angle=360, colour="black"),
      axis.text.y=element_text(face="bold", size=5, angle=360, colour="black"),
      axis.ticks=element_line(size=0.1), 
      axis.ticks.length=unit(.15, "cm"),
      legend.text=element_text(size=15),
      legend.title=element_text(size=25),
      panel.grid.major=element_blank(),       
      panel.grid.minor=element_blank(),
      panel.background=element_rect(colour="gray22", size=1, fill=NA),
      aspect.ratio=1 
    )
    
  )
  
  if( !is.null(scalebr.v) & length(scalebr.v)==4 ){
    
    base.lst <- c(base.lst, 
                  geom_rect(colour="black", fill="black", 
                            aes(xmin=scalebr.v["xmin"], xmax=scalebr.v["xmax"], 
                                ymin=scalebr.v["ymin"], ymax=scalebr.v["ymax"]))
                  )
    
  }
  
  test.TF <- symmetric & length(unique(metric.v))==2  
  if(test.TF){
    
    p <- ggplot() + 
      geom_raster(data=df.lst[[2]], aes(x=j, y=i, fill=value)) + 
      col.lst[[2]] + 
      # geoms below will use another color scale
      ggnewscale::new_scale_fill() +
      # Second scale should have the legend that will stay
      geom_raster(data=df.lst[[1]], aes(x=j, y=i, fill=value)) +
      col.lst[[1]] +
      geom_abline(slope=-1, intercept=0, colour="gray20") +
      base.lst
    
    pDown <- ggplot() + 
      geom_raster(data=df.lst[[1]], aes(x=j, y=i, fill=value)) + 
      col.lst[[1]] +
      base.lst
    
    legDown <- cowplot::get_legend(pDown)
    rm(pDown)
      
  } else {
    
    df <- do.call("rbind.data.frame", df.lst)
    rownames(df) <- NULL
    rm(df.lst)
    gc()
    
    p <- ggplot() +
      geom_raster(data=df, aes(x=j, y=i, fill=value)) + 
      col.lst[[1]] + 
      base.lst
       
  }
  
  if(test.TF){
    
    p <- cowplot::plot_grid(p + 
                              guides(fill="none") + 
                              theme(legend.position="left"), 
                            legDown,
                            ncol=2, rel_widths=c(.85, .10))
    
  } 
  
  return(p)
  
}
################################################################################
# rm(list=ls()); gc()

#df <- expand.grid(i=1:3, j=1:3)
#df <- df[df$j>df$i,]
#df$value <- 7:9
#df <- rbind(df, cbind(i=df$j, j=df$i, value=1:3))
#df$group <- c(rep("a",times=3), rep("b", times=3))

#p0 <- ggplot(data=df, aes(x=j, y=i)) +
#  scale_y_continuous(trans="reverse") + 
#  geom_raster(data=df[4:6,], aes(fill=value)) + 
#  scale_fill_viridis_c(option="D") +
#  # geoms below will use another color scale
#  ggnewscale::new_scale_fill() +
#  geom_raster(data=df[1:3,], aes(fill=value)) +
#  # Second scale should have the legend that will stay
#  scale_fill_viridis_c(option="A") # pinkish

# Should have the legend that will stay
#pDown <- ggplot(data=df[1:3,], aes(x=j, y=i)) +
#  geom_raster(aes(fill=value)) + 
#  scale_fill_viridis_c(option="A")
#legDown <- get_legend(pDown)

#cowplot::plot_grid(p0 + 
#          guides(fill = "none") + 
#          theme(legend.position="left"), 
#          legDown, 
#          ncol = 2, rel_widths = c(.85, .15))


#https://stackoverflow.com/questions/68369581/how-to-place-legends-at-different-sides-of-plot-bottom-and-right-side-with-ggp/68369737#68369737
