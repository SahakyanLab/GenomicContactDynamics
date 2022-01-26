################################################################################
# Map colour per metric (Cs, Cp, C||, simulation). Function can handle 0s and
# colours them black. Colour legend will show all expected factors/categories.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(viridis)
# library(ggplot2)
# library(scales)
# library(stats)
### FUNCTION ###################################################################
processForMap <- function(x, metric, contProb){
  
  # For checking whether conversion to factor didn't introduce NAs
  NA.TF <- is.na(x)
  
  col.lst <- list()
  if(metric=="Cs.raw-categorised"){
    
    x <- cut(
      
      x=as.numeric(x),
      # Based on Cs values of long-range contacts
      breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 10, Inf), 
      include.lowest=FALSE, right=TRUE,
      labels=c("0", "1", "2", "3", "4", "5", "(5,10]", ">10")
      
    )
    #x <- factor(x=as.character(x), levels=c("0", "1", "2", "3", "4", "5", "(5,10]", ">10"))
    coul <- setNames(object=c("black", brewer.pal(name="RdYlBu", n=11)[c(11:8,5,3,1)]),
                     nm=levels(x))
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
  
  }  else if(metric=="Cp"){
    
    usedlvl.v <- sort(unique(x), decreasing=FALSE)
    usedlvl.v <- usedlvl.v[!is.na(usedlvl.v)]
    usedlvl.v <- usedlvl.v[usedlvl.v!=0]

    x <- factor( x=as.character(x),  levels=as.character(c(0, usedlvl.v)) ) 
    coul <- setNames(object=c( "black", viridis::viridis_pal(option="viridis")(length(usedlvl.v)) ),
                     nm=levels(x))
    #coul <- setNames(object=c( "black", colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(length(usedlvl.v)) ),
    #                 nm=levels(x))
    #coul <- setNames(object=c( "black", colorRampPalette( rev(brewer.pal(11, "Spectral")) )(length(usedlvl.v)) ),
    #                 nm=levels(x))
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
    
  } else if( grepl(x=metric, pattern="CII.disc.", fixed=TRUE) ){
    
    x <- factor( as.character(x),  levels=c("-1", "0", "1") )
    coul <- setNames(object=viridis::viridis_pal(option="viridis")(3), 
                     nm=c("-1", "0", "1"))
    #coul <- setNames(object=brewer.pal(n=11, name="RdYlBu")[c(11, 5, 1)], 
    #                 nm=c("-1", "0", "1"))
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
  
  } else if( grepl(x=metric, pattern="Cs.raw|Cs.norm|SIM.|CII.cont.") ){
    
    x <- as.numeric(x)
    med.x <- median(x=x, na.rm=TRUE)
    onemad.x <- mad(x=x, center=med.x, na.rm=TRUE)
    twomad.x <- 2*onemad.x
    
    if(contProb){
      
      col.lst[[1]] <- scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(11, 8, 6, 4, 2, 1)], 
                                           na.value="white", limits=c(0,1)) # whole5
      
    } else {
      
      #col.lst[[1]] <- scale_fill_gradientn(colours=viridis::viridis_pal(option="viridis")(5), 
      #                                     na.value="white", oob=scales::squish, 
      #                                     limits=boxplot.stats(x=x, coef=1.5)$stats[c(2,4)])
      col.lst[[1]] <- scale_fill_gradientn(colours=viridis::viridis_pal(option="viridis")(5), 
                                           na.value="white", oob=scales::squish, 
                                           limits=c(med.x-onemad.x, med.x+twomad.x))
      #col.lst[[1]] <- scale_fill_gradientn(colours=c("#313695", brewer.pal(n=11, name="RdYlBu")[c(9, 6, 4, 1)]), 
      #                                     # scales::squish() for squishing out of bounds values into range.
      #                                     na.value="white", oob=scales::squish,
      #                                     limits=boxplot.stats(x=x, coef=1.5)$stats[c(2,4)]) 
      
      # Not color-blind-friendly
      #col.lst[[1]] <- scale_fill_gradientn(colours=c("#313695", brewer.pal(n=11, name="RdYlGn")[c(9, 6, 4, 1)]), 
      #                                     na.value="white", oob=scales::squish, #breaks=c(0.1, 0.2, 0.3, 0.4),
      #                                     limits=boxplot.stats(x=x, coef=1.5)$stats[c(2,4)]) 
      
    }
    
  }
  
  if( !identical(NA.TF, is.na(x)) ){
    stop("processForMap(): Conversion to factor introduced new NAs.")
  }
  
  return(list(value=x, col=col.lst))
  
}
################################################################################

