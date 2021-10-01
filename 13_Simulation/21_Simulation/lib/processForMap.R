################################################################################
# Map colour per metric (Cs, Cp, C||, simulation). Function can handle 0s and
# colours them black. Colour legend will show all expected factors/categories.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(ggplot2)
### FUNCTION ###################################################################
processForMap <- function(x=x, metric=metric){
  
  # For checking whether conversion to factor didn't introduce NAs
  NA.TF <- is.na(x)
  
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
    col.lst <- list(
      scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
    )
    
  }  else if(metric=="Cp"){
    
    usedlvl.v <- sort(unique(x), decreasing=FALSE)
    usedlvl.v <- usedlvl.v[!is.na(usedlvl.v)]
    usedlvl.v <- usedlvl.v[usedlvl.v!=0]

    x <- factor( x=as.character(x),  levels=as.character(c(0, usedlvl.v)) ) 
    coul <- setNames(object=c( "black", colorRampPalette( rev(brewer.pal(11, "Spectral")) )(length(usedlvl.v)) ),
                     nm=levels(x))
    col.lst <- list(
      scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
    )
    
  } else if( grepl(x=metric, pattern="CII.disc.", fixed=TRUE) ){
    
    x <- factor( as.character(x),  levels=c("-1", "0", "1") )
    coul <- setNames(object=c("#4292C6", "#FDC776","#9E0142"), 
                     nm=c("-1", "0", "1"))
    col.lst <- list(
      scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
    )
    
  } else if( grepl(x=metric, pattern="Cs.raw|Cs.norm|SIM.|CII.cont.") ){
    
    x <- as.numeric(x)
    col.lst <- list(
      
      #scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(n=5, name="RdYlBu"))(5)), 
      #                     na.value="white")
      #scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n=11, name="RdYlBu")[c(10,5,3,1)])(5), 
      #                     na.value="white")
      #scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(10, 5:3, 1)], 
      #                     na.value="white") # whole1
      #scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(11, 9, 5:3, 1)], 
      #                     na.value="white") # whole2
      #scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(11, 9, 6, 5, 2, 1)], 
      #                     na.value="white")  # whole3
      #scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(11, 6, 4, 2, 1)], 
      #                     na.value="white") # whole4
      scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(11, 8, 6, 4, 2, 1)], 
                           na.value="white") # whole5
      
    )
    
  }
  
  if( !identical(NA.TF, is.na(x)) ){
    stop("processForMap(): Conversion to factor introduced new NAs.")
  }
  
  return(list(value=x, col=col.lst))
  
}
################################################################################

