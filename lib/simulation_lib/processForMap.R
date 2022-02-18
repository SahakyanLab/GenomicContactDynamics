################################################################################
# Map colour per metric (Cs, Cp, C||, simulation). Function cannot distinguish
# 0s from NAs.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(viridis)
# library(ggplot2)
# library(scales)
### FUNCTION ###################################################################
processForMap <- function(x, metric, is.contProb, is.scaleContactByDist, species.id){
  
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
                     rnm=levels(x))
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value="white")
  
  }  else if(metric=="Cp"){
    
    usedlvl.v <- sort(unique(x), decreasing=FALSE)
    usedlvl.v <- usedlvl.v[!is.na(usedlvl.v)]
    usedlvl.v <- as.character(usedlvl.v[usedlvl.v!=0])

    x <- factor(x=as.character(x), levels=usedlvl.v)
    #coul <- setNames(object=c( "black", viridis::viridis_pal(option="viridis")(length(usedlvl.v)) ),
    #                 nm=levels(x))
    #coul <- setNames(object=c( "black", colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(length(usedlvl.v)) ),
    #                 nm=levels(x))
    coul <- setNames(object=colorRampPalette( rev(brewer.pal(11, "Spectral")) )(length(usedlvl.v)), nm=usedlvl.v)
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=FALSE, na.value="white",
                                      guide=guide_legend(ncol=1), name=bquote(c[p]))
    
  } else if( grepl(x=metric, pattern="CII.disc.", fixed=TRUE) ){
    
    x <- factor(x=as.character(x), levels=c("-1", "0", "1"))
    #coul <- setNames(object=viridis::viridis_pal(option="viridis")(3), 
    #                 nm=c("-1", "0", "1"))
    coul <- setNames(object=brewer.pal(n=11, name="RdYlBu")[c(11, 5, 1)], 
                     nm=c("-1", "0", "1"))
    coul[["0"]] <- "ivory" #"papayawhip"
    #coul[["0"]] <- adjustcolor(col=coul[["0"]], alpha.f=0.2)
    
    leg.lab <- strsplit(x=metric, split=".", fixed=T)[[1]][3]
    leg.lab <- bquote(c['||']~.(leg.lab))
    
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=FALSE, 
                                      na.value="white", name=leg.lab)
    
  #} else if( grepl(x=metric, pattern="CII.cont.", fixed=TRUE) ){
    
  #  x <- as.numeric(x)
  #  col.lst[[1]] <- scale_fill_gradientn(colours=c("#313695", brewer.pal(n=11, name="RdYlBu")[c(9, 6, 4, 1)]), 
  #                                       # scales::squish() for squishing out of bounds values into range.
  #                                       na.value="white", oob=scales::squish,
  #                                       limits=boxplot.stats(x=x, coef=1.5)$stats[c(1,5)]) 
    
  } else if( grepl(x=metric, pattern="Cs.raw|Cs.norm|SIM.|CII.cont", fixed=FALSE) ){
    
    x <- as.numeric(x)
    
    if( grepl(x=metric, pattern="Cs.", fixed=TRUE) ){
      
      leg.lab <- strsplit(x=metric, split=".", fixed=T)[[1]][2]
      leg.lab <- bquote(c['f']~.(leg.lab))
      
    } else if( grepl(x=metric, pattern="CII.", fixed=TRUE) ){
      
      leg.lab <- strsplit(x=metric, split=".", fixed=T)[[1]][3]
      leg.lab <- bquote(c['||']~.(leg.lab))
      
    } else if( grepl(x=metric, pattern="SIM.", fixed=TRUE) ){
      leg.lab <- bquote(c['sim'])
    } else {
      stop("processForMap(): Invalid continuous metric.")
    }
    
    col.v <- c("#313695", brewer.pal(n=11, name="RdYlBu")[c(9, 6, 4, 1)])
    if( grepl(x=metric, pattern="CII.cont.G.", fixed=TRUE) ){
      col.v <- rev(col.v)
    }
    
    mean.val <- mean(x, na.rm=T)
    sd.val <- sd(x, na.rm=T)
    
    if(is.contProb){
      
      #col.lst[[1]] <- scale_fill_gradientn(colours=brewer.pal(n=11, name="RdYlBu")[c(11, 8, 6, 4, 2, 1)], 
      #                                     na.value="white", limits=c(0,1), name=leg.lab) # whole5
      col.lst[[1]] <- scale_fill_gradientn(colours=col.v, na.value="white", 
                                           limits=c(0,1), name=leg.lab) 
    
    } else if( is.scaleContactByDist & grepl(x=metric, pattern="CII.cont", fixed=TRUE) ){
      
      # Option 1, chr17 final
      lim.min <- mean.val - (sd.val / 3)
      lim.max <- mean.val + (2 * sd.val)
      
      # Option 2 (so so for chr17 minmax scaling)
      #lim.min <- mean.val - (sd.val / 5)
      #lim.max <- mean.val + (3 * sd.val)
      
      # NEW, chr1
      #lim.min <- mean.val - (sd.val / 2)
      #lim.max <- mean.val + (2 * sd.val)
      
      # For chr21
      #lim.min <- mean.val - (sd.val / 3)
      #lim.max <- mean.val + (2 * sd.val)
        
      col.lst[[1]] <- scale_fill_gradientn(colours=col.v, na.value="white",
                                           # scales::squish() for squishing out of bounds values into range.
                                           oob=scales::squish, name=leg.lab,
                                           limits=c(lim.min, lim.max)
                                           )
      
    } else if( species.id%in%c("osa") ){
      
      lim.min <- mean.val - (sd.val / 2)
      lim.max <- mean.val + (3 * sd.val)
      
      col.lst[[1]] <- scale_fill_gradientn(colours=col.v, na.value="white",
                                           # scales::squish() for squishing out of bounds values into range.
                                           oob=scales::squish, name=leg.lab,
                                           limits=c(lim.min, lim.max)
      )
      
    } else {
      
      col.lst[[1]] <- scale_fill_gradientn(colours=col.v, na.value="white",
                                           # scales::squish() for squishing out of bounds values into range.
                                           oob=scales::squish, name=leg.lab,
                                           limits=c(boxplot.stats(x=x, coef=1.5)$stats[c(1,5)])
                                           ) 
      
    }
    
  } 
  
  if( !identical(NA.TF, is.na(x)) ){
    stop("processForMap(): Conversion to factor introduced new NAs.")
  }
  
  return(list(value=x, col=col.lst))
  
}
################################################################################

