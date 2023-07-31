################################################################################
# Map colour per metric (Cs, Cp, C||, simulation). Function cannot distinguish
# 0s from NAs.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(ggplot2)
# library(scales)
### FUNCTION ###################################################################
processForMap <- function(x, metric, is.contProb, is.scaleContactByDist, species.id,
                          na.colour="white"){
  
  # For checking whether conversion to factor didn't introduce NAs
  NA.TF <- is.na(x)
  
  col.lst <- list()
  if(metric=="Cs.raw-categorised"){
    
    x <- cut(
      
      x=as.numeric(x),
      # Breaks determined manually based on Cf distribution
      breaks=c(-Inf, 0, 1, 2, 3, 4, 5, 10, Inf), 
      include.lowest=FALSE, right=TRUE,
      labels=c("0", "1", "2", "3", "4", "5", "(5,10]", ">10")
      
    )
    
    coul <- setNames(object=c("black", brewer.pal(name="RdYlBu", n=11)[c(11:8,5,3,1)]), nm=levels(x))
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value=na.colour)
  
  }  else if(metric=="Cp"){
    
    ## Generic approach, regardless of length of Cp dimension, used for the atlas and manuscript figures
    #usedlvl.v <- sort(as.numeric(unique(x)), decreasing=FALSE)
    #usedlvl.v <- usedlvl.v[ !is.na(usedlvl.v) & usedlvl.v != 0 ] 
    #usedlvl.v <- as.character(usedlvl.v) 
    
    # For Cp dimension 1 to 21
    usedlvl.v <- as.character(1:21)
    
    coul <- colorRampPalette( rev(brewer.pal(11, "Spectral")) )(length(usedlvl.v))
    if(0 %in% x){ # Contacts not in tissue assigned a Cp=0
      usedlvl.v <- c("0", usedlvl.v)
      coul <- c("black", coul)
    }
    
    coul <- setNames(object=coul, nm=usedlvl.v)
    x <- factor(x=as.character(x), levels=usedlvl.v)
    
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, na.value=na.colour,
                                      guide=guide_legend(ncol=1), name=bquote(c[p]))
    
  } else if( grepl(x=metric, pattern="CII.disc.", fixed=TRUE) ){
    
    x <- factor(x=as.character(x), levels=c("-1", "0", "1"))
    coul <- setNames(object=brewer.pal(n=11, name="RdYlBu")[c(11, 5, 1)], nm=c("-1", "0", "1"))
    coul[["0"]] <- "ivory" 
    
    leg.lab <- strsplit(x=metric, split=".", fixed=T)[[1]][3]
    leg.lab <- bquote(c['||']~.(leg.lab))
    
    col.lst[[1]] <- scale_fill_manual(values=coul[levels(x)], na.translate=TRUE, 
                                      na.value=na.colour, name=leg.lab)
    
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
     
      lim.min <- 0
      lim.max <- 1
    
    } else if( is.scaleContactByDist & grepl(x=metric, pattern="CII.cont", fixed=TRUE) ){
      
      # Option 1, chr1 and chr17 first choice
      #lim.min <- mean.val - (sd.val / 3)
      #lim.max <- mean.val + (2 * sd.val)
      
      # Try
      lim.min <- mean.val - (sd.val / 4)
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
        
    } else if( species.id%in%c("osa") ){
      
      lim.min <- mean.val - (sd.val / 2)
      lim.max <- mean.val + (3 * sd.val)
      
    } else {
      
      bp.stat <- boxplot.stats(x=x, coef=1.5)$stats
      lim.min <- bp.stat[1]
      lim.max <- bp.stat[5]
      
    }
    
    col.lst[[1]] <- scale_fill_gradientn(colours=col.v, na.value=na.colour, name=leg.lab,
                                         # scales::squish() for squishing out of bounds values into range
                                         oob=scales::squish, limits=c(lim.min, lim.max))
    
  } 
  
  if( !identical(NA.TF, is.na(x)) ){
    stop("processForMap(): Conversion to factor introduced new NAs.")
  }
  
  return(list(value=x, col=col.lst))
  
}
################################################################################

