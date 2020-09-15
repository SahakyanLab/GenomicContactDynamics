################################################################################
# Map colour per metric (Cs, Cp, C||, simulation)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# library(ggplot2)
### FUNCTION ###################################################################
processForMap <- function(x=x, metric=metric){
  
  #-------------------raw Cs
  if(metric=="Cs.raw"){
    
    x <- cut(
      x=as.numeric(x), 
      breaks=c(-Inf, 1, 2, 3, 4, 5, 10, Inf), 
      include.lowest=FALSE, right=TRUE,
      labels=c("1", "2", "3", "4", "5", "(5,10]", ">10")
    )
    x <- factor(x=x, levels=c("1", "2", "3", "4", "5", "(5,10]", ">10"))
    
    # 5 blues from brewer.pal(name="Blues") plus yellow and red
    coul <- c(brewer.pal(name="Blues", n=9)[c(2,3,5,7,9)], "#FDC776","#9E0142")
    names(coul) <- c("1", "2", "3", "4", "5", "(5,10]", ">10")
    coul <- coul[names(coul)%in%levels(x)]
    col.lst <- list(
      scale_fill_manual(values=coul, na.translate=TRUE, na.value="white")
    )
    
  }  else if(metric=="Cp"){
    
    x <- factor( as.character(x),  levels=as.character(sort(unique(x))) )
    coul <- colorRampPalette( rev(brewer.pal(11, "Spectral")) )(21)
    col.lst <- list(
      scale_fill_manual(values=coul, na.translate=TRUE, na.value="white")
    )
    
  } else if( grepl(x=metric, pattern="CII.disc.", fixed=TRUE) ){
    
    x <- factor( as.character(x),  levels=as.character(sort(unique(x))) )
    coul <- c("#4292C6", "#FDC776","#9E0142")
    names(coul) <- c("-1", "0", "1")
    col.lst <- list(
      scale_fill_manual(values=coul, na.translate=TRUE, na.value="white")
    )
    
  } else if( grepl(x=metric, pattern="Cs.norm|SIM.|CII.cont.") ){
    
    x <- as.numeric(x)
    coul <- brewer.pal(name="Reds", n=9)[c(5,9)]
    col.lst <- list(
      scale_fill_gradient(low=coul[1], high=coul[2], na.value="white")
    )
    
  }
  
  return(list(value=x, col=col.lst))
  
}
################################################################################
